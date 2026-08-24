"""Pipeline specification utilities for Fiber processing workflows.

:Authors:
    - Sergio Diaz-Meza
    - Jonas Pätzel

:Contributors:
    - Christopher Wollin

"""

import pandas as pd
try:
    from fobench.core.fiber import Fiber
except ModuleNotFoundError:
    from ...fobench.core.fiber import Fiber

def pipeline_check(pipeline: list, allow_callables: bool = True) -> list[dict]:
    """Methods to validate a pipeline specification for parallel computing.
    Crucial to avoid fundamental errors before long parallel executions.

    Parameters
    ----------
    pipeline : list
        Sequence of pipeline steps. Each step can be:

        - ``dict`` with keys:
          ``op``, ``kwargs``, ``name``, ``capture``, ``attribute``,
          ``store_as``.
        - ``str`` shorthand for ``{"op": <step>}``.
        - callable shorthand for ``{"op": <callable>}`` when
          ``allow_callables=True``.

    allow_callables : bool, optional
        If ``True``, callables are accepted as step operators.

    Returns
    -------
    list[dict]
        Normalized pipeline with explicit keys:
        ``name``, ``op``, ``kwargs``, ``capture``, ``attribute``,
        ``store_as``.

    Raises
    ------
    TypeError
        ``pipeline`` or one of the steps has invalid type.
    ValueError
        Step definition is incomplete or inconsistent.

    """

    # Basic input validation before looping over steps.
    if not isinstance(pipeline, (list, tuple)) or (len(pipeline) == 0):
        raise TypeError('"Pipeline" must be a list or tuple of steps.')

    allowed_keys = {"name", "op", "kwargs", "capture", "attribute", "store_as"}
    allowed_capture = {"none", "return", "attribute"}

    normalized = []
    used_store_keys = set()

    # Loop over each step and normalize it to a common dictionary format.
    for i, raw_step in enumerate(pipeline):

        # Allow short definitions to keep user input simple.
        if isinstance(raw_step, str):
            step = {"op": raw_step}
        elif callable(raw_step):
            if not allow_callables:
                raise TypeError(f"Step {i}: callables are disabled for this validator.")
            step = {"op": raw_step}
        elif isinstance(raw_step, dict):
            step = dict(raw_step)
        else:
            raise TypeError(
                f"Step {i}: each step must be dict, str, or callable. Received: {type(raw_step).__name__}."
            )

        # Check unsupported keys early to prevent silent typos.
        unknown = set(step.keys()) - allowed_keys
        if unknown:
            raise ValueError(f"Step {i}: unknown key(s): {sorted(unknown)}.")

        if "op" not in step:
            raise ValueError(f'Step {i}: missing required key "op".')

        op = step["op"]
        if isinstance(op, str):
            if op.strip() == "":
                raise ValueError(f'Step {i}: "op" string can not be empty.')
        elif callable(op):
            if not allow_callables:
                raise TypeError(f"Step {i}: callable operators are disabled.")
        else:
            raise TypeError(
                f'Step {i}: "op" must be str or callable. Received: {type(op).__name__}.'
            )

        # Ensure kwargs always exists as dictionary.
        kwargs = step.get("kwargs", {})
        if kwargs is None:
            kwargs = {}
        if not isinstance(kwargs, dict):
            raise TypeError(f'Step {i}: "kwargs" must be a dict.')

        # Capture mode defines what is saved from each step.
        capture = step.get("capture", "none")
        if capture not in allowed_capture:
            raise ValueError(
                f'Step {i}: invalid "capture" value "{capture}". Use one of {sorted(allowed_capture)}.'
            )

        # Attribute is required only when capture mode points to Fiber attributes.
        attribute = step.get("attribute", None)
        if capture == "attribute":
            if not isinstance(attribute, str) or attribute.strip() == "":
                raise ValueError(
                    f'Step {i}: "capture=attribute" requires non-empty string "attribute".'
                )
        else:
            attribute = None

        # Optional human-readable name.
        name = step.get("name", None)
        if name is not None:
            if not isinstance(name, str) or name.strip() == "":
                raise TypeError(f'Step {i}: "name" must be a non-empty string when provided.')

        # Output key handling for captured results.
        store_as = step.get("store_as", None)
        if capture == "none":
            store_as = None
        else:
            if store_as is None:
                if name is not None:
                    store_as = name
                elif isinstance(op, str):
                    store_as = op
                else:
                    store_as = f"step_{i}"

            if not isinstance(store_as, str) or store_as.strip() == "":
                raise TypeError(f'Step {i}: "store_as" must be a non-empty string.')

            if store_as in used_store_keys:
                raise ValueError(f'Step {i}: duplicate "store_as" key "{store_as}".')
            used_store_keys.add(store_as)

        # Store normalized step so execution code can stay simple.
        normalized.append({
            "name": name,
            "op": op,
            "kwargs": kwargs,
            "capture": capture,
            "attribute": attribute,
            "store_as": store_as,
        })

    return normalized


def run_file_pipeline(task_row, pipeline : list, fiber_kwargs : dict = None,
    validate_pipeline : bool = False):
    """Run a pipeline over one file and return a compact result.

    Parameters
    ----------
    task_row : dict | pandas.Series
        Row-like object containing at least ``"file"``.
        Optional keys such as ``"file_index"`` and ``"window_id"``
        are propagated to output.
    pipeline : list
        Pipeline specification (raw or normalized).
    fiber_kwargs : dict, optional
        Keyword arguments for ``Fiber`` initialization.
    validate_pipeline : bool, optional
        If ``True``, validates pipeline before running.

    Returns
    -------
    dict
        Dictionary containing ``status``, ``error``, file identifiers,
        and ``outputs``.

    """

    # Default empty Fiber kwargs if not provided.
    # For single-file executions, Fiber keeps its own default show_progress=True.
    fiber_kwargs = {} if fiber_kwargs is None else dict(fiber_kwargs)

    # option for accepting both dict rows and pandas Series rows.
    if isinstance(task_row, pd.Series):
        row = task_row.to_dict()
    elif isinstance(task_row, dict):
        row = dict(task_row)
    else:
        raise TypeError('"task_row" must be dict or pandas.Series.')

    if "file" not in row:
        raise KeyError('"task_row" must contain "file".') # crash if this doesn'T work.

    # Validation is optional; for parallel runs this should be done once outside.
    steps = pipeline_check(pipeline, allow_callables=True) if validate_pipeline else pipeline

    # Standard output structure for downstream aggregation.
    result = {
        "status": "ok",
        "error": None,
        "file": row.get("file"),
        "file_index": row.get("file_index"),
        "window_id": row.get("window_id"),
        "outputs": {},
    }

    # Context dictionary can be used by custom callables to share state.
    context = {}
    fiber = None

    try:

        # Initialize Fiber object for the current file.
        fiber = Fiber(row["file"], **fiber_kwargs)

        # Apply pipeline in sequence.
        for i, step in enumerate(steps):

            # Support shorthand pipeline entries also at runtime.
            if isinstance(step, str):
                step = {"op": step}
            elif callable(step):
                step = {"op": step}

            # Read normalized step fields.
            op = step["op"]
            kwargs = {} if step.get("kwargs", None) is None else step.get("kwargs", {})
            capture = step.get("capture", "none")
            store_as = step.get("store_as", None)
            attribute = step.get("attribute", None)

            # Execute Fiber method or user callable.
            if isinstance(op, str):
                if not hasattr(fiber, op):
                    raise AttributeError(f'Fiber has no method "{op}" at step {i}.')
                out = getattr(fiber, op)(**kwargs)
            elif callable(op):
                out = op(fiber=fiber, task=row, context=context, **kwargs)
            else:
                raise TypeError(f'Step {i}: "op" must be str or callable.')

            # Save output depending on selected capture mode.
            if capture == "return":
                key = store_as if store_as is not None else f"step_{i}"
                result["outputs"][key] = out
            elif capture == "attribute":
                if not isinstance(attribute, str) or attribute.strip() == "":
                    raise ValueError(f'Step {i}: "attribute" must be set for capture="attribute".')
                key = store_as if store_as is not None else attribute
                result["outputs"][key] = getattr(fiber, attribute)
            elif capture == "none":
                pass
            else:
                raise ValueError(f'Step {i}: invalid capture mode "{capture}".')

            # Keep last output in shared context for custom callables.
            context["last_output"] = out

    # Keep worker robust by return error info instead of crashing batch execution.
    except Exception as exc:
        result["status"] = "error"
        result["error"] = str(exc)

    finally:
        # Explicit delete to release file/data resources early.
        if fiber is not None:
            del fiber

    return result


def run_batch_pipeline(batch, pipeline : list, fiber_kwargs : dict = None,
    validate_pipeline : bool = False)-> list[dict]:
    """Run pipeline processing over a batch of task rows.

    Parameters
    ----------
    batch : list
        List of task rows. Each row must be ``dict`` or ``pandas.Series``
        and contain at least ``"file"``.
    pipeline : list
        Pipeline specification (raw or normalized).
    fiber_kwargs : dict, optional
        Keyword arguments for ``Fiber`` initialization.
    validate_pipeline : bool, optional
        If ``True``, validates pipeline once before batch processing.

    Returns
    -------
    list[dict]
        List of outputs from ``run_file_pipeline`` (one dictionary per row).

    """

    # Basic batch check before processing.
    if not isinstance(batch, (list, tuple)):
        raise TypeError('"batch" must be a list or tuple of task rows.')

    # Validate once for the whole batch to avoid repeated checks.
    steps = pipeline_check(pipeline, allow_callables=True) if validate_pipeline else pipeline

    # Batch execution is usually called from Parallel.submit.
    # We silence per-file read bars here to keep terminal output clean.
    fiber_kwargs = {} if fiber_kwargs is None else dict(fiber_kwargs)
    if "show_progress" not in fiber_kwargs:
        fiber_kwargs["show_progress"] = False

    outputs = []

    # Process each row independently and collect results.
    for task_row in batch:
        outputs.append(
            run_file_pipeline(
                task_row=task_row,
                pipeline=steps,
                fiber_kwargs=fiber_kwargs,
                validate_pipeline=False
            )
        )

    return outputs


def window_map_tasks(file_window_map, columns=None):
    """Build worker task rows from a window-file mapping table.

    Parameters
    ----------
    file_window_map : pandas.DataFrame
        Mapping table generated by windowing methods.
    columns : list of str, optional
        Columns to keep in output task rows. If ``None``, a practical
        default set is used.

    Returns
    -------
    list[dict]
        Task rows ready to be consumed by ``run_batch_pipeline`` or
        ``Parallel.submit`` tasks.

    """

    if not isinstance(file_window_map, pd.DataFrame):
        raise TypeError('"file_window_map" must be a pandas.DataFrame.')

    if file_window_map.empty:
        return []

    # Default columns needed for most file- to window processing workflows.
    if columns is None:
        columns = ["file", "file_index", "window_id", "overlap_s", "file_weight", "window_weight"]

    if "file" not in file_window_map.columns:
        raise KeyError('"file_window_map" must contain "file" column.')

    # Keep only available requested columns to avoid crashes with lean maps.
    cols_available = [col for col in columns if col in file_window_map.columns]
    if "file" not in cols_available:
        cols_available = ["file"] + cols_available

    tasks_df = file_window_map[cols_available].copy()

    # Convert to list of dictionaries for backend-friendly serialization.
    return tasks_df.to_dict(orient="records")