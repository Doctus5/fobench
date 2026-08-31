"""Shared processing methods for applying ``Fiber`` workflows over ``Dataset/Interrogator``
windows.

:Authors:
	- Sergio Diaz-Meza
	- Jonas Pätzel

:Contributors:
	- Christopher Wollin

"""

# Necessary packages to import
import gc
import os
import numpy as np

from .parallel import Parallel
from .utils.pipeline import pipeline_check, run_batch_pipeline
from .utils.results_io import (
    init_zarr_store,
    init_weighted_accumulators,
    flush_weighted_updates,
    finalize_weighted_mean,
)

def _auto_n_cores(cpu_ratio: float = 0.85):
    """Compute a safe automatic number of cores.

    Parameters
    ----------
    cpu_ratio : float, optional
        Fraction of available CPU cores to use. Defaults to ``0.85``.

    Returns
    -------
    int
        Suggested number of cores (at least 1).
    """

    total = os.cpu_count()
    total = 1 if total is None else total

    cpu_ratio = max(0.1, min(float(cpu_ratio), 1.0))

    return max(1, int(np.floor(total * cpu_ratio)))


def _normalize_pipeline(task, output_key: str = None):
    """Normalize task input into a validated pipeline and output key.

    Parameters
    ----------
    task : list, tuple, dict, str, or callable
        Task specification:

        - list/tuple: full pipeline.
        - dict with ``"pipeline"`` or ``"steps"``.
        - str/callable: shorthand single-step pipeline.

    output_key : str, optional
        Output key to capture and read from pipeline outputs.

    Returns
    -------
    tuple
        ``(pipeline, output_key)``.
    """

    # Convert supported input styles to a pipeline object.
    if isinstance(task, dict):
        if ("pipeline" in task) and ("steps" in task):
            raise ValueError('Use only one of "pipeline" or "steps" in task.')
        pipeline = task.get("pipeline", task.get("steps", None))
        if pipeline is None:
            raise ValueError('Dictionary task must contain "pipeline" or "steps".')
    elif isinstance(task, (list, tuple)):
        pipeline = task
    elif isinstance(task, str) or callable(task):
        pipeline = [{"op": task, "capture": "return", "store_as": output_key if output_key is not None else "result"}]
    else:
        raise TypeError('"task" must be pipeline list/tuple, dict, str, or callable.')

    pipeline = pipeline_check(pipeline, allow_callables=True)

    # Infer output key from first captured step if not given.
    if output_key is None:
        captured = [step["store_as"] for step in pipeline if step.get("capture", "none") != "none"]
        if len(captured) == 0:
            raise ValueError('Pipeline must capture at least one output ("capture" != "none").')
        output_key = captured[0]

    return pipeline, output_key


def _default_output_adapter(value):
    """Convert an operation output into one 1D numeric vector.

    Parameters
    ----------
    value : object
        Output value from a pipeline step.

    Returns
    -------
    numpy.ndarray
        1D vector.

    Raises
    ------
    ValueError
        If output can not be interpreted as one vector.
    """

    # Handle tuple/list outputs by trying first element recursively.
    if isinstance(value, (tuple, list)):
        if len(value) == 0:
            raise ValueError("Empty tuple/list output can not be converted to vector.")
        return _default_output_adapter(value[0])

    arr = np.asarray(value, dtype=float)
    if arr.ndim == 0:
        return np.array([float(arr)], dtype=float)
    if arr.ndim == 1:
        return arr
    if arr.ndim == 2 and arr.shape[0] == 1:
        return arr[0]

    raise ValueError(f"Unsupported output shape for vector adapter: {arr.shape}")


def run_window_pipeline_to_zarr(source, task, output_store: str, output_key: str = None,
    time_range: tuple = None, window_size: float = None, step: float = None,
    include_overlap: bool = True, min_overlap_s: float = 0.0,
    group_cols: list[str] = None, merge_datasets: bool = True,
    fiber_kwargs: dict = None, parallel_params: dict = None,
    submit_chunk_size: int = None, cpu_ratio: float = 0.85,
    show_progress: bool = True, reducer: str = "mean",
    output_adapter=None):
    """Apply a ``Fiber`` task over file windows and save results incrementally to zarr.

    Parameters
    ----------
    source : ``Dataset`` or ``Interrogator``-like object
        Source object that provides ``window_map``.
    task : list, tuple, dict, str, or callable
        Task definition to run over each file through pipeline execution.
    output_store : str
        Output zarr path.
    output_key : str, optional
        Key to extract from ``result["outputs"]``. If ``None``, inferred
        from first captured pipeline step.
    time_range : tuple, optional
        Time range as ``(start_time, end_time)``.
    window_size : float
        Window size in seconds.
    step : float, optional
        Window step in seconds. If ``None``, non-overlapping windows.
    include_overlap : bool, optional
        Keep overlapping files in time range if ``True``.
    min_overlap_s : float, optional
        Minimum overlap in seconds for file-window relations.
    group_cols : list[str], optional
        Grouping columns for compatibility grouping.
    merge_datasets : bool, optional
        Only used for Interrogator-like sources. If ``True``, map datasets together.
    fiber_kwargs : dict, optional
        Extra keyword arguments for Fiber initialization.
    parallel_params : dict, optional
        Parallel backend parameters.
    submit_chunk_size : int, optional
        Number of files processed per wave.
    cpu_ratio : float, optional
        Automatic CPU fraction for ``n_cores`` when not explicitly set.
    show_progress : bool, optional
        Show parallel submit progress bars.
    reducer : str, optional
        Reduction mode. Currently supports ``"mean"``.
    output_adapter : callable, optional
        Adapter ``adapter(value) -> 1D vector`` for operation output.

    Returns
    -------
    dict
        Summary dictionary with execution counters and output path.
    """

    if window_size is None:
        raise ValueError('"window_size" must be defined for windowed processing.')

    reducer = str(reducer).lower()
    if reducer != "mean":
        raise ValueError('Only reducer="mean" is currently supported.')

    pipeline, output_key = _normalize_pipeline(task=task, output_key=output_key)
    output_adapter = _default_output_adapter if output_adapter is None else output_adapter

    # Build complete map with metadata because file paths are needed for Fiber opening.
    window_kwargs = dict(
        time_range=time_range,
        window_size=window_size,
        step=step,
        include_overlap=include_overlap,
        min_overlap_s=min_overlap_s,
        group_cols=group_cols,
        include_meta=True,
        return_windows=True
    )
    if hasattr(source, "datasets"):
        window_kwargs["merge_datasets"] = merge_datasets

    file_window_map, windows_df = source.window_map(**window_kwargs)
    if file_window_map.empty:
        raise RuntimeError("Window map is empty. Nothing to process.")

    # Build one task per unique file to avoid duplicate Fiber computations.
    files_df = file_window_map[["file_index", "file"]].drop_duplicates().sort_values("file_index")
    task_rows = files_df.to_dict(orient="records")
    if len(task_rows) == 0:
        raise RuntimeError("No files available to process.")

    # Configure parallel backend with safe defaults.
    parallel_params = {} if parallel_params is None else dict(parallel_params)
    parallel_params.setdefault("mode", "multiprocessing")
    if "n_cores" not in parallel_params or parallel_params["n_cores"] is None:
        parallel_params["n_cores"] = _auto_n_cores(cpu_ratio=cpu_ratio)

    hpc = Parallel(params=parallel_params)
    n_cores = int(hpc.n_cores)

    # Determine wave size for chunked processing.
    if submit_chunk_size is None:
        submit_chunk_size = n_cores
    submit_chunk_size = max(1, int(submit_chunk_size))

    fiber_kwargs = {} if fiber_kwargs is None else dict(fiber_kwargs)
    zroot, n_windows = init_zarr_store(output_store, windows_df, overwrite=True)
    channel_offset = int(file_window_map["channel_offset"].iloc[0]) if "channel_offset" in file_window_map.columns else 0

    # Lazy init of output arrays after first valid output vector.
    n_channels = None
    z_sum = None
    z_weight = None

    # Summary counters.
    n_total = len(task_rows)
    n_ok = 0
    n_error = 0
    n_shape_skip = 0

    # Process file tasks in waves and flush each wave to zarr.
    for i0 in range(0, n_total, submit_chunk_size):
        i1 = min(i0 + submit_chunk_size, n_total)
        wave_tasks = task_rows[i0:i1]

        wave_batches = hpc.submit(
            task=run_batch_pipeline,
            var_args=wave_tasks,
            fixed_args=(pipeline, fiber_kwargs, False),
            to_return=True,
            submission_mode="single_submit",
            show_progress=show_progress
        )
        wave_results = [row for batch in wave_batches for row in batch]

        # Convert successful outputs to vectors keyed by file index.
        file_vectors = {}
        for row in wave_results:
            if row.get("status") != "ok":
                n_error += 1
                continue

            try:
                file_idx = int(row["file_index"])
                value = row["outputs"][output_key]
                vec = output_adapter(value)
                file_vectors[file_idx] = np.asarray(vec, dtype=float)
                n_ok += 1
            except Exception:
                n_error += 1

        if len(file_vectors) == 0:
            gc.collect()
            continue

        # Initialize output arrays once channel count is known.
        if n_channels is None:
            n_channels = int(next(iter(file_vectors.values())).size)
            z_sum, z_weight = init_weighted_accumulators(
                zroot=zroot,
                n_windows=n_windows,
                n_channels=n_channels,
                prefix=output_key,
                window_chunk=256,
                channel_offset=channel_offset
            )

        # Filter map relations only for files in this wave.
        wave_file_idx = set(file_vectors.keys())
        wave_map = file_window_map.loc[
            file_window_map["file_index"].isin(wave_file_idx),
            ["file_index", "window_id", "window_weight"]
        ]

        # In-memory accumulation for the current wave.
        window_sum = {}
        window_w = {}
        for rel in wave_map.itertuples(index=False):
            file_idx = int(rel.file_index)
            win_idx = int(rel.window_id)
            weight = float(rel.window_weight)

            if not np.isfinite(weight) or weight <= 0:
                continue

            vec = file_vectors.get(file_idx, None)
            if vec is None:
                continue
            if vec.size != n_channels:
                n_shape_skip += 1
                continue

            if win_idx not in window_sum:
                window_sum[win_idx] = np.zeros(n_channels, dtype=np.float64)
                window_w[win_idx] = np.zeros(n_channels, dtype=np.float64)

            valid = np.isfinite(vec)
            if not np.any(valid):
                continue

            window_sum[win_idx][valid] += vec[valid] * weight
            window_w[win_idx][valid] += weight

        flush_weighted_updates(
            z_sum=z_sum,
            z_weight=z_weight,
            window_sum=window_sum,
            window_w=window_w
        )

        # Free wave-level objects.
        del wave_batches, wave_results, file_vectors, wave_map, window_sum, window_w
        gc.collect()

    if n_channels is None:
        raise RuntimeError("No valid task outputs were produced.")

    finalize_weighted_mean(
        zroot=zroot,
        z_sum=z_sum,
        z_weight=z_weight,
        out_name=f"{output_key}_mean",
        window_chunk=256,
        dtype="f4",
        fill_value=np.nan
    )

    return {
        "output_store": output_store,
        "output_key": output_key,
        "n_files_total": n_total,
        "n_files_ok": n_ok,
        "n_files_error": n_error,
        "n_shape_skip": n_shape_skip,
        "n_windows": int(n_windows),
        "n_channels": int(n_channels),
        "parallel_mode": hpc.mode,
        "n_cores": int(n_cores),
        "submit_chunk_size": int(submit_chunk_size),
    }