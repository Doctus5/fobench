"""
Windowing methods for handling windows over files in terms of times or space.

Created on 2026-04-06 16:05:46
Last modification on 2026-04-06 16:05:46

:author:
	- Sergio Diaz-Meza (sergioad@gfz-potsdam.de)
:contributors:
	- Jonas Pätzel (jonas.patzel@ulb.be)
	- Christopher Wollin (wollin@gfz-potsdam.de)
:license:

"""

# Necessary packages to import
import numpy as np
import pandas as pd


def build_time_windows(start_time, end_time, window_size, step=None):
    """Build fixed windows between two times.

    Parameters
    ----------
    start_time : object
        Start time of the full range.
    end_time : object
        End time of the full range.
    window_size : float or int
        Window size in seconds.
    step : float or int, optional
        Step between consecutive windows. If ``None``, non-overlapping windows
        are used with ``step = window_size``.

    Returns
    -------
    pandas.DataFrame
        Window table with columns:
        ``window_id``, ``start_time``, ``end_time``, ``center_time``,
        ``duration_s``.

    Raises
    ------
    ValueError
        If the provided time range is invalid or empty.
    """

    t0 = pd.to_datetime(start_time)
    tf = pd.to_datetime(end_time)

    if t0 >= tf:
        raise ValueError('"start_time" must be earlier than "end_time".')

    if not isinstance(window_size, (int, float)):
        raise TypeError('"window_size" must be int or float (seconds).')
    win = pd.to_timedelta(float(window_size), unit="s")
    if win <= pd.Timedelta(0):
        raise ValueError(f'"window_size" must be positive. Received: {window_size}.')

    if step is None:
        step = win
    else:
        if not isinstance(step, (int, float)):
            raise TypeError('"step" must be int or float (seconds).')
        step = pd.to_timedelta(float(step), unit="s")

    if step <= pd.Timedelta(0):
        raise ValueError(f'"step" must be positive. Received: {step}.')

    windows = []
    current = t0
    win_id = 0

    while current < tf:

        window_end = min(current + win, tf)

        if window_end <= current:
            break

        windows.append({
            "window_id": win_id,
            "start_time": current,
            "end_time": window_end,
            "center_time": current + (window_end - current) / 2,
            "duration_s": (window_end - current).total_seconds(),
        })

        current = current + step
        win_id += 1

    return pd.DataFrame(windows)


def overlap_seconds(start_a, end_a, start_b, end_b):
    """Calculate the overlap duration between two time intervals.

    Parameters
    ----------
    start_a : object
        Start time of interval A.
    end_a : object
        End time of interval A.
    start_b : object
        Start time of interval B.
    end_b : object
        End time of interval B.

    Returns
    -------
    float
        Overlap duration in seconds. Returns ``0.0`` when intervals do not
        overlap.
    """

    sa, ea = pd.to_datetime(start_a), pd.to_datetime(end_a)
    sb, eb = pd.to_datetime(start_b), pd.to_datetime(end_b)

    latest_start = max(sa, sb)
    earliest_end = min(ea, eb)

    overlap = (earliest_end - latest_start).total_seconds()
    return max(float(overlap), 0.0)


def map_files_to_windows(files_df, windows_df, start_label="start_time",
    end_label="end_time", min_overlap_s=0.0):
    """Map file intervals to windows using temporal overlap.

    Parameters
    ----------
    files_df : pandas.DataFrame
        File metadata table containing at least start/end time columns.
    windows_df : pandas.DataFrame
        Window metadata table from :func:`build_time_windows`.
    start_label : str, optional
        Column name for file start time.
    end_label : str, optional
        Column name for file end time.
    min_overlap_s : float or int, optional
        Minimum overlap in seconds to keep a file-window relation.

    Returns
    -------
    pandas.DataFrame
        Long-format mapping with columns:
        ``file_index``, ``window_id``, ``overlap_s``, ``file_weight``,
        ``window_weight``.

    Notes
    -----
    ``file_weight`` represents overlap divided by full file duration.
    ``window_weight`` represents overlap divided by full window duration.
    """

    if not isinstance(min_overlap_s, (int, float)):
        raise TypeError('"min_overlap_s" must be int or float (seconds).')
    min_overlap_s = float(min_overlap_s)

    required_file_cols = {start_label, end_label}
    if not required_file_cols.issubset(files_df.columns):
        missing = sorted(required_file_cols - set(files_df.columns))
        raise KeyError(f"Missing required file columns: {missing}")

    required_window_cols = {"window_id", "start_time", "end_time", "duration_s"}
    if not required_window_cols.issubset(windows_df.columns):
        missing = sorted(required_window_cols - set(windows_df.columns))
        raise KeyError(f"Missing required window columns: {missing}")

    files = files_df[[start_label, end_label]].copy()
    files[start_label] = pd.to_datetime(files[start_label])
    files[end_label] = pd.to_datetime(files[end_label])
    files["file_index"] = files.index
    files["file_duration_s"] = (files[end_label] - files[start_label]).dt.total_seconds()
    files = files.sort_values(by=[start_label, end_label]).reset_index(drop=True)

    windows = windows_df[["window_id", "start_time", "end_time", "duration_s"]].copy()
    windows["start_time"] = pd.to_datetime(windows["start_time"])
    windows["end_time"] = pd.to_datetime(windows["end_time"])
    windows = windows.sort_values(by=["start_time", "end_time"]).reset_index(drop=True)

    relations = []
    n_files = len(files)
    file_pointer = 0

    for _, win in windows.iterrows():

        win_start = win["start_time"]
        win_end = win["end_time"]
        win_duration = float(win["duration_s"])

        while file_pointer < n_files and files.loc[file_pointer, end_label] <= win_start:
            file_pointer += 1

        file_scan = file_pointer
        while file_scan < n_files and files.loc[file_scan, start_label] < win_end:

            file_start = files.loc[file_scan, start_label]
            file_end = files.loc[file_scan, end_label]
            file_duration = float(files.loc[file_scan, "file_duration_s"])

            ovlp = overlap_seconds(file_start, file_end, win_start, win_end)
            if ovlp > min_overlap_s:

                relations.append({
                    "file_index": int(files.loc[file_scan, "file_index"]),
                    "window_id": int(win["window_id"]),
                    "overlap_s": ovlp,
                    "file_weight": ovlp / file_duration if file_duration > 0 else np.nan,
                    "window_weight": ovlp / win_duration if win_duration > 0 else np.nan,
                })

            file_scan += 1

    return pd.DataFrame(relations)


def map_files_to_windows_grouped(files_df, windows_df, group_cols=None,
    start_label="start_time", end_label="end_time", min_overlap_s=0.0):
    """Map files to windows with optional grouping by acquisition parameters.

    Parameters
    ----------
    files_df : pandas.DataFrame
        File metadata table containing at least start/end time columns and
        optionally grouping columns.
    windows_df : pandas.DataFrame
        Window metadata table from :func:`build_time_windows`.
    group_cols : list of str, optional
        Columns used to split files into compatibility groups before mapping.
        If ``None`` or empty, behavior is identical to
        :func:`map_files_to_windows`.
    start_label : str, optional
        Column name for file start time.
    end_label : str, optional
        Column name for file end time.
    min_overlap_s : float or int, optional
        Minimum overlap in seconds to keep a file-window relation.

    Returns
    -------
    pandas.DataFrame
        File-window mapping. When grouping is enabled, output includes
        ``group_id`` and the grouping columns.
    """

    if group_cols is None or len(group_cols) == 0:
        return map_files_to_windows(
            files_df=files_df,
            windows_df=windows_df,
            start_label=start_label,
            end_label=end_label,
            min_overlap_s=min_overlap_s
        )

    missing = [col for col in group_cols if col not in files_df.columns]
    if missing:
        raise KeyError(f"Missing required grouping columns: {missing}")

    groups = files_df.groupby(group_cols, dropna=False)
    outputs = []
    expected_cols = ["file_index", "window_id", "overlap_s", "file_weight", "window_weight", "group_id", *group_cols]

    for group_id, (group_values, group_df) in enumerate(groups):

        group_map = map_files_to_windows(
            files_df=group_df,
            windows_df=windows_df,
            start_label=start_label,
            end_label=end_label,
            min_overlap_s=min_overlap_s
        )

        if group_map.empty:
            continue

        if not isinstance(group_values, tuple):
            group_values = (group_values,)

        group_map["group_id"] = int(group_id)
        for col_name, col_value in zip(group_cols, group_values):
            group_map[col_name] = col_value

        outputs.append(group_map)

    if not outputs:
        return pd.DataFrame(columns=expected_cols)

    return pd.concat(outputs, ignore_index=True)


def build_window_file_map(files_df, time_range, window_size, step=None,
    include_overlaps=True, min_overlap_s=0.0, group_cols=None,
    start_label="start_time", end_label="end_time", include_meta=True,
    return_windows=False):
    """Build a complete file-to-window map for a requested time range.

    Parameters
    ----------
    files_df : pandas.DataFrame
        File metadata table containing at least start/end time columns.
    time_range : tuple
        Two-value tuple ``(start_time, end_time)`` for the requested range.
    window_size : float or int
        Window size in seconds.
    step : float or int, optional
        Step between consecutive windows in seconds. If ``None``,
        non-overlapping windows are used.
    include_overlaps : bool, optional
        If ``True``, keep files that overlap the requested range.
        If ``False``, keep only files fully contained in the range.
    min_overlap_s : float or int, optional
        Minimum overlap in seconds for file-window relations.
    group_cols : list of str, optional
        Columns used to split files into compatibility groups before mapping.
    start_label : str, optional
        Column name for file start time.
    end_label : str, optional
        Column name for file end time.
    include_meta : bool, optional
        If ``True``, merge file and window metadata into the returned map.
    return_windows : bool, optional
        If ``True``, return both map and windows table.

    Returns
    -------
    pandas.DataFrame or tuple[pandas.DataFrame, pandas.DataFrame]
        File-window map. If ``return_windows=True``, returns
        ``(file_window_map, windows_df)``.
    """

    if len(time_range) != 2:
        raise ValueError('"time_range" must contain exactly two values: (start_time, end_time).')

    if start_label not in files_df.columns or end_label not in files_df.columns:
        raise KeyError(f'files_df must contain "{start_label}" and "{end_label}" columns.')

    t0 = pd.to_datetime(time_range[0])
    tf = pd.to_datetime(time_range[1])
    if t0 > tf:
        raise ValueError('"time_range[0]" must be earlier than or equal to "time_range[1]".')

    files = files_df.copy()
    files[start_label] = pd.to_datetime(files[start_label])
    files[end_label] = pd.to_datetime(files[end_label])

    if include_overlaps:
        keep = (files[start_label] <= tf) & (files[end_label] >= t0)
    else:
        keep = (files[start_label] >= t0) & (files[end_label] <= tf)

    selected = files.loc[keep].copy()
    selected = selected.sort_values(by=[start_label, end_label])

    windows_df = build_time_windows(start_time=t0, end_time=tf, window_size=window_size, step=step)

    if selected.empty:
        empty_map = pd.DataFrame(columns=["file_index", "window_id", "overlap_s", "file_weight", "window_weight"])
        if return_windows:
            return empty_map, windows_df
        return empty_map

    file_window_map = map_files_to_windows_grouped(
        files_df=selected,
        windows_df=windows_df,
        group_cols=group_cols,
        start_label=start_label,
        end_label=end_label,
        min_overlap_s=min_overlap_s
    )

    if include_meta and not file_window_map.empty:

        file_window_map = file_window_map.merge(
            selected,
            left_on="file_index",
            right_index=True,
            how="left",
            suffixes=("", "_file")
        )
        file_window_map = file_window_map.merge(
            windows_df[["window_id", "start_time", "end_time", "duration_s"]],
            on="window_id",
            how="left",
            suffixes=("_file", "_window")
        )

    if return_windows:
        return file_window_map, windows_df

    return file_window_map


def build_dataset_window_map(dataset, time_range, window_size, step=None,
    include_overlap=True, min_overlap_s=0.0, group_cols=None,
    include_meta=True, return_windows=False):
    """Build file-to-window mapping from a Dataset instance.

    Parameters
    ----------
    dataset : object
        Dataset-like object with ``database`` attribute.
    time_range : tuple
        Requested time range as ``(start_time, end_time)``.
    window_size : float or int
        Window size in seconds.
    step : float or int, optional
        Step between consecutive windows in seconds. If ``None``,
        non-overlapping windows are used.
    include_overlap : bool, optional
        If ``True``, keep files that overlap the requested range.
        If ``False``, keep only files fully contained in the range.
    min_overlap_s : float or int, optional
        Minimum overlap in seconds for file-window relations.
    group_cols : list of str, optional
        Columns used to split files into compatibility groups before mapping.
    include_meta : bool, optional
        If ``True``, merge file and window metadata into the output map.
    return_windows : bool, optional
        If ``True``, return both map and windows table.

    Returns
    -------
    pandas.DataFrame or tuple[pandas.DataFrame, pandas.DataFrame]
        Window mapping for the dataset. If ``return_windows=True``,
        returns ``(map_df, windows_df)``.
    """

    if not hasattr(dataset, "database"):
        raise TypeError('"dataset" must define a "database" attribute.')

    if dataset.database is None:
        raise RuntimeError("Dataset has no database. Build or load metadata first.")

    return build_window_file_map(
        files_df=dataset.database,
        time_range=time_range,
        window_size=window_size,
        step=step,
        include_overlaps=include_overlap,
        min_overlap_s=min_overlap_s,
        group_cols=group_cols,
        include_meta=include_meta,
        return_windows=return_windows
    )


def _resolve_group_columns(split_by_acquisition=False, group_cols=None):
    """Resolve compatibility grouping columns for file-window mapping.

    Parameters
    ----------
    split_by_acquisition : bool, optional
        If ``True``, default compatibility columns are used when
        ``group_cols`` is ``None``.
    group_cols : list of str, optional
        User-defined grouping columns.

    Returns
    -------
    list of str or None
        Grouping columns to use in mapping.
    """

    default_group_cols = [
        "sampling_frequency",
        "total_channels",
        "spatial_interval",
        "gauge_length",
        "channel_offset"
    ]

    if split_by_acquisition:
        return default_group_cols if group_cols is None else group_cols

    if group_cols is not None:
        return group_cols

    return None


def build_unit_window_map(unit, time_range, window_size, step=None,
    include_overlap=True, min_overlap_s=0.0, merge_datasets=True,
    split_by_acquisition=False, group_cols=None, include_meta=True,
    return_windows=False):
    """Build file-to-window mapping from a Unit instance.

    Parameters
    ----------
    unit : object
        Unit-like object with ``datasets`` attribute.
    time_range : tuple
        Requested time range as ``(start_time, end_time)``.
    window_size : float or int
        Window size in seconds.
    step : float or int, optional
        Step between consecutive windows in seconds. If ``None``,
        non-overlapping windows are used.
    include_overlap : bool, optional
        If ``True``, keep files that overlap the requested range.
        If ``False``, keep only files fully contained in the range.
    min_overlap_s : float or int, optional
        Minimum overlap in seconds for file-window relations.
    merge_datasets : bool, optional
        If ``True``, all datasets are mapped into one single timeline.
        If ``False``, each dataset is mapped independently and results are
        concatenated.
    split_by_acquisition : bool, optional
        If ``True``, files are split into compatibility groups before mapping.
    group_cols : list of str, optional
        Columns to define compatibility groups.
    include_meta : bool, optional
        If ``True``, merge file and window metadata into the output map.
    return_windows : bool, optional
        If ``True``, return both map and windows table.

    Returns
    -------
    pandas.DataFrame or tuple[pandas.DataFrame, pandas.DataFrame]
        Window mapping for the unit. If ``return_windows=True``,
        returns ``(map_df, windows_df)``.
    """

    if not hasattr(unit, "datasets"):
        raise TypeError('"unit" must define a "datasets" attribute.')

    group_cols_final = _resolve_group_columns(
        split_by_acquisition=split_by_acquisition,
        group_cols=group_cols
    )

    if merge_datasets is False:

        maps = []
        windows_ref = None

        for dataset_id, dataset in enumerate(unit.datasets):

            if getattr(dataset, "database", None) is None or dataset.database.empty:
                continue

            ds_map, ds_windows = build_dataset_window_map(
                dataset=dataset,
                time_range=time_range,
                window_size=window_size,
                step=step,
                include_overlap=include_overlap,
                min_overlap_s=min_overlap_s,
                group_cols=group_cols_final,
                include_meta=include_meta,
                return_windows=True
            )

            if windows_ref is None:
                windows_ref = ds_windows

            if ds_map.empty:
                continue

            ds_map = ds_map.copy()
            if "dataset_id" not in ds_map.columns:
                ds_map["dataset_id"] = dataset_id

            maps.append(ds_map)

        if not maps:
            empty = pd.DataFrame()
            if windows_ref is None:
                windows_ref = build_time_windows(
                    start_time=time_range[0],
                    end_time=time_range[1],
                    window_size=window_size,
                    step=step
                )
            if return_windows:
                return empty, windows_ref
            return empty

        merged = pd.concat(maps, ignore_index=True)
        if return_windows:
            return merged, windows_ref
        return merged

    parts = []
    for dataset_id, dataset in enumerate(unit.datasets):

        if getattr(dataset, "database", None) is None or dataset.database.empty:
            continue

        db = dataset.database.copy()
        db["dataset_id"] = dataset_id
        parts.append(db)

    if not parts:
        empty = pd.DataFrame()
        if return_windows:
            windows_ref = build_time_windows(
                start_time=time_range[0],
                end_time=time_range[1],
                window_size=window_size,
                step=step
            )
            return empty, windows_ref
        return empty

    files_df = pd.concat(parts, ignore_index=True)
    return build_window_file_map(
        files_df=files_df,
        time_range=time_range,
        window_size=window_size,
        step=step,
        include_overlaps=include_overlap,
        min_overlap_s=min_overlap_s,
        group_cols=group_cols_final,
        include_meta=include_meta,
        return_windows=return_windows
    )


def aggregate_file_vectors(file_vectors, file_window_map, n_windows=None,
    mode="sum", weight_label="overlap_s", fill_value=np.nan):
    """Aggregate per-file vectors into a windowed matrix.

    Parameters
    ----------
    file_vectors : array-like of shape (n_files, n_channels)
        Per-file vectors to aggregate (for example RMSA vectors).
    file_window_map : pandas.DataFrame
        Mapping generated by :func:`map_files_to_windows`.
    n_windows : int, optional
        Number of windows. If ``None``, it is inferred from
        ``file_window_map["window_id"]``.
    mode : {"rms", "mean", "sum"}, optional
        Aggregation mode:
        - ``"rms"``: overlap-weighted RMS combination
        - ``"mean"``: overlap-weighted mean
        - ``"sum"``: simple sum over mapped file vectors
    weight_label : str, optional
        Weight column from ``file_window_map`` used for ``"rms"`` and
        ``"mean"`` modes.
    fill_value : float, optional
        Value used for windows/channels with no contributions.

    Returns
    -------
    numpy.ndarray
        Matrix of shape ``(n_windows, n_channels)``.
    """

    values = np.asarray(file_vectors, dtype=float)
    if values.ndim != 2:
        raise ValueError('"file_vectors" must be a 2D array with shape (n_files, n_channels).')

    if file_window_map.empty:
        n_w = int(n_windows) if n_windows is not None else 0
        return np.full((n_w, values.shape[1]), fill_value, dtype=float)

    if n_windows is None:
        n_windows = int(file_window_map["window_id"].max()) + 1

    n_channels = values.shape[1]
    accum = np.zeros((n_windows, n_channels), dtype=float)
    weights = np.zeros((n_windows, n_channels), dtype=float)

    if mode not in {"rms", "mean", "sum"}:
        raise ValueError('"mode" must be one of: "rms", "mean", "sum".')

    if mode in {"rms", "mean"} and weight_label not in file_window_map.columns:
        raise KeyError(f'Missing weight column "{weight_label}" in file_window_map.')

    for rel in file_window_map.itertuples(index=False):

        file_idx = int(getattr(rel, "file_index"))
        win_idx = int(getattr(rel, "window_id"))
        vec = values[file_idx]
        valid = np.isfinite(vec)

        if not np.any(valid):
            continue

        if mode == "sum":

            accum[win_idx, valid] += vec[valid]
            weights[win_idx, valid] += 1.0

        else:

            w = float(getattr(rel, weight_label))
            if not np.isfinite(w) or w <= 0:
                continue

            if mode == "mean":
                accum[win_idx, valid] += vec[valid] * w
            elif mode == "rms":
                accum[win_idx, valid] += (vec[valid] ** 2) * w

            weights[win_idx, valid] += w

    result = np.full((n_windows, n_channels), fill_value, dtype=float)
    valid_out = weights > 0

    if mode == "sum":
        result[valid_out] = accum[valid_out]
    elif mode == "mean":
        result[valid_out] = accum[valid_out] / weights[valid_out]
    elif mode == "rms":
        result[valid_out] = np.sqrt(accum[valid_out] / weights[valid_out])

    return result
