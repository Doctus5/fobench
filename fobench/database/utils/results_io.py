"""
Methods for writing and finalizing chunked processing results to zarr.

Created on 2026-04-09 17:47:00
Last modification on 2026-04-09 17:47:00

:author:
	- Sergio Diaz-Meza (sergioad@gfz-potsdam.de)
	- Jonas Pätzel (jonas.patzel@ulb.be)
:contributors:
	- Christopher Wollin (wollin@gfz-potsdam.de)
:license:

"""

# Necessary packages to import
import shutil
import numpy as np
import zarr


def init_zarr_store(store_path: str, windows_df, overwrite: bool = True):
    """Initialize a zarr store and save window coordinate vectors.

    Parameters
    ----------
    store_path : str
        Output zarr path.
    windows_df : pandas.DataFrame
        DataFrame with at least ``start_time`` and ``end_time`` columns.
    overwrite : bool, optional
        If ``True``, remove existing store before creation.

    Returns
    -------
    tuple
        Pair ``(zroot, n_windows)`` where ``zroot`` is the opened
        zarr group and ``n_windows`` is the number of windows.
    """

    if overwrite:
        shutil.rmtree(store_path, ignore_errors=True)

    zroot = zarr.open_group(store_path, mode="w")
    n_windows = int(windows_df.shape[0])

    # Save basic window coordinates in unix seconds.
    win_start_unix = windows_df["start_time"].astype("int64").to_numpy(dtype=np.int64) / 1e9
    win_end_unix = windows_df["end_time"].astype("int64").to_numpy(dtype=np.int64) / 1e9
    zroot.create_dataset("window_start_unix", data=win_start_unix, shape=win_start_unix.shape, dtype="f8")
    zroot.create_dataset("window_end_unix", data=win_end_unix, shape=win_end_unix.shape, dtype="f8")

    return zroot, n_windows


def init_weighted_accumulators(zroot, n_windows: int, n_channels: int,
    prefix: str = "result", window_chunk: int = 256, channel_offset: int = 0):
    """Create weighted sum/weight arrays and channel coordinates.

    Parameters
    ----------
    zroot : zarr.Group
        Open zarr group.
    n_windows : int
        Number of time windows.
    n_channels : int
        Number of channels.
    prefix : str, optional
        Prefix used for accumulator arrays.
    window_chunk : int, optional
        Chunk size along window axis.
    channel_offset : int, optional
        Channel index offset to save as metadata.

    Returns
    -------
    tuple
        Pair ``(z_sum, z_weight)`` zarr arrays.
    """

    chunk_0 = min(int(window_chunk), int(n_windows)) if n_windows > 0 else 1

    z_sum = zroot.create_dataset(
        f"{prefix}_sum_weighted",
        shape=(n_windows, n_channels),
        chunks=(chunk_0, n_channels),
        dtype="f8",
        fill_value=0.0,
    )
    z_weight = zroot.create_dataset(
        f"{prefix}_weight",
        shape=(n_windows, n_channels),
        chunks=(chunk_0, n_channels),
        dtype="f8",
        fill_value=0.0,
    )
    zroot.create_dataset(
        "channel_index",
        data=np.arange(n_channels, dtype=np.int32),
        shape=(n_channels,),
        dtype="i4",
    )

    # Keep useful metadata at group level.
    zroot.attrs["n_windows"] = int(n_windows)
    zroot.attrs["n_channels"] = int(n_channels)
    zroot.attrs["channel_offset"] = int(channel_offset)

    return z_sum, z_weight


def flush_weighted_updates(z_sum, z_weight, window_sum: dict, window_w: dict):
    """Flush one wave of in-memory weighted updates to zarr arrays.

    Parameters
    ----------
    z_sum : zarr.Array
        Weighted sum array with shape ``(n_windows, n_channels)``.
    z_weight : zarr.Array
        Weight array with shape ``(n_windows, n_channels)``.
    window_sum : dict
        Mapping ``window_id -> weighted_sum_vector``.
    window_w : dict
        Mapping ``window_id -> weight_vector``.

    Returns
    -------
    None
        Updates arrays in-place.
    """

    # Single-writer update loop (safe for streaming accumulation).
    for win_idx, s_add in window_sum.items():
        s_cur = z_sum[win_idx, :]
        w_cur = z_weight[win_idx, :]
        s_cur += s_add
        w_cur += window_w[win_idx]
        z_sum[win_idx, :] = s_cur
        z_weight[win_idx, :] = w_cur


def finalize_weighted_mean(zroot, z_sum, z_weight, out_name: str = "result_mean",
    window_chunk: int = 256, dtype: str = "f4", fill_value: float = np.nan):
    """Finalize weighted mean array from weighted sum and weight arrays.

    Parameters
    ----------
    zroot : zarr.Group
        Open zarr group where output will be written.
    z_sum : zarr.Array
        Weighted sum array with shape ``(n_windows, n_channels)``.
    z_weight : zarr.Array
        Weight array with shape ``(n_windows, n_channels)``.
    out_name : str, optional
        Output dataset name.
    window_chunk : int, optional
        Chunk size along window axis for output.
    dtype : str, optional
        Output dtype.
    fill_value : float, optional
        Fill value for windows/channels without data.

    Returns
    -------
    zarr.Array
        Output mean array.
    """

    n_windows, n_channels = z_sum.shape
    chunk_0 = min(int(window_chunk), int(n_windows)) if n_windows > 0 else 1

    z_mean = zroot.create_dataset(
        out_name,
        shape=(n_windows, n_channels),
        chunks=(chunk_0, n_channels),
        dtype=dtype,
        fill_value=fill_value,
    )

    for win_idx in range(n_windows):
        s_row = z_sum[win_idx, :]
        w_row = z_weight[win_idx, :]
        out = np.full(n_channels, fill_value, dtype=np.dtype(dtype))
        valid = w_row > 0
        if np.any(valid):
            out[valid] = (s_row[valid] / w_row[valid]).astype(np.dtype(dtype))
        z_mean[win_idx, :] = out

    return z_mean
