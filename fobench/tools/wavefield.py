"""Functions to analyze wavefield characteristics."""

import numpy as np
from warnings import warn
from tqdm import trange
from obspy.signal.cross_correlation import correlate
from fobench.tools import signals
from fobench.plotting import plotting_pyqt as plot_pyqt
from fobench.plotting.plotting_mpl import plot_acfs


def spatial_coherence_matrix(data: np.ndarray, max_lag: int, fs: int, distances: np.ndarray,
                             channel_nums: np.ndarray = None, plot_mode: str = "pyqt",
                             results: bool = False, vmin: float = None, vmax: float = None) ->  np.ndarray:
    """Computes pairwise maximum cross-correlation between channels as a function of lag.

    Parameters
    ----------
    data : np.ndarray
        Array containing data to correlate, (n_channels, n_samples).
    max_lag : int
        Maximum shift in seconds.
    fs : int
        Sampling frequency of signal.
    distances : np.ndarray
        Channel distances.
    channel_nums : np.ndarray, optional
        Channel numbering.
    plot_mode : str
        ``"pyqt"``, ``"mpl"`` or no plotting.
    results : bool, optional
        Toggles return of computed values.
    vmin, vmax : float, optional
        Minimum and maximum limits of colorbar.

    Returns
    -------
    coherence_matrix : np.ndarray
        Correlation matrix.

    """

    n_ch, _ = data.shape
    data = (data - data.mean(axis=1, keepdims=True)) / data.std(axis=1, keepdims=True)
    coherence_matrix = np.zeros((n_ch, n_ch))
    for i in trange(n_ch, desc="Computing Correlation Matrix", leave=False):
        for j in range(i+1, n_ch):
            ccf = correlate(data[i], data[j], shift=int(max_lag*fs))
            max_corr = np.max(np.abs(ccf))
            coherence_matrix[i, j] = max_corr
            coherence_matrix[j, i] = max_corr
    np.fill_diagonal(coherence_matrix, 1.0)

    if plot_mode == "mpl":
        warn("⚠️ matplotlib plotting not implemented for this method, "
             "plotting using pyqtgraph instead")
        plot_mode = "pyqt"

    if plot_mode == "pyqt":
        if not channel_nums: channel_nums = np.arange(0, n_ch)
        vmin = vmin if vmin is not None else -1
        vmax = vmax if vmax is not None else 1
        plot_pyqt.plot_2d_distance(distances=distances, data=coherence_matrix,
                         y_ticks=channel_nums, channels_num=channel_nums,
                         y_label = "Channel #", x_label = "Channel #",
                         title = "Spatial Coherence Matrix",
                         cmap = "viridis", cbar_label = "Correlation Coefficient",
                         vmin=vmin, vmax=vmax)

    if results:
        return coherence_matrix

def autocorrelation_profile(data: np.ndarray, max_shift: int, axis: int, plot_mode: str,
                            deconvolve: bool, total_channels: int, distances: list, channels_num: list,
                            fs:int, window_size: int = None, vmin:float = None, vmax: float = None, **imshow_kwargs)-> np.ndarray:
    """Computes the autocorrelation either for each channel or each time sample and
    optionally deconvolves the autocorrelation source term using a moving window.
    Deconvoltion is performed by substracting the average autocorrelation in a window or of the full record.

    Parameters
    ----------
    data : np.ndarray
        Input data.
    max_shift : int
        The maximum shift to use for the autocorrelation, when applied in time
        represents time samples, in space number of channels.
    axis : int
        Axis where to apply the operation
    plot_mode : str
        ``"pyqt"``, ``"mpl"`` or no plotting.
    deconvolve : bool
        Deconvolution is applied if ``True``.
    window_size : int, optional
        Size of the moving window to compute the average autocorrelation, if not passed
        full record is used for deconvolution.
    total_channels : int
        Total number of channels.
    distances : list, np.ndarray
        Optical distances of channels.
    channels_num : list, np.ndarray
        Channel numbering.
    fs : int
        Sampling frequency of data.
    vmin, vmax : float, optional
        Minimum and maximum limits of colorbar.
    **imshow_kwargs :
        Additional kwargs passed to imshow.

    Returns
    -------
    result : np.ndarray
        Autcorrelation profile.

    """

    autocorrelate = lambda x: correlate(x, x, max_shift)[max_shift:]
    result = np.apply_along_axis(autocorrelate, axis=axis, arr=data)

    if deconvolve:
        if window_size is None:
            avg_acf = np.mean(result, axis=1)
            avg_acf = avg_acf[:, np.newaxis]
        else:
            avg_acf = []
            half_window = window_size // 2
            for i in trange(total_channels, desc="Deconvolving", leave=False):
                start_ch = max(0, i - half_window)
                end_ch = min(total_channels, i + half_window + 1)
                avg_acf.append(np.mean(result[:, start_ch:end_ch], axis=1))
            avg_acf = np.array(avg_acf).T
        result -= avg_acf


    if plot_mode == "pyqt":
        if vmin is None: vmin = result.min()
        if vmax is None: vmax = result.max()
        plot_pyqt.plot_2d_distance(distances=distances, data=np.rot90(result),
                         y_ticks=np.arange(0, max_shift)/fs, y_label = "Lag/TWT [s]",
                         title = "Autocorrelation Profile", cmap = "viridis",
                         channels_num=channels_num, cbar_label = "Correlation Coefficient",
                         invert_y=True, vmin=vmin, vmax=vmax)

    elif plot_mode == "mpl":
        plot_acfs(acfs=result, distances=distances, fs=fs,
                  max_shift=max_shift, **imshow_kwargs)

    return result

def rmsa(data: np.ndarray, axis: int, window: int, dim: str, times: np.ndarray,
         channels_num: list, distances: np.ndarray, plot_mode: str, vmin: float = None,
         vmax: float = None) -> np.ndarray:
    """Computes the root mean square amplitude of data. Data is split into windows
    before, pass window equal to data_length to compute one RMSA vector for full record.

    Parameters
    ----------
    data : np.ndarray
        Input data.
    axis : int
        Axis along which to compute RMSA.
    window : int
        Window lengths in time or space samples.
    dim : str
        Dimension along which to compute RMSA.
    times : np.ndarray
        Array containing time stamps of data.
    channels_num: list
        List of channel numbers.
    distances : np.ndarray
        Array containing optical distances.
    plot_mode : str
        ``"pyqt"``, ``"mpl"`` or no plotting.
    vmin, vmax : float, float
        Minimum and maximum of colorbar

    Returns
    -------
    np.ndarray
        Array containing RMSA values.

    """
    def compute_rmsa(data, data_length, window, axis):
        """Computes simple RMSA along an axis"""
        data = np.array_split(data, int(data_length / window), axis=axis)
        return np.array([np.sqrt(np.mean(window_data**2, axis=axis)) for window_data in data])

    if window is None:
        data_length = len(times) if dim == "t" else len(channels_num)
        result = compute_rmsa(data, data_length, data_length, axis)
    elif window:
        data_length = len(times) if dim == "t" else len(channels_num)
        result = compute_rmsa(data, data_length, window, axis)

    if plot_mode == "mpl":
        warn("⚠️ matplotlib plotting not implemented for this method, "
                      "plotting using pyqtgraph instead")
        plot_mode = "pyqt"

    if plot_mode == "pyqt":
        if window is None:
            if dim == "t":
                plot_pyqt.plot_distance(distances=distances, data=result[0,:], channels_num=channels_num,
                              title="RMS Amplitude Profile", y_label="RMS Amplitude")
            elif dim == "d":
                plot_pyqt.plot_timeseries(timestamps=times, data=result[0,:], dt=times[1]-times[0],
                                             y_label="RMS Amplitude", title="RMS Amplitude over Time")
        elif window:
            p95 = np.percentile(result, 95)
            if vmin is None: vmin = -p95
            if vmax is None: vmax = p95
            if dim == "t":
                timestamps = np.array_split(times, int(len(times)/window))
                timestamps = np.array([timestamp[int(len(timestamp)/2)] for timestamp in timestamps])
                plot_pyqt.plot_2d_timeseries(timestamps=timestamps, data=result, y_ticks=channels_num, y_label="Channel",
                                   dt=timestamps[1]-timestamps[0], title="RMS Amplitude", distances=distances,
                                   cmap="inferno", cbar_label="RMS Amplitude", vmin=vmin, vmax=vmax)
            elif dim == "d":
                chans = np.array_split(channels_num, int(len(channels_num)/window))
                chans = np.array([chan[int(len(chan)/2)] for chan in chans])
                dists = np.array_split(distances, int(len(distances)/window))
                dists = np.array([dist[int(len(dist)/2)] for dist in dists])
                plot_pyqt.plot_2d_timeseries(timestamps=times, data=result.T, y_ticks=chans, y_label="Channel",
                                   dt=times[1]-times[0], title="RMS Amplitude", distances=dists,
                                   cmap="inferno", cbar_label="RMS Amplitude", vmin=vmin, vmax=vmax)

    return result

def peak_to_peak_amp(data: np.ndarray, fs: int, axis: int)-> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Finds peak-to-peak amplitude of data.

    Parameters
    ----------
    data : np.ndarray
        Input data.
    fs : int
        Sampling frequency of data.
    axis : int
        Axis along which to apply operation.

    Returns
    -------
    pp_amp : np.ndarray
        Peak-to-peak amplitudes.
    up_index : np.ndarray
            Indices of maxima.
    down_index : np.ndarray
            Indices of minima.
    """

    data_dim = data.ndim
    peak_up, up_index = data.max(axis=axis), np.argmax(data, axis=axis)
    peak_down, down_index = data.min(axis=axis), np.argmin(data, axis=axis)
    pp_amp = peak_up - peak_down

    bad_picking = np.abs((up_index - down_index)) > fs/2
    if data_dim > 1:
        bad_picking = list(np.where(bad_picking)[0])

    if bad_picking:
        windows = [j for j in range(0, data.shape[0], int(fs/4))]
        pp =  np.zeros(len(bad_picking)) if data_dim > 1 else 0

        for i in range(len(windows)-1):         #for pos in bad_picking:
            index, index_1 = windows[i], windows[i+1]

            if data_dim > 1:
                new_pp = np.ptp(data[index:index_1,bad_picking], axis=axis)
                pp[new_pp > pp] = new_pp[new_pp > pp]
            else:
                new_pp = np.ptp(data[index:index_1], axis=axis)
                pp = max(pp, new_pp)
            pp_amp[bad_picking] = pp if data_dim > 1 else pp_amp

    return pp_amp, up_index, down_index

def frequency_content(data: np.ndarray, fs:int, order: int, nfft: int, norm: bool,
                      axis: int):
    """Computes the frequency-domain amplitude spectrum of detrended data via FFT.

    Parameters
    ----------
    data : np.ndarray
        Data to transform.
    fs : int
        Sampling frequency.
    order : int
        Order of polynomial for detrending.
    nfft : int
        Number of points to use for fft.
    norm : bool
        Normalize each spectrum by its maximum value for easier plotting.
    axis : int
        Axis along which to  perform the operation.

    Returns
    -------
    amp : np.ndarray
        FFT amplitude values.
    freq : np.ndarray
        Frequencies.
    """

    data = signals.filt_preprocess(data, axis=axis, order=order)

    n = data.shape[axis] if nfft is None else nfft
    fft = np.fft.rfft(data, n=n, axis=axis)
    freq = np.fft.rfftfreq(n, 1/fs)
    amp = np.abs(fft) / fft.max() if norm else np.abs(fft)

    return amp, freq