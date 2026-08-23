'''Functions to analyze wavefield characteristics'''

import numpy as np
from warnings import warn
from tqdm import trange
from obspy.signal.cross_correlation import correlate
from scipy.signal import correlate as scipy_correlate, correlation_lags
from fobench.tools import signals
from fobench.plotting import plotting_pyqt as plot_pyqt
from fobench.plotting.plotting_mpl import plot_acfs


def spatial_coherence_matrix(data: np.ndarray, max_lag: int, fs: int, distances: np.ndarray,
                             channel_nums: np.ndarray = None, plot_mode:str = 'pyqt',
                             results: bool = False, vmin: float = None, vmax: float = None) ->  np.ndarray:
    '''

    Parameters
    ----------
    data : np.ndarray
        array containing data to correlate, (n_channels, n_samples)
    max_lag : int
        maximum shift in seconds
    fs : int
        sampling frequency of signal
    distances : np.ndarray
        channel distances
    channel_nums : np.ndarray, optional
        channel numbering. The default is None.
    plot : bool, optional
        shows plot if True. The default is True.
    results : bool, optional
        whether to return the resulting values. The default is False.
    vmin, vmax : float, optional
        minimum and maximum limits of colorbar. The default is None.

    Returns
    -------
    coherence_matrix : np.ndarray
        correlation matrix

    '''

    n_ch, _ = data.shape
    data = (data - data.mean(axis=1, keepdims=True)) / data.std(axis=1, keepdims=True)
    coherence_matrix = np.zeros((n_ch, n_ch))
    for i in trange(n_ch, desc='Computing Coherence Matrix', leave=False):
        for j in range(i+1, n_ch):
            ccf = correlate(data[i], data[j], shift=int(max_lag*fs))
            max_corr = np.max(np.abs(ccf))
            coherence_matrix[i, j] = max_corr
            coherence_matrix[j, i] = max_corr
    np.fill_diagonal(coherence_matrix, 1.0)

    if plot_mode == 'mpl':
        warn('⚠️ matplotlib plotting not implemented for this method, '
             'plotting using pyqtgraph instead')
        plot_mode = 'pyqt'

    if plot_mode == 'pyqt':
        if not channel_nums: channel_nums = np.arange(0, n_ch)
        vmin = vmin if vmin is not None else -1
        vmax = vmax if vmax is not None else 1
        plot_pyqt.plot_2d_distance(distances=distances, data=coherence_matrix,
                         y_ticks=channel_nums, channels_num=channel_nums,
                         y_label = 'Channel #', x_label = 'Channel #',
                         title = 'Spatial Coherence Matrix',
                         cmap = 'viridis', cbar_label = 'Correlation Coefficient',
                         vmin=vmin, vmax=vmax)

    if results:
        return coherence_matrix

def autocorrelation_profile(data: np.ndarray, max_shift: int, axis: int, plot_mode: str,
                            deconvolve: bool, total_channels: int, distances: list, channels_num: list,
                            fs:int, window_size: int = None, vmin:float = None, vmax: float = None, **imshow_kwargs)-> np.ndarray:
    '''
    Computes the autocorrelation either for each channel or each time sample and
    optionally deconvolves the autocorrelation source term using a moving window.
    Deconvoltion is performed by substracting the average autocorrelation in a window or of the full record

    Parameters
    ----------
    data : np.ndarray
        input data.
    max_shift : int
        the maximum shift to use for the autocorrelation, when applied in time
        represents time samples, in space number of channels.
    axis : int
        axis where to apply the operation
    plot_mode : str
        pyqt or mpl plotting.
    deconvolve : bool
        deconvolution is applied if True.
    window_size : int, optional
        size of the moving window to compute the average autocorrelation, if not passed
        full record is used for deconvolution
    total_channels : int
        totalnumber of channels.
    distances : list, np.ndarray
        optical distances of channels.
    channels_num : list, np.ndarray
        channel numbering
    fs : int
        sampling frequency of data.
    vmin, vmax : float, optional
        minimum and maximum limits of colorbar. The default is None.
    **imshow_kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    result : np.ndarray
        autcorrelation profile

    '''

    autocorrelate = lambda x: correlate(x, x, max_shift)[max_shift:]
    result = np.apply_along_axis(autocorrelate, axis=axis, arr=data)

    if deconvolve:
        if window_size is None:
            avg_acf = np.mean(result, axis=1)
            avg_acf = avg_acf[:, np.newaxis]
        else:
            avg_acf = []
            half_window = window_size // 2
            for i in trange(total_channels, desc='Deconvolving', leave=False):
                start_ch = max(0, i - half_window)
                end_ch = min(total_channels, i + half_window + 1)
                avg_acf.append(np.mean(result[:, start_ch:end_ch], axis=1))
            avg_acf = np.array(avg_acf).T
        result -= avg_acf


    if plot_mode == 'pyqt':
        if vmin is None: vmin = result.min()
        if vmax is None: vmax = result.max()
        plot_pyqt.plot_2d_distance(distances=distances, data=np.rot90(result.T) if axis else np.rot90(result),
                         y_ticks=np.arange(0, max_shift)/fs, y_label = 'Lag/TWT [s]',
                         title = 'Autocorrelation Profile', cmap = 'viridis',
                         channels_num=channels_num, cbar_label = 'Correlation Coefficient',
                         invert_y=True, vmin=vmin, vmax=vmax)

    elif plot_mode == 'mpl':
        plot_acfs(acfs=np.rot90(result, k=-1) if axis else np.fliplr(result), distances=distances, fs=fs,
                  max_shift=max_shift, **imshow_kwargs)

    return result

def rmsa(data: np.ndarray, axis: int, window: int, dim: str, times: np.ndarray,
         channels_num: list, distances: np.ndarray, plot_mode: str, vmin: float = None,
         vmax: float = None) -> np.ndarray:
    '''
    computes the root mean square amplitude of data, data is split into windows
    before, pass window equal to data_length to compute one RMSA vector for full record

    Parameters
    ----------
    data : np.ndarray
        input data.
    axis : int
        axis along which to compute RMSA.
    window : int
        window lengths in time or space samples.
    dim : str
        dimension along which to compute
    times : np.ndarray
        array containing time stamps of data
    channels_num: list
        list of channel numbers
    distances : np.ndarray
        array containing optical distances
    plot_mode : str
        whether and hwo to plot
    vmin, vmax : float, float
        minimum and maximum of colorbar

    Returns
    -------
    np.ndarray
        array containing RMSA values.

    '''
    def compute_rmsa(data, data_length, window, axis):
        data = np.array_split(data, int(data_length / window), axis=axis)
        return np.array([np.sqrt(np.mean(window_data**2, axis=axis)) for window_data in data])

    if window is None:
        data_length = len(times) if dim == 't' else len(channels_num)
        result = compute_rmsa(data, data_length, data_length, axis)
    elif window:
        data_length = len(times) if dim == 't' else len(channels_num)
        result = compute_rmsa(data, data_length, window, axis)

    if plot_mode == 'mpl':
        warn('⚠️ matplotlib plotting not implemented for this method, '
                      'plotting using pyqtgraph instead')
        plot_mode = 'pyqt'

    if plot_mode == 'pyqt':
        if window is None:
            if dim == 't':
                plot_pyqt.plot_distance(distances=distances, data=result[0,:], channels_num=channels_num,
                              title='RMS Amplitude Profile', y_label='RMS Amplitude')
            elif dim == 'd':
                plot_pyqt.plot_timeseries(timestamps=times, data=result[0,:], dt=times[1]-times[0],
                                             y_label='RMS Amplitude', title='RMS Amplitude over Time')
        elif window:
            p95 = np.percentile(result, 95)
            if vmin is None: vmin = -p95
            if vmax is None: vmax = p95
            if dim == 't':
                timestamps = np.array_split(times, int(len(times)/window))
                timestamps = np.array([timestamp[int(len(timestamp)/2)] for timestamp in timestamps])
                plot_pyqt.plot_2d_timeseries(timestamps=timestamps, data=result, y_ticks=channels_num, y_label='Channel',
                                   dt=timestamps[1]-timestamps[0], title='RMS Amplitude', distances=distances,
                                   cmap='inferno', cbar_label='RMS Amplitude', vmin=vmin, vmax=vmax)
            elif dim == 'd':
                chans = np.array_split(channels_num, int(len(channels_num)/window))
                chans = np.array([chan[int(len(chan)/2)] for chan in chans])
                dists = np.array_split(distances, int(len(distances)/window))
                dists = np.array([dist[int(len(dist)/2)] for dist in dists])
                plot_pyqt.plot_2d_timeseries(timestamps=times, data=result.T, y_ticks=chans, y_label='Channel',
                                   dt=times[1]-times[0], title='RMS Amplitude', distances=dists,
                                   cmap='inferno', cbar_label='RMS Amplitude', vmin=vmin, vmax=vmax)

    return result

def peak_to_peak_amp(data: np.ndarray, fs: int, axis: int)-> tuple[np.ndarray, np.ndarray, np.ndarray]:
    '''
    finds peak-to-peak amplitude of data

    Parameters
    ----------
    data : np.ndarray
        input data.
    fs : int
        sampling frequency of data
    axis : int
        axis along which to apply operation.

    Returns
    -------
    pp_amp : np.ndarray
        peak-to-peak amplitudes.
    up_index : np.ndarray
            indices of minima.
    up_index : np.ndarray
            indices of maxima.
    '''

    data_dim = data.ndim
    peak_up, up_index = data.max(axis=axis), np.argmax(data, axis=axis)
    peak_down, down_index = data.min(axis=axis), np.argmin(data, axis=axis)
    pp_amp = peak_up - peak_down

    bad_picking = np.abs((up_index - down_index)) > fs/2
    if data_dim > 1:
        bad_picking = list(np.where(bad_picking)[0])

    if bad_picking:
        windows = [j for j in range(0, data.shape[axis], int(fs/4))]
        pp =  np.zeros(len(bad_picking)) if data_dim > 1 else 0

        for i in range(len(windows)-1):         #for pos in bad_picking:
            index, index_1 = windows[i], windows[i+1]

            if data_dim > 1:
                segment = data[index:index_1, bad_picking] if axis == 0 else data[bad_picking, index:index_1]
                new_pp = np.ptp(segment, axis=axis)
                pp[new_pp > pp] = new_pp[new_pp > pp]
            else:
                new_pp = np.ptp(data[index:index_1], axis=axis)
                pp = max(pp, new_pp)
            pp_amp[bad_picking] = pp if data_dim > 1 else pp_amp

    return pp_amp, up_index, down_index

def frequency_content(data: np.ndarray, fs:int, order: int, nfft: int, norm: bool,
                      axis: int):
    '''
    Parameters
    ----------
    data : np.ndarray
        data to transform.
    fs : int
        sampling frequency.
    order : int
        order of polynomial for detrending.
    nfft : int
        number of points to use for fft.
    norm : bool
        normalized each spectrum by its maximum value for better plotting
    axis : int
        axis along which to  perform the operation.

    Returns
    -------
    amp : np.ndarray
        fft amplitude values.
    freq : np.ndarray
        frequencies.
    '''

    data = signals.filt_preprocess(data, axis=axis, order=order)

    n = data.shape[axis] if nfft is None else nfft
    fft = np.fft.rfft(data, n=n, axis=axis)
    freq = np.fft.rfftfreq(n, 1/fs)
    amp = np.abs(fft) / fft.max() if norm else np.abs(fft)

    return amp, freq

def x_correlate(signal1: np.ndarray, signal2: np.ndarray, axis: int = -1,
                mode: str = 'full', demean: bool = True, normalize: bool = True,
                normalization: str = 'global', min_overlap: int = 1,
                gpu_use: bool = False) -> tuple[np.ndarray, np.ndarray]:
    '''
    Cross-correlates two signals (1D to N-D) along a selected axis.

    Parameters
    ----------
    signal1, signal2 : np.ndarray
        Input arrays. They can have any number of dimensions. Along ``axis``,
        lengths may differ. All other dimensions must be broadcast-compatible.
    axis : int, optional
        Axis along which cross-correlation is computed. The default is -1.
    mode : str, optional
        Correlation output mode:
        - ``full``: all overlaps, longest output (``n1 + n2 - 1`` samples).
        - ``same``: output length equals ``signal1.shape[axis]``.
        - ``valid``: only complete overlaps, shortest output.
        The default is ``full``.
    demean : bool, optional
        Removes mean from each trace before correlation. The default is True.
    normalize : bool, optional
        If True, correlation is divided by ``||x|| * ||y||``, keeping scores in
        ``[-1, 1]``. The default is True.
    normalization : str, optional
        Normalization strategy when ``normalize=True``:
        - ``global``: uses one denominator per trace (faster; can show edge spikes in
          ``mode='full'`` with different lengths).
        - ``windowed``: computes denominator per lag using only overlapping samples
          (Pearson-like per lag, more robust at edges).
        The default is ``global``.
    min_overlap : int, optional
        Minimum number of overlapping samples required at a lag. Lags with less
        overlap are set to 0. Useful to suppress unstable edge values in
        ``mode='full'``. The default is 1.
    gpu_use : bool, optional
        Kept for backwards compatibility. Currently unused.

    Returns
    -------
    lags : np.ndarray
        Lags (in samples) corresponding to ``corr`` along ``axis``.
    corr : np.ndarray
        Cross-correlation values. Shape equals the broadcasted input shape with
        length ``len(lags)`` along ``axis``.
    '''

    if isinstance(axis, bool):
        # Backward compatibility with the old signature x_correlate(s1, s2, gpu_use)
        gpu_use = axis
        axis = -1
    if gpu_use:
        warn('gpu_use is currently not implemented; running on CPU.')

    signal1 = np.asarray(signal1, dtype=float)
    signal2 = np.asarray(signal2, dtype=float)

    if signal1.ndim == 0 or signal2.ndim == 0:
        raise ValueError('signal1 and signal2 must be at least 1D arrays.')
    if signal1.ndim != signal2.ndim:
        raise ValueError('signal1 and signal2 must have the same number of dimensions.')
    if mode not in {'full', 'same', 'valid'}:
        raise ValueError("mode must be one of: 'full', 'same', 'valid'.")
    if normalization not in {'global', 'windowed'}:
        raise ValueError("normalization must be one of: 'global', 'windowed'.")
    if min_overlap < 1:
        raise ValueError('min_overlap must be >= 1.')

    if axis < -signal1.ndim or axis >= signal1.ndim:
        raise np.AxisError(axis, ndim=signal1.ndim)
    axis = axis % signal1.ndim
    data1 = np.moveaxis(signal1, axis, -1)
    data2 = np.moveaxis(signal2, axis, -1)

    try:
        batch_shape = np.broadcast_shapes(data1.shape[:-1], data2.shape[:-1])
        data1 = np.broadcast_to(data1, batch_shape + (data1.shape[-1],))
        data2 = np.broadcast_to(data2, batch_shape + (data2.shape[-1],))
    except ValueError as exc:
        raise ValueError(
            'signal1 and signal2 are not broadcast-compatible outside the selected axis.'
        ) from exc

    n1 = data1.shape[-1]
    n2 = data2.shape[-1]
    lags = correlation_lags(n1, n2, mode=mode)

    flat1 = data1.reshape(-1, n1)
    flat2 = data2.reshape(-1, n2)
    corr = np.empty((flat1.shape[0], lags.size), dtype=float)

    for i in range(flat1.shape[0]):
        x = flat1[i]
        y = flat2[i]

        if normalize and normalization == 'windowed':
            cc = np.zeros(lags.size, dtype=float)
            for j, lag in enumerate(lags):
                start = max(0, lag)
                stop = min(n1, n2 + lag)
                overlap = stop - start

                if overlap < min_overlap:
                    continue

                x_seg = x[start:stop]
                y_seg = y[start-lag:stop-lag]

                if demean:
                    x_seg = x_seg - x_seg.mean()
                    y_seg = y_seg - y_seg.mean()

                denom = np.linalg.norm(x_seg) * np.linalg.norm(y_seg)
                cc[j] = np.dot(x_seg, y_seg) / denom if denom > 0 else 0.0
        else:
            if demean:
                x = x - x.mean()
                y = y - y.mean()

            cc = scipy_correlate(x, y, mode=mode, method='auto')

            if normalize:
                denom = np.linalg.norm(x) * np.linalg.norm(y)
                cc = cc / denom if denom > 0 else np.zeros_like(cc)

        corr[i] = cc

    corr = corr.reshape(batch_shape + (lags.size,))
    corr = np.moveaxis(corr, -1, axis)

    return lags, corr