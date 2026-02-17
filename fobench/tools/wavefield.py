import numpy as np
from obspy.signal.cross_correlation import correlate
from tqdm import trange
from fobench.tools import signals
from fobench.plotting import plotting_pyqt as plot_pyqt
from fobench.plotting.plotting_mpl import plot_acfs
from warnings import warn



def spatial_coherence_matrix(data: np.ndarray, max_lag: int, fs: int, distances: np.ndarray,
                             channel_nums: np.ndarray = None, plot_mode:str = 'pyqt',
                             result: bool = False, vmin: float = None, vmax: float = None) ->  np.ndarray:
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
    result : bool, optional
        whether to return the resulting values. The default is False.
    vmin, vmax : float, optional
        minimum and maximum limits of colorbar. The default is None.

    Returns
    -------
    coherence_matrix : np.ndarray
        correlation matrix

    '''

    n_ch, n_samp = data.shape
    coherence_matrix = np.zeros((n_ch, n_ch))

    for i in trange(n_ch, desc='Computing Coherence Matrix', leave=False):
        for j in range(i, n_ch):
            ch_i = data[i]
            ch_j = data[j]

            # normalize signals
            ch_i = (ch_i - np.mean(ch_i)) / np.std(ch_i)
            ch_j = (ch_j - np.mean(ch_j)) / np.std(ch_j)
            ccf = correlate(ch_i, ch_j, shift=int(max_lag*fs))
            max_corr = np.max(np.abs(ccf))

            coherence_matrix[i, j] = max_corr
            coherence_matrix[j, i] = max_corr

    if plot_mode == 'pyqt':
        if not channel_nums: channel_nums = np.arange(0, n_ch)
        if vmin is None: vmin = coherence_matrix.min()
        if vmax is None: vmax = coherence_matrix.max()
        plot_pyqt.plot_2d_distance(distances=distances, data=np.rot90(coherence_matrix),
                         y_ticks=channel_nums, channels_num=channel_nums,
                         y_label = 'Channel #', x_label = 'Channel #',
                         title = 'Spatial Coherence Matrix',
                         cmap = 'viridis', cbar_label = 'Correlation Coefficient',
                         vmin=vmin, vmax=vmax)

    elif plot_mode == 'mpl':
        warn('matplotlib plotting not implemented for this method!')

    if result:
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
        plot_pyqt.plot_2d_distance(distances=distances, data=np.rot90(result),
                         y_ticks=np.arange(0, max_shift)/fs, y_label = 'Lag/TWT [s]',
                         title = 'Autocorrelation Profile', cmap = 'viridis',
                         channels_num=channels_num, cbar_label = 'Correlation Coefficient',
                         invert_y=True, vmin=vmin, vmax=vmax)

    elif plot_mode == 'mpl':
        plot_acfs(acfs=result, distances=distances, fs=fs,
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

    if plot_mode == 'pyqt':
        if window is None:
            if dim == 't':
                plot_pyqt.plot_distance(distances=distances, data=result[0,:], channels_num=channels_num,
                              title='RMS Amplitude Profile')
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

    elif plot_mode == 'mpl':
        warn('matplotlib plotting not implemented for this method!')

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