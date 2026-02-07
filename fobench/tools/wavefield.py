import numpy as np
from obspy.signal.cross_correlation import correlate
from tqdm import tqdm
from fobench.tools import signals
from fobench.plotting.plotting_pyqt import plot_2d_distance
from fobench.plotting.plotting_mpl import plot_acfs


def spatial_coherence_matrix(data: np.ndarray, max_lag: int, fs: int, distances: np.ndarray,
                             channel_nums: np.ndarray = None, plot: bool = True,
                             result: bool = False) ->  np.ndarray:
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
        DESCRIPTION. The default is False.

    Returns
    -------
    coherence_matrix : np.ndarray
        correlation matrix

    '''

    n_ch, n_samp = data.shape
    coherence_matrix = np.zeros((n_ch, n_ch))

    for i in tqdm(range(n_ch), desc='Computing Coherence Matrix', leave=False):
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

    if plot:
        if not channel_nums: channel_nums = np.arange(0, n_ch)
        plot_2d_distance(distances=distances, data=np.rot90(coherence_matrix),
                         y_ticks=channel_nums, channels_num=channel_nums,
                         max_value = None,
                         y_label = 'Channel #',
                         x_label = 'Channel #',
                         title = 'Spatial Coherence Matrix',
                         cmap = 'viridis',
                         cbar_label = 'Cross Correlation Coefficient')
    if result:
        return coherence_matrix

def autocorrelation_profile(data: np.ndarray, max_shift: int, axis: int, plot_mode: str,
                            deconvolve: bool, total_channels: int, distances: list, channels_num: list,
                            fs:int, window_size: int = None, **imshow_kwargs)-> np.ndarray:
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
            for i in range(total_channels):
                start_ch = max(0, i - half_window)
                end_ch = min(total_channels, i + half_window + 1)
                avg_acf.append(np.mean(result[:, start_ch:end_ch], axis=1))
            avg_acf = np.array(avg_acf).T
        result -= avg_acf


    if plot_mode == 'pyqt':
        plot_2d_distance(distances=distances, data=np.rot90(result),
                         y_ticks=np.arange(0, max_shift)/fs,
                         max_value = None,
                         y_label = 'Lag/TWT [s]',
                         title = 'Autocorrelation Profile',
                         cmap = 'viridis',
                         channels_num=channels_num,
                         cbar_label = 'Correlation Coefficient',
                         invert_y=True)

    elif plot_mode == 'mpl':
        plot_acfs(acfs=result, distances=distances, fs=fs,
                  max_shift=max_shift, **imshow_kwargs)

    return result

def rmsa(data: np.ndarray, axis: int, data_length: int, window: int) -> np.ndarray:
    '''
    computes the root mean square amplitude of data, data is split into windows
    before, pass window equal to data_length to compute one RMSA vector for full record

    Parameters
    ----------
    data : np.ndarray
        input data.
    axis : int
        axis along which to compute RMSA.
    data_length : int
        length of data in seconds.
    window : int
        window lengths.

    Returns
    -------
    np.ndarray
        array containing RMSA values.

    '''

    data = np.array_split(data, int(data_length / window), axis=axis)
    return np.array([np.sqrt(np.mean(window_data**2, axis=axis)) for window_data in data])


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