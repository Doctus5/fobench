import numpy as np
from obspy.signal.cross_correlation import correlate
from tqdm import tqdm
from fobench.plotting.plotting_pyqt import plot_2d_distance
from fobench.plotting.plotting_mpl import plot_acfs


def spatial_coherence_matrix(data: np.ndarray, max_lag: int, fs: int, channel_nums: np.ndarray = None,
                      plot: bool = True, result: bool = False) ->  np.ndarray:
    '''
    Co-authors: Jonas Pätzel
    Description: 
        computes and plots the maximum cross correlation coefficient for all 
        channel combinations
    :Params:
        - data(type: numpy): array containing data to correlate, (n_channels, n_samples)
        - max_lag(type: int): maximum shift in seconds
        - fs(type: int): sampling frequency of signal
        - result(type: bool): return result, default is False.
    :Return:
        - coherence_matrix(type: np.ndarray): spatial coherence matrix
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
        plot_2d_distance(distances=channel_nums, data=np.rot90(coherence_matrix), 
                         y_ticks=channel_nums,
                         max_value = None, 
                         y_label = 'Channel #',
                         x_label = 'Channel #',
                         title = 'Spatial Coherence Matrix',
                         cmap = 'viridis',
                         cbar_label = 'Cross Correlation Coefficient')
    if result:
        return coherence_matrix

def autocorrelation_profile(data: np.ndarray, max_shift: int, axis: int, dim:str, plot_mode: str, deconvolve: bool, 
                            window_size: int, total_channels: int, distances: list, fs:int, **imshow_kwargs)-> np.ndarray:

    autocorrelate = lambda x: correlate(x, x, max_shift)[max_shift:]    
    result = np.apply_along_axis(autocorrelate, axis=axis, arr=data)

    if deconvolve and dim=='t':
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
    elif deconvolve and dim=='d':
            raise ValueError('Source deconvolution only available in time dimension!')
    
    if plot_mode == 'pyqt':
        plot_2d_distance(distances=distances, data=np.rot90(result), 
                         y_ticks=np.arange(0, max_shift)/fs,
                         max_value = None, 
                         y_label = 'Lag/TWT [s]',
                         title = 'Autocorrelation Profile',
                         cmap = 'viridis',
                         cbar_label = 'Correlation Coefficient',
                         invert_y=True)
    return result

def fk_transform(data: np.ndarray, dt: float, dx: int, plot:str ='pyqt') -> tuple[np.ndarray, np.ndarray, np.ndarray]:

    n_samp, n_ch = data.shape
    fk_f= np.fft.rfft(data, axis=0)
    fk = np.fft.fftshift(np.fft.fft(fk_f, axis=1), axes=1)
    f = np.fft.rfftfreq(n_samp, dt)
    k = np.fft.fftshift(np.fft.fftfreq(n_ch, dx))
    
    if plot == 'pyqt':
        fk_plot(fk, f, k)
    
    return fk, f, k

def fk_plot(fk: np.ndarray, f: np.ndarray, k: np.ndarray)-> None: 
    pass