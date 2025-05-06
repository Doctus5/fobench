import numpy as np
from obspy.signal.cross_correlation import correlate
from tqdm import tqdm
from fobench.plotting.plotting_pyqt import plot_2d_distance

def spatial_coherence(data: np.ndarray, max_lag: int, fs: int, channel_nums: np.ndarray = None,
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

    n_channels, n_samples = data.shape
    coherence_matrix = np.zeros((n_channels, n_channels))

    for i in tqdm(range(n_channels), desc='Computing Coherence Matrix'):
        for j in range(i, n_channels):
            ch_i = data[i]
            ch_j = data[j]
            
            # normalize signals
            ch_i = (ch_i - np.mean(ch_i)) / np.std(ch_i)
            ch_j = (ch_j - np.mean(ch_j)) / np.std(ch_j)

            ccf = correlate(ch_i, ch_j, shift=max_lag*fs)
            max_corr = np.max(np.abs(ccf))

            coherence_matrix[i, j] = max_corr
            coherence_matrix[j, i] = max_corr
    
    if plot:
        if not channel_nums: channel_nums = np.arange(0, n_channels)
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
