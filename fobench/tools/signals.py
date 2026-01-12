'''
signal processing functions except filtering
'''

import numpy as np
import scipy.signal as signal
import scipy.integrate as integrate

from tqdm import tqdm

def hilbert(data: np.ndarray, axis: int = 0) -> np.ndarray:
    '''
    computes Hilbert transform for 2D array of signals
    Parameters
    ----------
    data : np.ndarray
        2D array holding data.
    axis : int, optional
        axis along which to operate. The default is 0.

    Returns
    -------
    ht : np.ndarray
        array holding analytic signals.

    '''

    return signal.hilbert(data, axis=axis)

def envelope(data: np.ndarray, axis: int = 0)-> np.ndarray:
    '''
    computes envelopes for 2D array of signals
    Parameters
    ----------
    data : np.ndarray
        2D array holding data.
    axis : int, optional
        axis along which to operate. The default is 0.

    Returns
    -------
    env : np.ndarray
        array holding signal envelopes.

    '''

    ht = hilbert(data, axis=axis)
    env = np.sqrt(data**2 + np.real(np.conjugate(ht)*ht)) #version of Christopher Wollin
	#env = (data**2 + ht**2)**0.5 #obspy version

    return env

def integrate_signal(data: np.ndarray, dx: int, axis: int)-> np.ndarray:
    '''
    integrates data in space or time using the cumulative trapezoid method

    Parameters
    ----------
    data : np.ndarray
        data to integrate
    dx : int
        time or space sampling.
    axis : int
        axis along which to apply operation.

    Returns
    -------
    np.ndarray
        integrated and detrended data

    '''

    integ = integrate.cumulative_trapezoid(y=data, dx=dx, axis=axis, initial=0) #+ self.data[0,:]
    return signal.detrend(integ, axis=axis)


def differentiate_signal(data: np.ndarray, method: str, axis: int, dt: int)-> np.ndarray:
    '''
	differentiates data in time or space

    Parameters
    ----------
    data : np.ndarray
        data to differentate.
    method : str
        sets the prefered method for differentiation. can be 'gradient' or 'diff',
        when using 'diff' data is prepended with inital value along specified axis.
    axis : int
        axis along which to apply operation.
    dt : int
        sampling period of data.

    Returns
    -------
    diff : np.ndarray
        differentiated data.

    '''

    if method == 'gradient':
        diff = np.gradient(data, dt, axis=axis)
    elif method == 'diff':
        padding = np.take(data, [0], axis=axis)
        diff = np.diff(data, axis=axis, prepend=padding) / dt

    return diff

def detrend_signal(o_signal: np.ndarray, order: int, axis:int =-1)-> np.ndarray:
    '''
    detrends signal(s)

    Parameters
    ----------
    o_signal : np.ndarray
        array containing the the signal(s)
    order : int
        order number of the fitting curve used to apply the detrend.
    axis : int, optional
        axis  alogn which to apply the detrending. The default is -1.

    Returns
    -------
    new_signal : np.ndarray
        detrended signal(s).
    '''

    o_signal = np.asarray(o_signal, dtype=float)
    new_signal = np.empty_like(o_signal)

    if o_signal.ndim == 1:
        t = np.arange(len(o_signal))
        trend = np.polyval(np.polyfit(t, o_signal, order), t)
        new_signal = o_signal - trend
        return new_signal

    elif o_signal.ndim == 2:
        n = o_signal.shape[axis]
        t = np.arange(n)

    for i in tqdm(range(o_signal.shape[1-axis]), desc='Detrending', leave=False):
        if axis == 0:
            y = o_signal[:, i]
            trend = np.polyval(np.polyfit(t, y, order), t)
            new_signal[:, i] = y - trend
        else:
            y = o_signal[i, :]
            trend = np.polyval(np.polyfit(t, y, order), t)
            new_signal[i, :] = y - trend

    return new_signal


def demean_signal(o_signal: np.ndarray, axis: int = None)-> np.ndarray:
    '''
    removes mean from signal(s)

    Parameters
    ----------
    o_signal : np.ndarray
        signal(s) to demean.
    axis : int, optional
        axis along which to apply demeaning. The default is None.

    Returns
    -------
    demeaned signal(s)
    '''

    o_signal = np.asarray(o_signal)

    if o_signal.ndim == 1:
        return o_signal - np.mean(o_signal)
    elif o_signal.ndim == 2:
        return o_signal - np.mean(o_signal, axis=axis, keepdims=True)


def get_tukey_window(M: int, alpha: float, sym: bool)-> np.ndarray:
    '''
    returns Tukey window for filtering and tapering

    Parameters
    ----------
    M : int
        number of points of window.
    alpha : float
        shape parameter, fraction of window inside tapered region, [0-1].
    sym : bool
        returns symmetric window when True.

    Returns
    -------
    np.ndarray
        Tukey window array.

    '''
    return signal.windows.tukey(M, alpha, sym)

def taper_signal(data: np.ndarray, axis: int, alpha: float = 0.1, detaper: bool = False) -> np.ndarray:
    '''
    taper signal using Tukey window

    Parameters
    ----------
    data : np.ndarray
        data array on which to perform tapering.
    alpha : float
        shape parameter, fraction of window inside tapered region, [0-1].
    axis : int
        axis along which to perform tapering.
    detaper : bool
        option to remove taper from signal

    Returns
    -------
    np.ndarray
        tapered data array.
    '''
    M = data.shape[axis]

    taper = get_tukey_window(M, alpha, sym=True)

    taper = taper[:, None] if axis == 0 else taper[None, :]
    if not detaper:
        return np.multiply(data, taper)
    elif detaper:
        return np.divide(data, taper)

def filt_preprocess(io_signal: np.ndarray, axis:int = None, order:int = 1,
                    sym: bool = True, alpha: float = 0.1,
                    steps: tuple[bool, bool, bool] = (True, True, True)) -> np.ndarray:
    '''
	performs pre-processing of the signal before filtering. Steps detrend, demean, and taper
    with Tukey window

    Parameters
    ----------
    io_signal : np.ndarray
        original input signal.
    axis : int, optional
        DESCRIPTION. The default is None.
    order : int, optional
        order number of the fitting curve used to apply the detrend. The default is 1.
    sym : bool, optional
        symmetric window for tapering if True. The default is True.
    alpha : float, optional
        shape parameter of taper, fraction of window inside tapered region,
        [0-1].The default is 0.1.
    steps : tuple[bool, bool, bool], optional
        defines which steps to perform. The default is True for all three steps.

    Returns
    -------
    io_signal : np.ndarray
        pre-processed signals.

    '''

    do_detrend, do_demean, do_taper = steps

    if do_detrend:
        io_signal = detrend_signal(io_signal, order, axis=axis)
    if do_demean:
        io_signal = demean_signal(io_signal, axis=axis)
    if do_taper:
        io_signal = taper_signal(io_signal, alpha=alpha, axis=axis)

    return io_signal

def normalize_signal(data: np.ndarray, method:str = 'absolute max',
                     axis: int = None, ram_window: int = None, fs: int = None,
                     total_channels: int = None, num_points: int = None)-> np.ndarray:
    '''
    Parameters
    ----------
    data : np.ndarray
        data array to normalize_signal.
    method : str, optional
        normalization method, options are:
			-'absolute max': with respect to the whole record (default)
			-'trace max': for each channel/timestep individually
			-'running mean': running absolute mean normalization
			-'1bit': 1-bit normalization.
        The default is 'absolute max'.
    axis : int, optional
        axis along which to perform operation. The default is None.
    ram_window : int, optional
        window length in seconds, only for running absolute mean normalization. The default is None.
    fs : int, optional
        sampling frequency of signal. The default is None.
    total_channels : int, optional
        total number of channels, only for running absolute mean normalization. The default is None.
    num_points : int, optional
        total number of time samples, only for running absolute mean normalization. The default is None.

    Raises
    ------
    TypeError
        running mean normalization chosen without window_length given.
    ValueError
        no valid normalization method chosen.

    Returns
    -------
    normalized_data : np.ndarray
        normalized data array
    '''

    if method == 'absolute max':
        normalized_data = (data - data.min()) / (data.max() - data.min())
        normalized_data = normalized_data * 2 -1

    elif method == 'trace max':
        channel_min = data.min(axis=axis, keepdims=True)
        channel_max = data.max(axis=axis, keepdims=True)
        normalized_data = (data - channel_min)/(channel_max - channel_min)
        normalized_data = normalized_data * 2 -1

    elif method == 'running mean':
        if ram_window is None:
            raise TypeError('please provide a window length for the running absolute mean normalization')

        normalized_data = data.copy()

        w_len = int(fs * ram_window)

        for i in tqdm(range(total_channels), desc='Running mean normalization', leave=False):
            for segment_start in range(num_points - w_len + 1):
                segment_end = segment_start + w_len
                segment = normalized_data[:,i][segment_start:segment_end]
                weight = np.mean(np.abs(segment))

                normalized_data[segment_start:segment_end, i] /= weight

    elif method == '1bit':
        normalized_data = np.sign(data).astype(np.float64)

    else:
        raise ValueError(f'"{method}" is not a valid normalization method')

    return normalized_data

def whiten_signal(data: np.ndarray, freq_min: int, freq_max: int, total_channels: int,
                  sampling_frequency: int)-> np.ndarray:
    '''
    adapted code from: https://github.com/seismo-live/seismo_live
    performs spectral whitening of all channels
    signals should be adequatly preprocessed


    Parameters
    ----------
    data : np.ndarray
        DESCRIPTION.
    freq_min,freq_max : int, int
        minimum and maximum of frequency band in which to perform spectral whitening.
    total_channels : int
        total number of channels.
    sampling_frequency : int
        sampling frequency of data

    Returns
    -------
    whitened_matrix : np.ndarray
        whitened data matrix

    '''

    whitened_matrix = np.zeros_like(data, dtype='float32')

    for i in tqdm(range(total_channels), desc='Whitening', leave=False):
  	  	channel = data[:, i]
  	  	n = len(channel)

  	  	f_range = float(freq_max) - float(freq_min)
  	  	nsmo = int(np.fix(min(0.01, 0.5 * f_range) * float(n) / sampling_frequency))
  	  	f = np.arange(n) * sampling_frequency / (n - 1.0)
  	  	JJ = ((f > float(freq_min)) & (f < float(freq_max))).nonzero()[0]

  	  	# channel FFT
  	  	FFTs = np.fft.fft(channel)
  	  	FFTsW = np.zeros(n, dtype=complex)

  	  	# apodization left
  	  	smo1 = (np.cos(np.linspace(np.pi / 2, np.pi, nsmo + 1)) ** 2)
  	  	FFTsW[JJ[0]:JJ[0] + nsmo + 1] = smo1 * np.exp(1j * np.angle(FFTs[JJ[0]:JJ[0] + nsmo + 1]))

  	  	# boxcar
  	  	FFTsW[JJ[0] + nsmo + 1:JJ[-1] - nsmo] = np.ones(len(JJ) - 2 * (nsmo + 1)) * \
  	  	np.exp(1j * np.angle(FFTs[JJ[0] + nsmo + 1:JJ[-1] - nsmo]))

  	  	# apodization to the right
  	  	smo2 = (np.cos(np.linspace(0.0, np.pi / 2.0, nsmo + 1)) ** 2)
  	  	espo = np.exp(1j * np.angle(FFTs[JJ[-1] - nsmo:JJ[-1] + 1]))
  	  	FFTsW[JJ[-1] - nsmo:JJ[-1] + 1] = smo2 * espo

  	  	# channel IFFT
  	  	whitedata = 2.0 * np.fft.ifft(FFTsW).real
  	  	whitened_matrix[:, i] = np.require(whitedata, dtype='float32')
    return whitened_matrix



def spectrum(o_signal, fs, pre_processing=True, order=1, pad=0, nfft=None,
             axis=None):
	'''
	Co-authors: --
	Description:
		Calculates de spectrum of a given signal with a specified sampling frequency.
	:Params:
		- signal(type:Numpy): original signal.
		- sampling_rate(type:Float): sampling rate of the signal.
		- pre_processing(type:Boolean): if True, signal will be detrend, demean, and tapered. Default is True.
		- order(type:int): polinomial order of the detrending curve.
		- pad(type:int): number of zeros to add to the signal before and after to increase num of points.
	:Return:
		- positive_freqs(type:Numpy): frequency values of the spectral curve.
        - magnitude(type:Numpy): amplitude values of spectral curve of the signal.
	'''

	if pre_processing:
		o_signal = filt_preprocess(io_signal=o_signal, order=order, axis=axis)

	o_signal = np.pad(o_signal, (pad-1, pad), mode='constant') if pad > 0 else o_signal # pad the signal to add points.

	n = len(o_signal) if nfft is None else nfft*len(o_signal)
	fft = np.fft.fft(o_signal, n=n)
	# Calculate the frequency axis
	freq_axis = np.fft.fftfreq(n, 1 / fs)

    # Take the positive frequencies and their corresponding magnitudes
	positive_freqs = freq_axis[:n//2]
	magnitude = 2/n * np.abs(fft)[:n//2]
	return positive_freqs, magnitude


def psd(o_signal, sampling_rate, pre_processing=True, order=1, n=None):
	'''
	Co-authors: --
	Description:
		Calculates de power spectral density based on the Welch method.
	:Params:
		- signal(type:Numpy): original signal.
		- sampling_rate(type:Float): sampling rate of the signal.
		- pre_processing(type:Boolean): if True, signal will be detrend, demean, and tapered. Default is True.
	:Return:
		- positive_freqs(type:Numpy): frequency values of the PSD curve.
		- magnitude(type:Numpy): amplitude values of PSD curve of the signal.
	'''

	if pre_processing == True:
		o_signal = filt_preprocess(o_signal, order)

	# We compute the PSD based on the Welch method.
	positive_freqs, magnitude = signal.welch(o_signal, sampling_rate, nperseg=n)

	return positive_freqs, magnitude