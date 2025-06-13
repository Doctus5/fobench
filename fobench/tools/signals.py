#Contains Signal Processing and Analysis Functionality except filtering

import numpy as np
import scipy.signal as signal
from tqdm import tqdm

def hilbert(data, axis=0):
	'''
	Computes the Hilbert transform in 2D.
	:Params(type):
		- data(type: numpy): matrix data (2D) of the DAS Class.
		- N(type: int): number of Fourier components. Default None: x.shape[axis].
		- axis(type: int): axis for where to operate. Default 0: along the columns of DAS (per channel).
	:Return(type):
		- ht(numpy): 1D analytic signal.  
	'''

	ht = signal.hilbert(data, axis=axis)
	return ht

def envelope(data, axis=0):
	'''
	Computes the envelope signal in 2D
	:Params(type):
		- data(numpy): matrix data (2D) of the DAS Class.
		- axis(int): axis for where to operate. Default 0: along the columns of DAS (per channel).
	:Return(type):
		- env(numpy): 2D analytic signal which is the envelope.
	'''

	ht = hilbert(data, axis=axis)
	env = np.sqrt(data**2 + np.real(np.conjugate(ht)*ht)) #version of Christopher Wollin
	#env = (data**2 + ht**2)**0.5 #obspy version
	
	return env

def detrend_signal(o_signal, order, axis=-1):
    
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


def demean_signal(o_signal, axis=None):

    o_signal = np.asarray(o_signal)
    
    if o_signal.ndim == 1:
        return o_signal - np.mean(o_signal)
    elif o_signal.ndim == 2:
        return o_signal - np.mean(o_signal, axis=axis, keepdims=True)

def get_tukey_window(N, alpha, sym): # alpha is the percentage of all the data contained/affected by the tapered window
    return signal.windows.tukey(N, alpha, sym)


def filt_preprocess(io_signal, axis=None, order=1, sym=True, percent=0.1, steps=(True, True, True)):
	'''
	Co-authors: --
	Description: 
		Do pre-processing of the signal for adecuate filtering. This includes detrend, demean, and tape in borders.
		User can define with 'steps' variable which processes want to do.
	:Params:
		- io_signal(type:Numpy): input original signal to be modified for output.
		- order(type:Int): degree or order of the polyfit for detrending the signal. Default is 1.
		- shape(type:String): Determines if the taper window has taper on both sides ('sym'=symmetrical),
		or just in one of the sides ('left' or 'right'). Default is 'sym'.
		- percent(type:Float): Percentage of the window/data where the taper is applied. Default is 0.1.
		- steps(Type:Tuple of Booleans): tuple or list of bolleans which determines by True or False the pre-processing to do,
		in the defined order -> (detrend, demean, taper). Default is True for all.
	:Return:
		- new_signal(type:Numpy): pre-processed signal (detrended, demenaed, tapered).
	'''

	if steps[0] == True:
		io_signal = detrend_signal(io_signal, order, axis=axis)
	if steps[1] == True:
		io_signal = demean_signal(io_signal, axis=axis)
	if steps[2] == True:
		window = get_tukey_window(io_signal.shape[0], alpha=percent, sym=sym)
		if io_signal.ndim == 2:
			window = window[:, np.newaxis]
		io_signal *= window
	return io_signal

def peak_to_peak_amp(data, sampl_freq, axis=0):
	'''
	Co-authors: --
	Description:
		Fiber Class method to find the peak to peak amplitude of the waveforms all across the channels.
	:Params:
		- data(type: numpy): numpy data fo DAS. rows are time steps and columns are channles/stations.
		- sampl_freq(type: float): sampling frequency of the data.
		- axis(type: int): axis to apply the method. Default is 0.
	:Return:
		- pp_amp(type: numpy): 1D array containing the peak to peak values per channel.
	'''

	peak_up, up_index = data.max(axis=axis), np.argmax(data, axis=axis)
	peak_down, down_index = data.min(axis=axis), np.argmin(data, axis=axis)
	pp_amp = peak_up - peak_down
 
	bad_picking = np.abs((up_index - down_index)) > sampl_freq/2
	bad_picking = list(np.where(bad_picking == True)[0])

	if bad_picking:
    
		windows = [j for j in range(0, data.shape[0], int(sampl_freq/4))]
		pp =  np.zeros(len(bad_picking))
    
		#for pos in bad_picking:
		for i in range(len(windows)-1):
    
			index, index_1 = windows[i], windows[i+1]
			new_pp = np.ptp(data[index:index_1,bad_picking], axis=axis)
			pp[new_pp > pp] = new_pp[new_pp > pp]

		pp_amp[bad_picking] = pp

	return pp_amp # try to return also indexes of maximum and minimum!

def spectrum(o_signal, sampling_rate, pre_processing=True, order=1, pad=0, nfft=None):
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

	if pre_processing == True:

		o_signal = filt_preprocess(o_signal, order)

	o_signal = np.pad(o_signal, (pad-1, pad), mode='constant') if pad > 0 else o_signal # pad the signal to add points.

	n = len(o_signal) if nfft is None else nfft*len(o_signal)
	fft = np.fft.fft(o_signal, n=n)

    # Calculate the frequency axis
	freq_axis = np.fft.fftfreq(n, 1 / sampling_rate)

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
