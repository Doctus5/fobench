import obspy as ob
from obspy.core.stream import Stream
#from obspy.core.trace import Trace
from obspy.core import UTCDateTime as UTC
from datetime import datetime


import numpy as np
import scipy.signal as signal
import scipy as scipy
import h5py as h5
#from numba import jit, cuda

#local paths
#from pyFunctions.PARAMS import *


#Function to sscan hierarchycally the HDF5 files.					
def scan_hdf5(path, recursive=True, tab_step=2):
    def scan_node(g, tabs=0):
        print(' ' * tabs, g.name)
        for k, v in g.items():
            if isinstance(v, h5.Dataset):
                print(' ' * tabs + ' ' * tab_step + ' -', v.name)
            elif isinstance(v, h5.Group) and recursive:
                scan_node(v, tabs=tabs + tab_step)
    with h5.File(path, 'r') as f:
        scan_node(f)
        
        
'''
TOOLS USED EXCLUSIVE FOR DAS CLASS
'''


# NEW DATA INSTRUMENT CORRECTION
#Does the simple instrumental correction for the DAS, Infrasound or BroadBand data by a simple per-count multiplication of instrumental factors. No instruemnt RESP file is needed.	Returns the corrected data.
def instr_corr(data=None, attributes=None, target='strain-rate'):
	'''
	Co-authors: --
	Description: 
		Applies an instrument correction to the data by converting the counts to the respective measuring unit depending on the instrument and data.
	:Params:
		- data(type:Numpy): original signal.
		- attributes(type:Dict): attributes dictionary of the Fiber Class.
	:Return:
		- data(type:Numpy): modified signal or data according to target unit.
	'''

	format, company, units = attributes['format'], attributes['company'], attributes['units']

	if format == 'tdms' and company == 'silixa': # Silixa TDMS
    
		if units == 'counts' and target == 'strain-rate':
    
			i_cst = 116E-9 # meters per radians.
			gauge_L = attributes['gauge_length'] # gauge lenght in meters.
			digital_N = 2**13 # magic number linked to the digitalization of the data. why not 2**16?
			fs = attributes['sampling_frequency'] # sampling frequency, which can be 1000 Hz for raw data.
			factor = i_cst*(fs/gauge_L)/digital_N # strain Rate per counts.
			data = np.multiply(data,factor)
     
            
    # if (format == 'h5' or format == 'hdf5') and company == 'febus': # FEBUS HDF5

	if (format == 'h5' or format == 'hdf5') and company == 'silixa': # Silixa HDF5

		if units == 'counts' and target == 'strain-rate':

			i_cst = 116E-9 # meters per radians.
			gauge_L = attributes['gauge_length'] # gauge lenght in meters.
			digital_N = 2**13 # magic number linked to the digitalization of the data. why not 2**16?
			fs = attributes['sampling_frequency'] # sampling frequency, which can be 1000 Hz for raw data.
			factor = i_cst*(fs/gauge_L)/digital_N # strain Rate per counts.
			data = np.multiply(data,factor)
        
    # if format == 'npy' and company == 'bam': # .npy format for BAM. This might fail always since the unit is NON-COMMERCIAL!

	if format == 'npz' and company == 'bam': # .npy format for BAM. This might fail always since the unit is NON-COMMERCIAL!
    
		if units == 'counts' and target == 'strain':
        
			factor = 1E-6 / 18.4 # strain per count (WEIRD!)
			data = np.multiply(data,factor)
    
	if (format == 'h5' or format == 'hdf5') and company == 'terra15': # Terra15 HDF5
        
		if units == 'm/s' and target == 'strain-rate':
    
			gauge_samples = int(round(attributes['gauge_length'] / attributes['spatial_interval']))
			print(gauge_samples)
			data = (data[:, gauge_samples:] - data[:, :-gauge_samples]) / (gauge_samples * attributes['spatial_interval'])
			n_decim = int(gauge_samples/2)
			attributes['channels'] = attributes['channels'][n_decim:-n_decim]
			attributes['channels_num'] = attributes['channels_num'][n_decim:-n_decim]
			attributes['total_channels'] = attributes['channels'].size

	if (format == 'h5' or format == 'hdf5') and company == 'asn': # ASN OptoDAS HDF5 (It can be a bit more complex, so I'm trying to make it simple!)

		if units == 'rad/(strain*m)' and target == 'strain-rate':

			data = data / attributes['conv_factor'] # divide by sensitivities. It seems they already provide the conversion factor

	if (format == 'h5' or format == 'hdf5') and company == 'quantx': # QuantX OptoaSense HDF5 (CHECK THIS!! WITH VERIFICATION OR CALIBRATION).

		if target == 'strain-rate':

			i_cst = 116E-9 # meters per radians.
			gauge_L = attributes['gauge_length'] # gauge lenght in meters.
			digital_N = int(attributes['units'][-4]) ** int(attributes['units'][-2:]) # magic number linked to the digitalization of the data. why not 2**ints
			fs = attributes['sampling_frequency'] # sampling frequency, which can be 1000 Hz for raw data.
			factor = i_cst*(fs/gauge_L)/digital_N # strain Rate per counts.
			data = np.multiply(data,factor)
		
	return data, target, attributes['channels'], attributes['channels_num'], attributes['total_channels']


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
	From Chirstopher Wollin version
	Computes the envelope of DAS in 2D.
	:Params(type):
		- data(numpy): matrix data (2D) of the DAS Class.
		- axis(int): axis for where to operate. Default 0: along the columns of DAS (per channel).
	:Return(type):
		- env(numpy): 2D analytic signal which is the envelope.  
	'''

	ht = hilbert(data, axis=axis)
	env = np.sqrt(data**2 + np.real(np.conjugate(ht)*ht)) #Chirstopher Wollin version
	#env = (data**2 + ht**2)**0.5 #obspy version
	
	return env


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


def detrend_signal(o_signal, order):
    
    t = np.arange(len(o_signal))
    new_signal = o_signal.astype(float)
    
    trend = np.polyval(np.polyfit(t, o_signal, order), t)
    new_signal -= trend
        
    return new_signal


def filt_preprocess(o_signal, order=1):
	'''
	Co-authors: --
	Description: 
		Do pre-processing of the signal for adecuate filtering. This includes detrend, demean, and tape in borders.
	:Params:
		- signal(type:Numpy): original signal.
		- order(type:Int): degree or order of the polyfit for detrending the signal. Default is 1.
	:Return:
		- new_signal(type:Numpy): pre-processed signal (detrended, demenaed, tapered).
	'''

	new_signal = detrend_signal(o_signal, order) # detrend signal
	new_signal -= new_signal.mean() # demean
	new_signal *= signal.windows.tukey(len(new_signal), alpha=0.05*2) # taper	

	return new_signal


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


def spatial_upsampling(das_class):
	'''
	Co-authors: --
	Description: Tool for upscaling spatialy the DDSS (DAS) data by double. Creates an interpolation between consecutive 
	channels to simulate an increase spatial resolution.
	:Params:
		- das_class(type:DAS): an initialized DAS Class with data.
	:Return:
		- new_data(type:numpy): 2D matrix containing the new spatial upsampled data.
		- new_channels_num(type_numpy): a list containing the new numbers of the channels, including the intermediate ones.  
	'''
	
	new_channels_num = [das_class.channels_num[0]]
	shape = (len(das_class.data[:,0]),1)
	new_data = das_class.data[:,0].reshape(shape) #reshaping is important to not affect the original dimensionality.
	
	for i in range(das_class.total_channels-1):
	
		first, second = new_data[:,-1].reshape(shape), das_class.data[:,i+1].reshape(shape)
		inter = (first + second) / 2
		new_data = np.concatenate((new_data, inter), axis=1)
		new_data = np.concatenate((new_data, second), axis=1)
				
		inter_num = (new_channels_num[-1] + das_class.channels_num[i+1]) / 2
		new_channels_num += ([inter_num, das_class.channels_num[i+1]])
		
	return new_data, new_channels_num
		
		
def spatial_downsampling(das_class):
	'''
	Co-authors: --
	Description: Tool for downscaling spatialy the DDSS (DAS) data by half. Erase one channel between consecutive 
	channels to simulate a decrease spatial resolution.
	:Params:
		- das_class(type:DAS): an initialized DAS Class with data.
	:Return:
		- new_data(type:numpy): 2D matrix containing the new spatial downsampled data.
		- new_channels_num(type_numpy): a list containing the new numbers of the channels, where the inermediate ones are eliminated.  
	'''
	
	new_data = das_class.data[:,::2]
	new_channels_num = das_class.channels_num[::2] #only if the label of the channel wants to be fixed (0,2,4,6,...,N)
	#new_channels_num = [i for i in range(0,int(len(das_class.channels_num)/2))] #channel numbers change due to the downsampling (0,1,2,3,...,N/2)
	#new_channels_num = das_class.channels_num[:int(len(das_class.channels_num)/2)] #channel numbers change due to the downsampling (0,1,2,3,...,N/2)
		
	return new_data, new_channels_num


















	
	
	
	
	


	 
