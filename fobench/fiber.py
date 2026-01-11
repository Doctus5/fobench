"""
Class "Fiber" for creating, storing and manipulating fiber optic sensing data.
So far it recieves TDMS format (Silixa) and H5 format (Febus).

Created on 2022-08-19 12:07:17
Last modification on 2024-06-28 19:17:00

:authors:
	- Sergio Diaz-Meza (sergioad@gfz-potsdam.de)
	- Jonas Pätzel (jonas.patzel@ulb.be)

:contributors:
	- Christopher Wollin (wollin@gfz-potsdam.de)
:license:

"""

#Necessary packages to import
import sys
import copy
from warnings import warn
from tqdm import tqdm

import numpy as np
import scipy.signal as signal
import scipy.integrate as integrate
from scipy.fft import fftshift, ifftshift, fft2, ifft2

from obspy.core import UTCDateTime as UTC

import pyqtgraph as pg
from pyqtgraph.Qt import QtWidgets

# inner functions
from .tools import file_io, utils, filters, signals, wavefield
from .plotting import plotting_mpl as plot
from .plotting import plotting_pyqt as plot_pyqt
from .plotting.pyqt_explorer import Explorer


class Fiber(object):
	'''
	IMPORTANT INFO: Most of the methods perform changes within the class permanently. Therefore it is useful to make a copy of the fiber instance
	with the method .copy() before performing any processing or changes.
	'''

	#Creates the basic variables of the DAS object with its characteristics
	def __init__(self, filepath, company=None, range_ch=None, sensing='das', load_data=True):
		'''
		Co-authors: --
		Description:
			Initializes a DAS Class which is reading a DAS file and saving all variables and metadata.
			The basis for manipulating the data is numpy. Tools are inspired in Obspy, however using
			an obspy class for this takes long time in their processing tools.
		:Params:
			- filepath(type:String): complete path of the file to be read.
			- company(type:String): manufacturer or the instrument that generates the data. Currently supporting 'silixa', 'febus'', 'bam', 'aragon', 'quantx','asn' and 'terra15'
			- range_ch(type:Int or List): channel number(s) to load only in data. Method to avoid loading all the data.
			- sensing(type:String): specifies the type of fiber optic sensing technique of the data. Default is 'das'.
		:Return:
			- NA.
		'''
		if not company:
			raise ValueError(
                '\nNo company provided! Please choose one of:\n'
                ' -"silixa"\n -"febus"\n -"bam"\n -"aragon"\n -"quantx"\n -"asn"\n'
                ' -"terra15"\n -"sintela"'
            )

		# Private attributes
		self.__filepath__ = filepath

		# Public attributes
		self.company = company
		self.format = filepath.split('.')[-1]

		self.attributes = file_io.read_data(self.__filepath__, self.company, range_ch, self.format, load_data=load_data)

		self.basefiles = [self.attributes['basefile']]
		self.fiber = self.attributes['fiber']
		self.properties = self.attributes['properties']
		self.channels = self.attributes['chans']
		self.channels_num = self.attributes['chans_nums']
		self.total_channels = self.attributes['list_chans_num']
		self.sampling_frequency = self.attributes['sampling_frequency'] # sampling rate of the data.
		self.o_sampling_frequency = self.attributes['o_sampling_frequency'] if self.attributes['o_sampling_frequency'] != None else self.attributes['sampling_frequency'] # original sampling frequency. Important for conversion factor.
		self.dt = 1 / self.sampling_frequency # calculated time step.
		self.start_time = self.attributes['start_time'] # start time of the data in file.
		self.end_time = self.attributes['end_time'] # end time of the data in file.
		self.spatial_interval = self.attributes['spatial_interval'] # channel spacing or spatial interval between channels [m].
		self.time_length = self.end_time - self.start_time
		self.num_points = self.attributes['num_points'] # int(self.time_length/self.dt)
		self.gauge_length = self.attributes['gauge_length'] # gauge length used in the measurement [m].
		self.channel_offset = self.attributes['channel_offset'] # offset where measurement started. It will not always record at channel 0 or distance 0.
		self.data = self.attributes['data']
		self.corrected = False
		self.sensing = sensing
		self.units = self.attributes['units']
		self.conv_factor = self.attributes['conv_factor'] # Extra variables (ONLY FOR ASN HDF5)
		self.processing = [{'instance creation' : UTC.utcnow().ctime()}]
		self.distances = [(num * self.spatial_interval) + self.channel_offset for num in self.channels_num]

		# Attributed not initialized since beginning. Requires further processing to be initialized
		self.ch_coord = None # coordinates of channels.

	'''
	####################################################
	Internal methods
	####################################################
	'''


	def __str__(self):
		'''
		defines output of print(Fiber); overview of most important recording parameters
		'''

		attributes = ['units', 'start_time', 'end_time', 'num_points', 'total_channels',
					'spatial_interval', 'sampling_frequency', 'gauge_length']

		return ('Instance of Fiber class\n'
                'recording parameters:\n'
                f'{"-" * 65}\n'+ "\n".join(f"{attr.ljust(25)} = {getattr(self, attr)}" for attr in attributes))

	__repr__ = __str__

	def __iadd__(self, other):
		'''
		allows to concatenate two Fiber instances with self += other
		'''
		if not isinstance(other, Fiber):
			raise TypeError('Object to add must be instance of Fiber class')
		return self.concatenate(other, fill_gaps=0)


	# Translates the string inout dimension into numerical axis of numpy.
	# If there is a change in how data matrix is opperated from now on, can be regulated from here instead throug all methods.
	def __axis__(self, dim):
		'''
		Co-authors: --
		Description:
			Translates a string input into numerical axial value for numpy. Used for the other methods.
		:Params:
			- dim(type: String): 't', or 'd' to differentiate between time or distance respectively.
		:Return:
			- axis(type: Int): integer value which denotes in numpy dimension where is time and distance.
		'''

		# Axis 0 in a matrix is the row dimension (downwards) and the 1 is column-wise (rightwards, elements inside each sub-array.)
		axial = {'t':0, 'd':1}

		return axial[dim]


	'''
	####################################################
	Standard methods
	####################################################
	'''


	def metadata(self, meta_dict=False):
		'''
		Co-authors: --
		Description:
			Print out the metadata in an organized way.
		:Params:
			- meta_dict(type:Boolean): if True, metadata is returned as dictionary. Default is False.
		:Return:
			- NA.
		'''

		if meta_dict:
			metainfo = {key: value for key, value in vars(self).items() if not key.startswith('__')}
			return metainfo

		else:
			for prop, value in self.properties.items():
				print(f"{prop} = {value}")

	def copy(self):
		'''
		Co-authors: --
		Description: Returns a deep copy of the class in the moment of execution.
		:Params:
			- NA.
		:Return:
			- (type:DAS Class): Same DAS Class in the state when the method is called.
		'''

		return copy.deepcopy(self)


	#Performs instrument correction on the data to get the strain rate values.
	def instr_correct(self, target='strain-rate'):
		'''
		Co-authors: --
		Description:
			Originally all data comes in counts. This method calls tools from another file to correct them to strain-rate (default).
			In future this can change.
		:Params:
			- NA.
		:Return:
			- NA.
		'''

		if self.corrected == False:

			self.data, self.units, self.channels, self.channels_num, self.total_channels = utils.instr_corr(self.data, vars(self), target=target)
			self.corrected = True

		return self


	#Slice the data on time. Waring, this affect the original data and prior time-lenghts can not be retrieved. See ".copy()" function. Input can be in ISO-format (String)
	def trim(self, t0=None, tf=None):
		'''
		Co-authors: --
		Description:
			Cuts the data in time between a given start-time and end-time. Updates properties of the class.
		:Params:
			- t0(type:UTC or String): start-time in UTC Class or string in ISOformat style.
			- tf(type:UTC or String): end-time in UTC Class or string in ISOformat style.
		:Return:
			- NA.
		'''

		t0, tf = UTC(t0), UTC(tf)

		# in case one of the triming times is beyond the range of the start and end times of the data, it redefines the limits to the ones of the data.
		t0 = max(t0, self.start_time)
		tf = min(tf, self.end_time)

		if tf < t0: raise ValueError("End time (tf) must be after start time (t0).")

		t = self.times()
		t0_pos = max(0, np.searchsorted(t, t0, side='right') - 1)
		tf_pos = max(0, np.searchsorted(t, tf, side='right') - 1)

		self.data = self.data[t0_pos:tf_pos, :]
		self.start_time = t[t0_pos]
		self.end_time = t[tf_pos]
		self.time_length = self.end_time - self.start_time
		self.num_points = self.data.shape[0]

		return self


	#Slice the data spatially (ranges of channels) by establshing the intial and the final channel. Anyformat of the channel code is acceptable.
	def restrict_channels(self, ch0, chf):
		'''
		Co-authors: --
		Description:
			Trims the data spatially. Only one channel is selected in the first channel is only specified. If the second channel is specified
			then the channels between the first and second specified are selected. Updates properties of the class.
		:Params:
			- ch0(type:Int or String): first channel ID to be selected.
			- chf(type:Int or String, optional): second channel ID to be selected.
		:Return:
			- NA.
		'''

		ch0, chf = int(min(ch0, chf)), int(max(ch0, chf)) # in case ch0 and chf not ordered
		ch0, chf = self.channels_num.index(ch0), self.channels_num.index(chf)
		self.data = self.data[:,ch0:chf+1]
		self.channels = self.channels[ch0:chf+1]
		self.channels_num = self.channels_num[ch0:chf+1]
		self.distances = self.distances[ch0:chf+1]
		self.total_channels = len(self.channels_num)

		return self


	#Attach X and Y coordinates (generally lon and lat) to the Fibre class.
	def append_coord(self, n_ch, x_ch, y_ch, z_ch):
		'''
		Co-authors: --
		Description:
			Attach coordinates of specified channels for the class instance. Necessary for plotting the fibre path and other spatial operations.
		:Params:
			- n_ch(type:Numpy): 1D array of channel number.
			- x_ch(type:Numpy): 1D array of X (longitude) coordinates of the channels specified in "n_ch".
			- y_ch(type:Numpy): 1D array of Y (latitude) coordinates of the channels specified in "n_ch".
			- z_ch(type:Numpy): 1D array of Z (depth - meters) coordinates of the channels specified in "n_ch".
		:Return:
			- NA.
		'''

		x_ch = np.zeros(n_ch.size) if x_ch is None else x_ch
		y_ch = np.zeros(n_ch.size) if y_ch is None else y_ch
		z_ch = np.zeros(n_ch.size) if z_ch is None else z_ch

		ch_coord = np.zeros((n_ch.size, 4))
		ch_coord[:,0] = n_ch
		ch_coord[:,1] = x_ch
		ch_coord[:,2] = y_ch
		ch_coord[:,3] = z_ch

		self.ch_coord = ch_coord

		return self

	#Attach X and Y coordinates (generally lon and lat) to the Fibre class.
	def georeference(self, n_ch, x_ch, y_ch, z_ch, system='decimal', err=None):
		'''
		Co-authors: --
		Description:
			Georeferencing of channels of the data. If tap tests were done to geolocate specific channels, this ones can be used for
			georeferencing other channels by linear interpolations between the located ones (assuming straight paths).
			It automatically attach the new coordinates to the Fiber class.
		:Params:
			- n_ch(type:Numpy): 1D array of channel number (ID) that were located with tap tests.
			- x_ch(type:Numpy): 1D array of X (longitude) coordinates of the channels specified in "n_ch".
			- y_ch(type:Numpy): 1D array of Y (latitude) coordinates of the channels specified in "n_ch".
			- z_ch(type:Numpy): 1D array of Z altitude [meters] of the channels specified in "n_ch".
			- system(type:String): Defined the receiving coordinate systems for X and Y. It can be 'decimal' for decimal degrees
			or 'utm' for Universal Transverse Mercator. Default is 'decimal'.
			- err(type:Float - Optional): maximum accepted error from interpolation in decimals. In case is given, the method will evaluate
			the error of the calculated channel spacing vs. the original one from the metadata.
		:Return:
			- NA.
		'''

		x_ch = np.zeros(n_ch.size) if x_ch is None else x_ch
		y_ch = np.zeros(n_ch.size) if y_ch is None else y_ch
		z_ch = np.zeros(n_ch.size) if z_ch is None else z_ch

		n_ch, x_ch, y_ch, z_ch = utils.interpolate_channels(n_ch, x_ch, y_ch, z_ch, system, err, self.spatial_interval) # georeferencing new channels between the tap tests points.

		self.append_coord(n_ch, x_ch, y_ch, z_ch)

		return self


	#Return the data with a channel specified if it's wanted
	def get_data(self, channel=None):
		'''
		Co-authors: --
		Description:
			Returns the data in the same way as using the attribute "data", however, this has an option for selecting the data corresponding to one channel.
		:Params:
			- channel(type:Int or Float): Channel number to get the data from. If not specified, it will be as same as adquiring the attribute "data".
		:Return:
			- data_n(type:Numpy): The data with the specific channel of interest, or the entire dataset.
		'''

		if channel is not None:

			ch = int(channel)
			index = self.channels_num.index(ch)
			return self.data[:,index]

		return self.data


	def times(self, time_type='UTCDateTime'):
		'''
        returns array of sample times, can be 'UTCDateTime', 'isoformat', 'matplotlib'
        or 'unix´'
	    '''
		return utils.return_times(self, time_type)


	#Function to concatenate 2 DAS classes. The concatenation will be done on the same class where the function is called. The 2 DAS objects must have the same characcteristics (sampling frequency, channels)
	def concatenate(self, input_das=None, fill_gaps=0):
		'''
		Co-authors: --
		Description:
			Concatenates 2 different DAS Classes, the one with the method called, and the one that enters as parameter.
			The code will identify the order of concatenation based on each class start and end-times. If there is an overlap, the data will be filled
			with one of them and continue filling with the other class once the overlap is finished.
			Updates variables and properties as consequence.
			Updates the Class itself to contain the data of both Classes.
			It is assumed that the host Class and the input Class have the same sampling rate and channels.
		:Params:
			- input_das(type:DAS Class): The Class to concatenate with. No matter if the class start-time is before or after the one to concatenate with.
			- fill_gaps(type:Int): If there is a gap between the 2 DAS Classes, then the gap will be filled with any specified value. Default = 0.0. Can also be np.nan
		:Return:
			- NA.
		'''

		axis = self.__axis__('t')

		if self.start_time <= input_das.start_time:

			first, second = self, input_das

		else:

			first, second = input_das, self

		tf = first.end_time + first.dt
		num_t = int((second.start_time + second.dt - first.end_time) / first.dt) - 1

# 		if num_t < 0:
# 			num_t = abs(num_t)
# 			second.data = second.data[num_t:,:]
#
# 		if num_t > 0:
# 			fill = np.zeros((num_t, first.total_channels))
# 			fill[fill==0] = np.nan if fill_gaps == None else fill_gaps #Can also work for putting NonType values (NaN) if fill_gaps is None or any value.
# 			first.data = np.concatenate((first.data, fill), axis=axis)

		self.data = np.concatenate((first.data, second.data), axis=axis)
		self.start_time = first.start_time
		self.end_time = second.end_time
		self.num_points = self.data.shape[axis]
		self.time_length = self.end_time - self.start_time
		self.basefiles.extend(input_das.basefiles)

		return self


	def to_traces(self, t_type='obspy'):
		'''
        returns channels as 'obspy' or 'pyrocko' streams, see fiber.tools.utils.to_traces
        for more details
		'''
		return utils.to_traces(self, t_type)


	'''
	####################################################
	Signal Processing methods
	####################################################
	'''

	#function for upsampling spatially by double/half depending if it is upsampling or downsampling spatially.
	@utils._update_processing
	def spatial_resample(self, rs_type=None):
		'''
		Co-authors: --
		Description:
			Affects the spatial resolution by adding to the data artifical channels between each channel (upsampling),
			or erases them in an interleaved order (downsampling).
			In case of upsampling, the artificial channels are made by inporlating the values in between.
			This method duplicates or divides by half the number of channels, and therefore the data. In case of wanting more spatial resolution,
			or less, the method must be applied several times.
			For future this can be changed an addapted, maybe to a desired channel spacing.
			For reference, in both upsamplig or downsampling, he first channel is always present (not affected).
			Updates properties and variables afterwards.
		:Params:
			- rs_type(type:String): Selects the mode for spatial resampling. Only 2 options possible: 1) 'upsampling' or 'downsampling'.
		:Return:
			- NA.
		'''

		if rs_type == 'upsampling':

			print('Upsampling takes longer than downsampling. It might take a while...')
			new_data, new_channels_num = utils.spatial_upsampling(self)
			self.spatial_interval = self.spatial_interval / 2

		elif rs_type == 'downsampling':

			new_data, new_channels_num = utils.spatial_downsampling(self)
			self.spatial_interval = self.spatial_interval * 2

		#if (rs_type != 'downsampling') & (rs_type != 'upsampling'):
		else:

			raise ValueError('Spatial resampling type is not recognizable. Only (upsampling) or (downsampling).')

		self.data = new_data
		self.channels_num = new_channels_num
		self.total_channels = len(self.channels_num)

		# print(data.shape, )

		return self

	@utils._update_processing
	def detrend(self, order=1, dim='t'):
		'''
        detrend signal with specified order polynomial
        see fiber.tools.signals.detrend_signal for more details
	    '''

		axis = self.__axis__(dim)
		self.data = signals.detrend_signal(self.data, order=order, axis=axis)

		return self

	@utils._update_processing
	def demean(self, dim='t'):
		'''
        remove mean of signal along specified dimension
        see fiber.tools.signals.demean_signal for more details
		'''

		axis = self.__axis__(dim)
		self.data = signals.demean_signal(self.data, axis=axis)

		return self

	@utils._update_processing
	def taper(self, alpha=0.05, dim='t', detaper=False):
		'''
        taper data
        see fiber.tools.signals.taper_signal for more details
		'''

		axis = self.__axis__(dim)
		self.data = signals.taper_signal(data=self.data, axis=axis,
                                   alpha=alpha, detaper=detaper)

		return self


	#Function to decimate the data by any frequency below the original. Carefull when applying decimations with factors over or equal to 13, then better to call decimation twice (see scipy.signal.decimate for more...). The decimation function of scipy performs a pre-filtering process to avoid anti-aliasing on the signals.
	@utils._update_processing
	def decimate(self, new_freq=None, ftype='fir-remez'):
		'''
		Co-authors: --
		Description:
			Decimates the data by any frequency below the original (preferably a divisible one).
			Carefull when applying decimations with factors over or equal to 13, then better to call decimation twice (see scipy.signal.decimate for more...).
			The decimation function of scipy performs a pre-filtering process to avoid anti-aliasing on the signals.
		:Params:
			- new_frequency(type:Int or Float): the new sampling frequency or sampling rate of the decimated data.
			- ftype(type:String): There are 3 types of available pre-filters: 1) "fir-remez" Marius Isken adptative antialiasing filter, 2)
			"fir235" Javier Quinteros designed filter for DAS, and if "None", then a order 8 Chebyshev type I filter is used.
		:Return:
			- NA.
		'''

		axis = self.__axis__('t') # axis to apply, default is time.

		if self.sampling_frequency % new_freq != 0:
			warn(f'Decimation to {new_freq} Hz not possible! Decimating to {self.sampling_frequency / int(self.sampling_frequency / new_freq)} Hz instead')
		down_factor = int(self.sampling_frequency / new_freq)
		new_freq = self.sampling_frequency / down_factor

		#Check prefilter... which one is?
		if ftype is not None:

			new_data = filters.decimate(data=self.data, factor=down_factor, ftype=ftype, axis=axis)

		else:

			new_data = signal.decimate(x=self.data, q=down_factor, axis=axis)

		self.data = new_data
		self.sampling_frequency = new_freq
		self.dt = 1 / self.sampling_frequency
		self.num_points = self.data.shape[0]

		return self


	@utils._update_processing
	def normalize(self, method='absolute max', dim='d', ram_window=None):
		'''
	    normalize data, methods are 'absolute max', 'trace max', 'running mean' and '1bit'
        see fiber.tools.signals.normalize_signal for more details
	    '''

		axis = self.__axis__(dim)
		self.data = signals.normalize_signal(self.data, method=method,
                                       ram_window=ram_window, axis=axis,
                                       fs=self.sampling_frequency, total_channels=self.total_channels,
                                       num_points=self.num_points)

		return self

	@utils._update_processing
	def whiten(self, freq_min=0.01, freq_max=100):
	    '''
	    spectral whitening of data
        see fiber.tools.signals.whiten_signal for more details
	    '''
	    if not any('filter' in preprocessing for preprocessing in self.processing):
        	    	  warn('Data has possibly not been filtered before whitening! Check'
                     'preprocessing and results carefully!\ncontinuing...')

	    self.data = signals.whiten_signal(freq_min=freq_min, freq_max=freq_max,
                                    sampling_frequency=self.sampling_frequency,
                                    total_channels=self.total_channels)

	    return self

	# Function for filtering. Shall we also declare dimensionality options here?
	@utils._update_processing
	def filter(self, f_type=None, freq=None, pre_process=True, frac=0.05, order=1, **options):
		'''
		Co-authors: --
		Description:
			Filters the data based on a specified type of filter (lowpass, bandpass, highpass, bandstop, median) and the values. Filters are based on Obspy codes
			and so the multiple options are also.
		:Params:
			- f_type(type:String): type of filter to apply. Options are: 'lowpass', 'bandpass', 'highpass', 'bandstop' and 'median'
			- freq(type:Int, Float, Tuple): cut-off value for the filter in 'lowpass' and 'highpass'. If it for 'bandpass', then it must be a tuple containing the cut-offs of the bandwith.
			- pre_process(type:Boolean): to use to detren, demean and tape before filtering. Default is True. Automatically set to False when using median filter
			- frac(type:Float): see description in method "taper".
			- order(type:Int): order of polynomuial fit for detrending.
		:Return:
			- NA.
		'''
		if f_type == 'median':
			pre_process = False
		if pre_process:
			self.data = signals.filt_preprocess(self.data, axis=0)

		new_data = filters.point_filter(f_type=f_type, data=self.data, df=self.sampling_frequency, freq=freq, **options)
		self.data = new_data

		return self


	# Function to make a fk-filter based on the input parameters.
	@utils._update_processing
	def fk_filter(self, freq_min, freq_max, k_min, k_max):
		'''
		Co-authors: --
		Description:
			under construction. DON'T USE THIS METHOD!
		:Params:
			- param1(type:--): --.
		:Return:
			- return1(type:--): --.
		'''

		data_fk = fftshift(fft2(ifftshift(self.data)))

		# Define the frequency and wavenumber grids
		#num_rows, num_cols = self.data.shape
		freq_grid = np.fft.fftfreq(self.num_points)
		k_grid = 2 * np.pi * np.fft.fftfreq(self.total_channels, d=self.spatial_interval)#.reshape((1,self.total_channels))

		freq_mesh, k_mesh = np.meshgrid(freq_grid, k_grid, indexing='ij')

		# Define the filter mask
		mask = (np.abs(freq_mesh) >= freq_min) & (np.abs(freq_mesh) <= freq_max) & (np.abs(k_mesh) >= k_min) & (np.abs(k_mesh) <= k_max)
		filt_fk = data_fk * mask

		# Apply inverse 2D Fourier transform to obtain the filtered data
		filt_data = np.abs(fftshift(ifft2(ifftshift(filt_fk))))

		return np.abs(data_fk) #filt_data


	#Function for integrating the signal
	@utils._update_processing
	def integrate(self, method='cum_trapezoid', dim='t', taper=True):
		'''
		Co-authors: --
		Description:
			Integrates the data in time (dim='t') or in space (dim='d').
		:Params:
			- method(type: String): sets the prefered method for integration. Default and only one is "cum_trapezoid" as cumulative trapezoid
			- dim(type: String): dimension to where to apply the operation ('t' = time, 'd' = space). Default is 't'.
			- taper(type: Boolean): Decided wether to apply a taper to the signal before integration to avoid offsets or trends (recommended).
			Default is True.
		:Return:
			- NA.
		'''

		axis = self.__axis__(dim)
		dx = self.dt if axis == 0 else self.spatial_interval

		if taper:
			self.taper(dim=dim)

		if method == 'cum_trapezoid':
			result = integrate.cumulative_trapezoid(y=self.data, dx=dx, axis=axis, initial=0) #+ self.data[0,:]

		result = signal.detrend(result, axis=axis) #to detrend the signal

		#if taper == True:

		#	self.detaper(axis=axis)

		self.data = result

		return self


	#Function for differentiating the signal
	@utils._update_processing
	def differentiate(self, method='gradient', dim='t'):
		'''
		Co-authors: --
		Description:
			Differentiates the data in time (dim='t') or in space (dim='d').
		:Params:
			- method(type:String): sets the prefered method for differentiation. can be 'gradient' or 'diff', when using 'diff' data is prepended with inital value along specified axis
			- dim(type: String): dimension to where to apply the operation ('t' = time, 'd' = space). Default is 't'.
		:Return:
			- NA.
		'''

		axis = self.__axis__(dim)
		if method == 'gradient':
			result = np.gradient(self.data, self.dt, axis=axis)
		elif method == 'diff':
			padding = np.take(self.data, [0], axis=axis)
			result = np.diff(self.data, axis=axis, prepend=padding) / self.dt
		else:
			raise ValueError(f'Invalid method: "{method}". Choose "gradient" or "diff"')
		self.data = result

		return self

	#Function for calculating the Signal to Noise Ratio. The method is based on the simple SNR from scipy at version 0.4.0 (old version, not present in recent versions). Doing it with Power Spectral energies might be something in future. For now let's keep it simple.
	def SNR(self, dim='t'):
		'''
		Co-authors: --
		Description:
			under construction. DON'T USE THIS METHOD!
		:Params:
			- dim(type: String): dimension to where to apply the operation ('t' = time, 'd' = space). Default is 't'.
		:Return:
			- return1(type:--): --.
		'''

		axis = self.__axis__(dim)

		#m = self.data.mean(axis=0)
		sd = self.data.std(axis=axis)
		#result = np.where(sd == 0, 0, m/sd)
		#result = 20*np.log10(abs(result)) #For dB values
		#result = sd

		return sd


	#Function for plotting or returning the Root-Mean-Square amplitude (RMS-A) of the data.
	def rmsa(self, window=None, overlap=None, dim='t', make_plot=False):
		'''
		Co-authors: --
		Description:
			Calculates a RMS-Amplitude along the traces. Still under construction or need evaluation for approval.
		:Params:
			- window(type:Float): moving window length in seconds to use for the RMS-A calculation. Default is the time length of the data.
			- overlap(type:Float): overlapping time between windows. Still under construction. DO NOT USE.
			- dim(type: String): dimension to where to apply the operation ('t' = time, 'd' = space). Default is 't'.
			- make_plot(type: Bool): if set to True plot of the RMSA witll be generated
		:Return:
			- times(type:Numpy): array of the new times per each RMS-A value.
			- rms_a(type:Numpy): array containing the RMS-A values.
		'''

		axis = self.__axis__(dim)
		window = self.time_length if window == None else window
		times_d = np.array_split(self.times('matplotlib'), int(self.time_length/window))
		times = np.array([item[int(len(item)/2)] for item in times_d])
		data_d = np.array_split(self.data, int(self.time_length/window), axis=axis)
		rms_a = []

		for subdata in data_d:

			rms = np.sqrt(np.mean(subdata**2, axis=axis))
			rms_a.append(rms)

		if make_plot:
			plot.gen_DAS_plot(data=np.array(rms_a), t=times, channels=self.channels, cmap='inferno', title = f'RMSA for {window}s window')

		return times, np.array(rms_a)


	# Function for peak to peak amplitude calculation in every channel.
	def pp_amp(self, dim='t'):
		'''
		Co-authors: --
		Description:
			Calculate the peak to peak amplitude values per available channels.
		:Params:
			- dim(type: String): dimension to where to apply the operation ('t' = time, 'd' = space). Default is 't'.
		:Return:
			- ptp_amplitude(type: numpy): complete path of the resulting file.
		'''

		axis = self.__axis__(dim)
		ptp_amplitude = signals.peak_to_peak_amp(self.data, self.sampling_frequency, axis=axis)

		return ptp_amplitude



	'''
	####################################################
	Plotting methods below...
	Also proper from the Class
	####################################################
	'''


	#Function to plot spectrogram agains channels for an specific window defined by the actual length or start/end times of the DAS object. In order to avoid by computation time, please remember to trim first the DAS object to the time window of interest and/or restrict the number of channels before executing this function.
	def spectrogram(self, norm=False, max_value=None, order=1, nfft=None, figsize=None, show=True, cmap='viridis', results=False, file_name=None, where=None, plot_mode='pyqt', **kwargs):
		'''
		Co-authors: --
		Description:
			Plots the spectrogram of the DAS Class by channel instead of time dependent. The time window is defined by the actual start-time
			and end-time of the DAS Class. To avoid large computation times, better to use the trim() function of the class to first trim
			the data in time to the time-window of interest and then execute the spectrogram plot function.
		:Params:
			- norm(type:Boolean; optional): in case of True, each channel spectrum is normalized by its maximum value, so then the colorscale is not affected
			by the global maximum. Default = False.
			- order (type:Int): order number for detrending. Default is 1.
			- figsize(type:Tuple; optional): Tuple of 2 positions containing width and heigth of the figure. Default = None.
			- show(type:Boolean; optional): state if the plot must be shown. In case is False, the plot will not be shown, but the figure instance would be open
			so the user can add further changes. Default = True.
			- cmap(type:String; optional): name of the matplotlib colormap to use for the spectrogram. Default = 'viridis'.
			- results(type:Boolean): if set to True, the function will return the values for further manipulation (read Return section).
			Default = True.
   			- file_name(type:String; optional): in case the image want to be saved, this argument must be the name of the file, including the format
			(f.e.: "example.png"). Default = None.
			- where(type:String; optional): path of the directory where the plot wants to be saved.
		:Return:
			- NA.
		'''

		spectrogram = []

		for i in tqdm(range(self.total_channels), desc='Comptung Frequency Content', leave=False):

			o_signal = self.data[:,i] #- self.data[:,i].mean()
			#N = len(o_signal)

			#Pre-process
			#o_signal = tools.detrend_signal(o_signal) #detrend
			#o_signal -= o_signal.mean() #demean
			#o_signal *= signal.windows.tukey(M=self.num_points, alpha=0.05*2) #taper

			#powers = np.abs(rfft(o_signal)) #produce real and imaginary.

			#if i == 0:

			#	freqs = rfftfreq(N, 1 / self.sampling_frequency)

			#fft_values = np.flip(powers/powers.max()) if norm == True else np.flip(powers)
			#fft_values  = signal.savgol_filter(fft_values, 10, 2)

			freqs, fft_values = signals.spectrum(o_signal, self.sampling_frequency, True, order, int(self.num_points/4), nfft)
			fft_values = fft_values/fft_values.max() if norm == True else fft_values
			#fft_values  = signal.savgol_filter(fft_values, 10, 2) # to smooth the surve

			spectrogram.append(fft_values)
		spectrogram = np.array(spectrogram).T

		if results == True:
			return freqs, spectrogram

		else:
			if plot_mode == 'pyqt':
				plot_pyqt.plot_2d_distance(distances=self.distances,
                               y_ticks=freqs, data=np.flip(np.rot90(spectrogram, k=1), axis=0),
                               cmap=cmap, max_value=max_value,
                               y_label='Frequency [Hz]',
                               title='Frequency content over optical distance')

			elif plot_mode == 'mpl':
				plot.gen_spectrogram(spec_matrix=spectrogram[::-1], freqs=freqs,
                         x=self.channels_num, max_value=max_value, units_y=self.units,
                         figsize=figsize, title=self.start_time.isoformat()[:10],
                         cmap=cmap, show=show, file_name=file_name,
                         where=where, **kwargs)

	#Function to plot a spectrum (1D signal; freq vs Amplitude) of defined channel(s). Due to the label, it is recommended to not use many channels for plotting the spectrum, or can do it, but then legend must be turned off in the options (default = True).
	def spectrum(self, channels=None, norm=False, pre=True, order=1, pad=0, nfft=None, s_type='spectrum', figsize=None, show=True,
			  file_name=None, where=None, legend=True, results=False, **kwargs):
		'''
		Co-authors: --
		Description:
			Plots the spectrum of a specified channel, a list of them, or all the channels of the DAS Class. It is recommended not to use many channels
			with the legend option set as True.
		:Params:
			- channels(type:String or Int or Float or List; optional): channel to compute the spectrum.
			In case of a list, is all the channels specified in the list.
			In case is None, all spectrums of each channel would be computed. Default = None.
			- norm(type:Boolean; optional): in case of True, each channel spectrum is normalized by its maximum value. Default = False.
			- order (type:Int): order number for detrending. Default is 1.
			- pad(type:int): number of zeros to add to the signal before and after to increase num of points.
			- nfft(type:Int): number of samples in total Fast Fourier Transform.
			- s_type(type:String): Mode of spectral curve. 'spectrum' from normal spectral surve,
				'psd' for Power Spectral Density on Welch method. Default is 'spectrum'.
			- figsize(type:Tuple; optional): Tuple of 2 positions containing width and heigth of the figure. Default = None.
			- show(type:Boolean; optional): state if the plot must be shown. In case is False, the plot will not be shown,
			but the figure instance would be open so the user can add further changes. Default = True.
			- file_name(type:String; optional): in case the image want to be saved, this argument must be the name of the file, including the format
			(f.e.: "example.png"). Default = None.
			- where(type:String; optional): path of the directory where the plot wants to be saved.
			- legend(type:Boolean): sets if the legend would be shown or not. If many channels are being used, it's better to set this False.
			- results(type:Boolean): if set to True, the function will return the values for further manipulation (read Return section).
			Default = True.
		:Return:
			- spectrums(type:Numpy-Array): matrix spectral amplitude values of channels in order of the "channels" input variable.
			- freqs(type:Numpy-Array): frequencies used in the spectrum.
		'''

		spectrums = []

		# Evaluates how many chanels are comin as input. It can be one Int, or a list containing several Ints.
		if channels == None:

			channels = self.channels_num

		else:

			channels = [channels] if isinstance(channels,list) == False else channels

		for i in range(len(channels)):

			ch = int(channels[i])
			index = self.channels_num.index(ch)
			o_signal = self.data[:,index] #- self.data[:,index].mean()
			#N = len(o_signal)

			#Pre-process
			#o_signal = tools.detrend_signal(o_signal) #detrend
			#o_signal -= o_signal.mean() #demean
			#o_signal *= signal.windows.tukey(M=self.num_points, alpha=0.05*2) #taper

			#powers = np.abs(rfft(o_signal)) #produce real and imaginary.

			#if i == 0:

			#	freqs = rfftfreq(N, 1 / self.sampling_frequency)

			#fft_values = powers/powers.max() if norm == True else powers
			#fft_values = signal.savgol_filter(fft_values,15,2) smooth curve

			if s_type == 'spectrum':

				freqs, fft_values = signals.spectrum(o_signal, self.sampling_frequency, pre, order, pad, nfft)
				y_units = self.units

			if s_type == 'psd':

				freqs, fft_values = signals.psd(o_signal, self.sampling_frequency, pre, order, nfft)
				y_units = self.units.split(' ')[-1]
				y_units = y_units+'$^{2}$/Hz'

			fft_values = fft_values/fft_values.max() if norm == True else fft_values
			#fft_values = signal.savgol_filter(fft_values,15,2) smooth curve
			spectrums.append(fft_values)

		spectrums = np.array(spectrums)
		spectrums = spectrums[0] if spectrums.shape[0] == 1 and results == True else spectrums

		if results == True:
			return freqs, np.array(spectrums)
		else:
			plot.simple_spectrum(spectrums=np.array(spectrums), freqs=freqs, channels=channels, y_units=y_units, legend=legend, figsize=figsize,
						title=self.start_time.isoformat()[:10], show=show, file_name=file_name, where=where, **kwargs)


	def channel_plot(self, channel, max_value=None, figsize=None, show=True, file_name=None, where=None, plot_mode='pyqt', **kwargs):
		'''
		Co-authors: --
		Description:
			Plots the time-signal of a single selected channel.
		:Params:
			- channel(type:String or Int or Float): channel to plot.
			- max_value(type:Float; optional): maximum value of the y-axis. It will limit the plot in a range of -max_value to max_value. Default = None.
			- figsize(type:Tuple; optional): Tuple of 2 positions containing width and heigth of the figure. Default = None.
			- show(type:Boolean; optional): state if the plot must be shown. In case is False, the plot will not be shown,
			but the figure instance would be open so the user can add further changes. Default = True.
			- file_name(type:String; optional): in case the image want to be saved, this argument must be the name of the file, including the format
			(f.e.: "example.png"). Default = None.
			- where(type:String; optional): path of the directory where the plot wants to be saved.
		:Return:
			- NA.
		'''

		channel = int(channel)
		index = self.channels_num.index(channel)
		selected = self.data[:,index]

		if plot_mode=='pyqt':
			t = self.times(time_type='unix')
			plot_pyqt.plot_timeseries(data=selected, timestamps=t, y_label=self.units,
							 title='')

		elif plot_mode=='mpl':
			t = self.times('matplotlib')
			plot.simple_plot(data=selected, t=t, channel=str(channel),
					units_y=self.units, max_value=max_value, spectrogram=False,
					show=show, figsize=figsize,
					title=self.start_time.isoformat()[:10],
					file_name=file_name, where=where, **kwargs)


	def plot(self, max_value=None, figsize=None, show=True, cmap='seismic', file_name=None, where=None, add_data=None, plot_mode='pyqt', **kwargs):
		'''
		Co-authors: --
		Description:
			Plot the DAS Class data (matrix) as a colormap (Channel vs. Time).
		:Params:
			- max_value(type:Float; optional): maximum value of the colormap. It will limit the plot in a range of -max_value to max_value.
			All values above this will look saturated with the color limits of the colormap.
			- figsize(type:Tuple; optional): Tuple of 2 positions containing width and heigth of the figure. Default = None.
			- show(type:Boolean; optional): state if the plot must be shown. In case is False, the plot will not be shown,
			but the figure instance would be open so the user can add further changes. Default = True.
			- cmap(type:String; optional): name of the matplotlib colormap to use for the data. Default = 'seismic'.
			- file_name(type:String; optional): in case the image want to be saved, this argument must be the name of the file, including the format
			(f.e.: "example.png"). Default = None.
			- where(type:String; optional): path of the directory where the plot wants to be saved.
		:Return:
			- NA.
		'''
		if plot_mode == 'pyqt':
			t = self.times(time_type='unix')
			plot_pyqt.plot_2d_timeseries(timestamps=t, y_ticks=self.distances,
						data=self.data, y_label='Optical Distance [m]',
						title='', max_value=max_value, cbar_label=self.units)

		elif plot_mode == 'mpl':
			t = self.times(time_type='matplotlib')
			plot.gen_DAS_plot(data=self.data, t=t, channels=self.channels_num,
					 units_y=self.units, max_value=max_value, figsize=figsize,
					 show=show, title=self.start_time.isoformat()[:10], cmap=cmap,
					 file_name=file_name, where=where, add_data=add_data, **kwargs)

	def channel_spectrogram(self, channel, norm=False, trace=False, figsize=None, show=True, cmap='viridis', file_name=None, where=None, freq_lim=None, verbose=False, make_plot=True,
						 plot_mode='pyqt', **kwargs):
		'''
		Co-authors: --
		Description:
			Plots the spectrogram of a specific channel. This is different to the function spectrogram(), as this one is usng one channel,
			and its the spectrum shown in time, while spectrogram() shows it by space for a fixed time-window.
		:Params:
			- channel(type:String or Int or Float): channel to calculate the spectrogram.
			- norm(type:Boolean; optional): in case of True, each spectrogram window in time is normalized by its maximum value,
			so then the colorscale is not affected by the global maximum. Default = False
			- figsize(type:Tuple; optional): Tuple of 2 positions containing width and heigth of the figure. Default = None.
			- show(type:Boolean; optional): state if the plot must be shown. In case is False, the plot will not be shown,
			but the figure instance would be open so the user can add further changes. Default = True.
			- cmap(type:String; optional): name of the matplotlib colormap to use for the spectrogram. Default = 'viridis'.
			- file_name(type:String; optional): in case the image want to be saved, this argument must be the name of the file, including the format
			(f.e.: "example.png"). Default = None.
			- where(type:String; optional): path of the directory where the plot wants to be saved.
			- verbose(type:bool): if set to true, result (Spectrum, frequencies, time) is returned
		:Return:
			- NA.
		'''

		axis = self.__axis__('t')
		channel = int(channel)
		index = self.channels_num.index(channel)
		spec= self.data[:,index]

		nyquist = self.sampling_frequency/2
		#nfft, nperseg = nyquist*2, int(self.sampling_frequency/5) # TUNNING MUST BE DONE WITH PHILIPPE!
		nfft, nperseg = nyquist*2, int(self.sampling_frequency/5) # TUNNING MUST BE DONE WITH PHILIPPE!
		noverlap = int(nperseg/2)
		f, t, Sxx = signal.spectrogram(spec, self.sampling_frequency, nfft=nfft, nperseg=nperseg, noverlap=noverlap)
		Sxx = np.flip(Sxx,axis=axis)
		Sxx = Sxx / Sxx.max(axis=axis) if norm == True else Sxx
		trace = spec if trace == True else None

		if make_plot is True and plot_mode=='pyqt':
			t = self.times(time_type='unix')
			plot_pyqt.plot_2d_timeseries(timestamps=t, y_ticks=f,
						data=np.rot90(Sxx, k=-1), y_label='Frequency [Hz]',
						title=f'Spectrogram channel {channel}', cmap='viridis')

		elif make_plot is True and plot_mode=='mpl':
			t = self.times(time_type='matplotlib')
			plot.simple_spectrogram(data=Sxx, freq=f, t=t, units_y=self.units,
						trace=trace, figsize=figsize, cmap=cmap,
						title=self.start_time.isoformat()[:10]+'  '+'Ch:'+str(channel),
						show=show, file_name=file_name, where=where, freq_lim=freq_lim, **kwargs)

		if verbose is True:
			return Sxx, f, t
		#simple_spectrogram(spec_matrix=selected, freqs=self.sampling_frequency, x=t, units_x='time', figsize=figsize, cmap=cmap, title=self.start_time.isoformat()[:10], show=show, file_name=file_name, where=where, **kwargs)

	def plot_record_section(self, channels):
		'''
		Co-authors: --
		Description:
			Plot the DAS data as multiple seismograms in the same image (record section).
		:Params:
			- channels(type: tuple or list): channels to plot into the record section,
			if tuple all channels in the range will be plotted
			if list all channels given in the list will be plotted
		:Return:
			- NA.
		'''

		if isinstance(channels, tuple):
			ch0, chf = int(channels[0]), int(channels[1])
			ch0, chf = self.channels_num.index(ch0), self.channels_num.index(chf)
			das_data = self.data[:,ch0:chf+1]
			das_channels = self.channels_num[ch0:chf+1]

		elif isinstance(channels, list):
			channels = [int(channel) for channel in channels]
			ch_i = [self.channels_num.index(channel) for channel in channels]
			ch_i.sort()
			das_data = self.data[:, ch_i]
			das_channels = [self.channels_num[ind] for ind in ch_i]

		t = self.times('matplotlib')
		date = self.times()[0].isoformat()[:10]

		plot.plot_record_section(signals=das_data, t=t, channels=das_channels, date=date)


	def acf_profile(self, max_lag, plot_mode='pyqt', deconvolve=False,
                    window_size=None, result=False, **imshow_kwargs):
		'''
		computes autocorrelation profile, see fiber.tools.wavefield.autocorrelation_profile
        for more details

	    '''
		axis = self.__axis__('t')

		max_shift = int(max_lag*self.sampling_frequency)

		if max_shift >= self.num_points:
			raise ValueError('selected max_shift is too large')

		acf = wavefield.autocorrelation_profile(self.data, max_shift, axis, plot_mode,
                                                deconvolve, self.total_channels,
                                                self.distances, self.sampling_frequency,
                                                window_size=window_size, **imshow_kwargs)

		if result:
			return acf



	def spatial_coherence(self, max_lag, result=False, plot=True):
		'''
		computes sptial coherence matrix, see fiber.tools.wavefield.spatial_coherence_matrix
		for more details
	    '''

		coh = wavefield.spatial_coherence_matrix(data=self.data.T,
                                           max_lag=max_lag,
                                           fs=self.sampling_frequency,
                                           channel_nums=self.channels_num,
                                           plot=plot, result=result)
		if result:
			return coh


	def explore(self):
		'''
        launches the Fobench Data Explorer
		'''

		print(f'{"-"*65}\nStarting Fobench Data Explorer')
		app = QtWidgets.QApplication.instance()
		if app is None:
			app = QtWidgets.QApplication(sys.argv)

		self._explorer = Explorer(self)
		self._explorer.show()
		pg.exec()
		print(f'{"-"*65}')