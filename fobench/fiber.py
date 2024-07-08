"""
Class "Fiber" for creating, storing and manipulating fiber optic sensing data. 
So far it recieves TDMS format (Silixa) and H5 format (Febus).

Created on 2022-08-19 12:07:17
Last modification on 2024-06-28 19:17:00

:author:
	- Sergio Diaz-Meza (sergioad@gfz-potsdam.de)
:contributors:
	- Jonas Pätzel (jonas.patzel@ulb.be)
	- Christopher Wollin (wollin@gfz-potsdam.de)
:license:

"""




#Necessary packages to import
import numpy as np
import copy

import scipy.signal as signal
import scipy.integrate as integrate
from scipy.fft import rfft, rfftfreq, fftshift, ifftshift, fft2, ifft2

from obspy.core import UTCDateTime as UTC
from obspy.core.trace import Trace as oTrace
from obspy.core.stream import Stream

from pyrocko.util import str_to_time
from pyrocko.trace import Trace as pTrace

# inner functions
from . import read_data as read
from . import tools as tools
from . import filter as filter
from . import plotting as plot



class Fiber(object):
	'''
	IMPORTANT INFO: Most of the methods perform changes within the class permanently. Therefore is usefull to make a copy of the class
	with the method copy() before performing any processing or changes.
	'''
	
	#Creates the basic variables of the DAS object with its characteristics
	def __init__(self, filepath, company='silixa', range_ch=None, sensing='das', load_data=True):
		'''
		Co-authors: --
		Description: 
			Initializes a DAS Class which is reading a DAS file and saving all variables and metadata.
			So far it can read TDMS (Silixa) and H5 (Febus) formats.
			The basis for manipulating the data is numpy. Tools are inspired in Obspy, however using
			an obspy class for this takes long time in their processing tools.
		:Params:
			- filepath(type:String): compelte path of the file to be read.
			- company(type:String): manufacturer or the instrument that generates the data. Currently supporting "silixa" (Default), "febus", and "bam".
			- range_ch(type:Int or List): channel number(s) to load only in data. Method to avoid loading all the data.
			- sensing(type:String): specifies the type of fiber optic sensing technique of the data. Default is 'das'.
		:Return:
			- NA.  
		'''
	
		self.__filepath__ = filepath
		self.company = company
		self.format = filepath.split('.')[-1]
		# range_ch = [range_ch] if isinstance(range_ch,int) else range_ch
  
		attributes = read.read_data(self.__filepath__, self.company, range_ch, self.format, load_data=load_data)
		
		#file_file.close() # You might need more things from there afterwards.
		
		self.base = attributes.file
		self.fiber = attributes.fiber
		self.dataset = attributes.dataset
		self.properties = attributes.properties
		self.channels = attributes.chans
		self.channels_num = attributes.chans_nums
		self.total_channels = attributes.list_chans_num
		self.sampling_frequency = attributes.sampling_frequency # sampling rate of the data.
		self.dt = 1 / self.sampling_frequency # calculated time stamp.
		self.start_time = attributes.start_time # start time of the data in file.
		self.end_time = attributes.end_time # end time of the data in file.
		self.spatial_interval = attributes.spatial_interval # channel spaciong or spatial interval between channels [m].
		self.time_length = self.end_time - self.start_time
		self.num_points = attributes.num_points # int(self.time_length/self.dt)
		self.gauge_length = attributes.gauge_length # gauge length used in the measurement [m].
		self.channel_offset = attributes.channel_offset # offset where measurement started. It will not always record at channel 0 or distance 0.
		# self.data = self.__data__() if load_data == True else None
		self.data = attributes.data
		self.corrected = False
		self.sensing = sensing
		self.units = attributes.units
		self.conv_factor = attributes.conv_factor # Extra variables (ONLY FOR ASN HDF5)

		# Attributed not initialized since beginning. Required further processing to be initialized
		self.ch_coord = None # coordinates of channels.

		# Clean variables. Usually because h5py objects can't be copied with copy() function.
		del self.dataset
		del self.base


	'''
	####################################################
	Internal functions...
	####################################################
	'''


	def __str__(self):
		'''
		Co-authors: Jonas Pätzel
		Description: 
			defines output of print(Fiber); overview of most important recording parameters
		:Params:
			- NA.
		:Return:
			- NA.  
		'''
		return f'''Instance of Fiber class \n
recording parameters:
---------------------------------------------------------------------
{self.units = }
{self.start_time = }
{self.end_time = }
{self.total_channels = }
{self.spatial_interval = }
{self.sampling_frequency = }
{self.gauge_length = }
'''
	

	def __iadd__(self, other):
		'''
		Co-authors: Jonas Pätzel
		Description: 
			allows to concatenate two Fiber instances with self += other
		:Params:
			- other: Fiber instance.
		:Return:
			- NA.  
		'''
		if not isinstance(other, Fiber):
			raise TypeError
		return self.concatenate(other, fill_gaps=0)
	
	
	# #Loads the data of the tdms file into a numpy array. Axis 0 is the time, and axis 1 are the channels.
	# def __data__(self):
	# 	'''
	# 	Co-authors: --
	# 	Description: 
	# 		Extracts the data depending on the file type. This is done automatically during initialization of the class.
	# 	:Params:
	# 		- NA.
	# 	:Return:
	# 		- values(type:Numpy): 2D numpy matrix with values in time per channel. Axis 0 (rows) is time and axis 1 (columns) are the channels.  
	# 	'''
	
	# 	values = []
		
	# 	if self.format == 'tdms' and self.company == 'silixa':
			
	# 		values = np.array(self.channels).T # Old slow method.
	# 		# values = self.base['Measurement'].as_dataframe().to_numpy() # New way to load data. Cuts time by half.
			
	# 	if (self.format == 'h5' or self.format == 'hdf5') and self.company == 'febus':
		
	# 		dims = self.dataset.shape
	# 		values = self.dataset[:,:self.LAG,:].reshape(int(dims[0]*self.LAG),dims[2])
   
	# 	if (self.format == 'h5' or self.format == 'hdf5') and self.company == 'silixa':

	# 		values = np.array(self.dataset[:,self.channels])
			
	# 	if self.format == 'npy' and self.company == 'bam':
		
	# 		values = np.load(self.__filepath__)

	# 	if self.format == 'npz' and self.company == 'bam':
		
	# 		values = np.load(self.__filepath__)['ph']
   
	# 	if (self.format == 'h5' or self.format == 'hdf5') and self.company == 'terra15':

	# 		values = np.array(self.dataset['data'])

	# 	if (self.format == 'h5' or self.format == 'hdf5') and self.company == 'asn':

	# 		values = np.array(self.dataset)

	# 	if (self.format == 'h5' or self.format == 'hdf5') and self.company == 'quantx':

	# 		values = np.array(self.dataset['RawData'])
		
	# 	if (self.format == 'h5' or self.format == 'hdf5') and self.company == 'aragon':

	# 		values = np.array(self.dataset)*(10**-9)

	# 	return values.astype('float')


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
		
		axial = {'t':0, 'd':1}
    
		return axial[dim]


	'''
	####################################################
	Standard functions...
	####################################################
	'''


	# Returns the complete metadata of the file.
	def metadata(self):
		'''
		Co-authors: --
		Description: 
			Print out the metadata in an organized way.
		:Params:
			- NA.
		:Return:
			- NA.  
		'''
	
		for prop, value in self.properties.items():
			print(f"{prop} = {value}")
		
	
	#Return a deep copy of the object. Useful for instances where there is no wish to affect the original data while keeping notherone affected.	
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
		
			self.data, self.units, self.channels, self.channels_num, self.total_channels = tools.instr_corr(self.data, vars(self), target=target)
			
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
		t0 = self.start_time if t0 <= self.start_time else t0
		tf  = self.end_time if tf >= self.end_time else tf

		t = self.times()
		t0_new, tf_new, t0_pos, tf_pos = None, None, None, None
		i = 0

		for spec_time in t:

			if tf >= spec_time:
			
				tf_new = spec_time
				tf_pos = i
				
			if t0 >= spec_time:
			
				t0_new = spec_time
				t0_pos = i
			
			i += 1
								
		self.data = self.data[t0_pos:tf_pos,:]
		self.start_time = t0_new
		self.end_time = tf_new
		self.time_length = self.end_time - self.start_time
		self.num_points = self.data.shape[0]
		
		return self
			
	
	#Slice the data spatially (ranges of channels) by establshing the intial and the final channel. Anyformat of the channel code is acceptable.		
	def restrict_channels(self, ch):
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
		
		if isinstance(ch, tuple):
		
			ch0, chf = ch
		
		else:
		
			ch0 = ch
			chf = ch0
		
		ch0, chf = int(ch0), int(chf)
		ch0, chf = self.channels_num.index(ch0), self.channels_num.index(chf)
		self.data = self.data[:,ch0:chf+1]
		self.channels = self.channels[ch0:chf+1]
		self.channels_num = self.channels_num[ch0:chf+1]
		self.total_channels = len(self.channels_num)
		
		return self


	#Attach X and Y coordinates (generally lon and lat) to the Fibre class.		
	def append_coord(self, n_ch, x_ch, y_ch, z_ch):
		'''
		Co-authors: --
		Description:
			Attach coordinates of specified channels for the class instance. Necessary for plotting the fibre path and other spacial operations.
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

		n_ch, x_ch, y_ch, z_ch = tools.interpolate_channels(n_ch, x_ch, y_ch, z_ch, system, err, self.spatial_interval) # georeferencing new channels between the tap tests points.

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
			data_n = self.data[:,index]
			
		else:
		
			data_n = self.data
			
		return data_n 
		
	
	#Return a an array of timesin three different formats: UTCDateTime, ISOformat and matplotlib for plotting.	
	def times(self, time_type='UTCDateTime'):	
		'''
		Co-authors: --
		Description: 
			Return an array containing each time-step in the specified format option.
		:Params:
			- time_type(type:String): specific format of the time-steps. Options are: 1) UTCDateTime, 2) isoformat string, 3) matplotlib-dates (date-time)
			normally for plots. Default = 'UTCDateTime'
		:Return:
			- t(type:Numpy): a 1D array containing time-steps of the data in the specified format.  
		'''
		
		if time_type == 'UTCDateTime' or time_type == 'UTC':
		
			t = np.array([(self.start_time + (i*self.dt)) for i in range(self.data.shape[0])])
			
		elif time_type == 'isoformat':
		
			t = np.array([(self.start_time + (i*self.dt)).isoformat() for i in range(self.data.shape[0])])
			
		elif time_type == 'matplotlib':
		
			t = np.array([(self.start_time + (i*self.dt)).matplotlib_date for i in range(self.data.shape[0])])
			
		else:
		
			raise ValueError('Unrecognized time format. Please check the possible values.')
		
		return t


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
		
		if num_t < 0:
			
			num_t = abs(num_t)
			second.data = second.data[num_t:,:]
		
		if num_t > 0:
		
			fill = np.zeros((num_t, first.total_channels))
			fill[fill==0] = np.nan if fill_gaps == None else fill_gaps #Can also work for putting NonType values (NaN) if fill_gaps is None or any value.
			first.data = np.concatenate((first.data, fill), axis=axis)
			
		self.data = np.concatenate((first.data, second.data), axis=axis)
		self.start_time = first.start_time
		self.end_time = second.end_time
		self.num_points += second.num_points
		self.time_length = self.end_time - self.start_time
			
		'''
		#Old version of Code
		if self.start_time <= input_das.start_time:
		
			tf = self.end_time + self.dt
			num_t = int((input_das.start_time+input_das.dt - self.end_time)/self.dt)-1
			fill = np.zeros((num_t, self.total_channels))
			
			if fill.size != 0:
				
				self.data = np.concatenate((self.data,fill), axis=0)

			self.data = np.concatenate((self.data, input_das.data), axis=0)

			self.end_time = input_das.end_time
			self.num_points += input_das.num_points
			self.time_length = self.end_time - self.start_time
		'''
			
		return self


	#Creates a Stream object from obspy with traces inside and returns it for obspy-type manipulation. Each Trace class represents a channel from the DAS data. All Traces must have the same time range and number of points.
	def as_Traces(self, t_type='obspy'):
		'''
		Co-authors: --
		Description:
			Creates an obpsy/pyrocko Stream object and fill it with Traces in it. Each Trace would represent each channel of the DAS Class, including 
			the metadata which are attributes of the Trace Class. This is mainly done so users can have access to obspy tools with this data. However,
			it can be slower and memory demanding.
		:Params:
			- t_type(type:String): option wether to convert to pyrocko or obspy stream/traces.
		:Return:
			- stream(type:Stream Class): stream with traces representing channels of the DAS Class..  
		'''

		stream = Stream() if t_type == 'obpsy' else []

		for i in range(self.total_channels):
			
			if t_type == 'obspy':
   
				trace = oTrace(data=self.data[:,i])
				trace.stats.network = self.fiber
				trace.stats.station = str(self.channels_num[i]).zfill(5)
				trace.stats.npts = self.num_points + 1
				trace.stats.sampling_rate = self.sampling_frequency
				trace.stats.delta = self.dt
				trace.stats.starttime = self.start_time
				trace.stats.calib = tools.instr_corr(np.array(1), d_type='das')
				#trace.stats.endtime = self.end_time

				stream += trace
    
			if t_type == 'pyrocko':

				trace = pTrace(ydata=self.data[:,i])
				trace.network = self.fiber
				trace.station = str(self.channels_num[i]).zfill(5)
				trace.deltat = self.dt
				trace.tmin = str_to_time(self.start_time.isoformat().replace('T',' '))
				trace.tmax = str_to_time(self.end_time.isoformat().replace('T',' '))
				stream.append(trace)

		return stream
		

	'''
	####################################################
	Signal Processing functions...
	####################################################
	'''

		
	#function for upsampling spatially by double/half depending if it is upsampling or downsampling spatially.
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
			new_data, new_channels_num = tools.spatial_upsampling(self)
			self.spatial_interval = self.spatial_interval / 2
			
		elif rs_type == 'downsampling':
		
			new_data, new_channels_num = tools.spatial_downsampling(self)
			self.spatial_interval = self.spatial_interval * 2
		
		#if (rs_type != 'downsampling') & (rs_type != 'upsampling'):
		else:
			
			raise ValueError('Spatial resampling type is not recognizable. Only (upsampling) or (downsampling).')
		
		self.data = new_data
		self.channels_num = new_channels_num
		self.total_channels = len(self.channels_num)
		
		return self


	#Function for detrending the data
	def detrend(self, order=1, dim='t'):
		'''
		Co-authors: --
		Description:
			Detrends the data, taking any unwanted trend component in the data that might come artifacts such as temperature, instrument, very long period signal, etc.
		:Params:
			- order(type: Int): order number of the fitting curve used to apply the detrend.
			- dim(type: String): dimension to where to apply the operation ('t' = time, 'd' = space). Default is 't'.
		:Return:
			- NA.  
		'''

		axis = self.__axis__(dim)		

		M = self.total_channels if axis == 0 else self.num_points

		for i in range(M): # I think there is a way to do this matrix wise, and coefficients might appear per column. See numpy.polyfit()
    
			if axis == 0: # if it's time sample
    
				self.data[:,i] = tools.detrend_signal(self.data[:,i], order)
    
			if axis == 1: # if it's space sample
    
				self.data[i,:] = tools.detrend_signal(self.data[i,:], order)
		
		return self


	#Function for demeaning the data
	def demean(self, dim='t'):
		'''
		Co-authors: --
		Description:
			Demean the data, trying to reduce any trend outside of the 0 value line.
		:Params:
			- dim(type: String): dimension to where to apply the operation ('t' = time, 'd' = space). Default is 't'.
		:Return:
			- NA.  
		'''
		
		axis = self.__axis__(dim)
		self.data -= self.data.mean(axis=axis)
		
		return self
	
	#Function for tappering the data
	def taper(self, frac=0.05, dim='t'):
		'''
		Co-authors: --
		Description:
			Tapers the data in time (dim='t') or in space (dim='d'). The taper used is a tapered cosine window (Tukey).
		:Params:
			- frac(type:Float): it is the fraction of the taper applied to one side of the window. In total the tapered part of the data will be twice of the indicated in the parameter.
			- dim(type: String): dimension to where to apply the operation ('t' = time, 'd' = space). Default is 't'.
		:Return:
			- NA.  
		'''
		
		axis = self.__axis__(dim)
		M = self.num_points if axis == 0 else self.total_channels
		taper = signal.windows.tukey(M=M, alpha=frac*2)
		taper = taper[:,None] if axis == 0 else taper[None, :]
		
		self.data = np.multiply(self.data, taper)
		
		return self
		
	
	#Function to decimate the data by any frequency below the original. Carefull when applying decimations with factors over or equal to 13, then better to call decimation twice (see scipy.signal.decimate for more...). The decimation function of scipy performs a pre-filtering process to avoid anti-aliasing on the signals.
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
		down_factor = int(self.sampling_frequency / new_freq)
		new_freq = self.sampling_frequency / down_factor
		
		#Check prefilter... which one is?
		if ftype is not None:
			
			new_data = filter.decimate(data=self.data, factor=down_factor, ftype=ftype, axis=axis)
			
		else:
		
			new_data = signal.decimate(x=self.data, q=down_factor, axis=axis)
		
		self.data = new_data
		self.sampling_frequency = new_freq
		self.dt = 1 / self.sampling_frequency
		self.num_points = self.data.shape[0]

		return self


	# does it need to be also generalized for dimension option? (f.e.: if i want to do it in time or spatial?)
	def normalize(self, method='absolute max', dim='d'):
		'''
		Co-authors:Jonas Pätzel
		Description:
			min-max normalizes the data according to the chosen method
		:Params:
			- method (type:String): normalization either with respect to 
			the whole record ('absolute max', default) or for each channel/timestep individually ('trace max').
			- dim(type: String): dimension to where to apply the operation ('t' = time, 'd' = space). Default is 'd'.
		:Return:
			- NA.  
		'''

		axis = self.__axis__(dim)

		if method == 'absolute max':
      
			normalized_data = (self.data - self.data.min()) / (self.data.max() - self.data.min())
			normalized_data = normalized_data * 2 -1 
   
		elif method == 'trace max':
      
			channel_min = self.data.min(axis=axis, keepdims=True)
			channel_max = self.data.max(axis=axis, keepdims=True)
			normalized_data = (self.data - channel_min)/(channel_max - channel_min)
			normalized_data = normalized_data * 2 -1
		
		self.data = normalized_data
		
		return self


	# Function for filtering. Shall we also declare dimensionality options here?
	def filter(self, f_type=None, freq=None, pre_process=True, frac=0.05, order=1, **options):
		'''
		Co-authors: --
		Description:
			Filters the data based on a specified type of filter (lowpass, bandpass, highpass) and the values. Filters are based on Obspy codes 
			and so the multiple options are also.
		:Params:
			- f_type(type:String): type of filter to apply. Options are: 'lowpass', 'bandpass', 'highpass'.
			- freq(type:Int, Float, Tuple): cut-off value for the filter in 'lowpass' and 'highpass'. If it for 'bandpass', then it must be a tuple containing the cut-offs of the bandwith.
			- pre_process(type:Boolean): to use to detren, demean and tape before filtering. Default is True.
			- frac(type:Float): see description in method "taper".
			- order(type:Int): order of polynomuial fit for detrending.
		:Return:
			- NA.  
		'''

		if pre_process == True:		
  
			self.detrend(order=order)
			self.demean()
			self.taper(frac=frac)
    
		new_data = filter.point_filter(f_type=f_type, data=self.data, df=self.sampling_frequency, freq=freq, **options)
		self.data = new_data
		
		return self


	# Function to make a fk-filter based on the input parameters.
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
		
		if taper == True:
			self.taper(axis=axis)
		
		if method == 'cum_trapezoid':
			res = integrate.cumulative_trapezoid(y=self.data, dx=dx, axis=axis, initial=0) #+ self.data[0,:]
		
		res = signal.detrend(res, axis=axis) #to detrend the signal
		
		#if taper == True:
		
		#	self.detaper(axis=axis)
		
		self.data = res
		self.units = 'Strain'
		
		return self
	
	
	#Function for differentiating the signal
	def differentiate(self, method='gradient', dim='t'):
		'''
		Co-authors: --
		Description:
			Differentiates the data in time (dim='t') or in space (dim='d').
		:Params:
			- method(type:String): sets the prefered method for differentiation. Default and only one is "gradient".
			- dim(type: String): dimension to where to apply the operation ('t' = time, 'd' = space). Default is 't'.
		:Return:
			- NA.  
		'''

		axis = self.__axis__(dim)
		res = np.gradient(self.data, self.dt, axis=axis)
		self.data = res
		self.units = 'Strain acc. [1/s$^{2}$]'
		
		return self
		
		
	#Function to detaper...
	def detaper(self, frac=0.05, dim='t'):
		'''
		Co-authors: --
		Description:
			Tapers the data in time (dim='t') or in space (dim='d'). The taper used is a tapered cosine window (Tukey).
		:Params:
			- frac(type:Float): it is the fraction of the taper applied to one side of the window. In total the tapered part of the data will be twice of the indicated in the parameter.
			- dim(type: String): dimension to where to apply the operation ('t' = time, 'd' = space). Default is 't'.
		:Return:
			- NA.  
		'''
		
		axis = self.__axis__(dim)
		M = self.num_points if axis == 0 else self.total_channels
		taper = signal.windows.tukey(M=M, alpha=frac*2)
		taper = taper[:,None] if axis == 0 else taper[None, :]
		
		self.data = np.divide(self.data, taper)
		
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
	def rmsa(self, window=None, overlap=None, dim='t'):
		'''
		Co-authors: --
		Description:
			Calculates a RMS-Amplitude along the traces. Still under construction or need evaluation for approval.
		:Params:
			- window(type:Float): moving window length in seconds to use for the RMS-A calculation. Default is the time length of the data.
			- overlap(type:Float): overlapping time between windows. Still under construction. DO NOT USE.
			- dim(type: String): dimension to where to apply the operation ('t' = time, 'd' = space). Default is 't'.
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
		ptp_amplitude = tools.peak_to_peak_amp(self.data, self.sampling_frequency, axis=axis)

		return ptp_amplitude
		
		
		
	'''
	####################################################
	Plotting functions below...
	Also proper from the Class
	####################################################
	'''
	
	
	#Function to plot spectrogram agains channels for an specific window defined by the actual length or start/end times of the DAS object. In order to avoid by computation time, please remember to trim first the DAS object to the time window of interest and/or restrict the number of channels before executing this function. 
	def spectrogram(self, norm=False, max_value=None, order=1, nfft=None, figsize=None, show=True, cmap='viridis', file_name=None, where=None, **kwargs):
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
			- file_name(type:String; optional): in case the image want to be saved, this argument must be the name of the file, including the format 
			(f.e.: "example.png"). Default = None.
			- where(type:String; optional): path of the directory where the plot wants to be saved.
		:Return:
			- NA.  
		'''
	
		spectrogram = []
	
		for i in range(self.total_channels):
			
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

			freqs, fft_values = tools.spectrum(o_signal, self.sampling_frequency, True, order, int(self.num_points/4), nfft)
			fft_values = np.flip(fft_values/fft_values.max()) if norm == True else np.flip(fft_values)
			#fft_values  = signal.savgol_filter(fft_values, 10, 2) # to smooth the surve


			spectrogram.append(fft_values)
			
		spectrogram = np.array(spectrogram).T
			
		plot.gen_spectrogram(spec_matrix=spectrogram, freqs=freqs, x=self.channels_num, max_value=max_value, units_y=self.units, figsize=figsize, title=self.start_time.isoformat()[:10], cmap=cmap, show=show, file_name=file_name, where=where, **kwargs)
		
		
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
       
				freqs, fft_values = tools.spectrum(o_signal, self.sampling_frequency, pre, order, pad, nfft)
				y_units = self.units
    
			if s_type == 'psd':
       
				freqs, fft_values = tools.psd(o_signal, self.sampling_frequency, pre, order, nfft)
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



	#Function for plotting a 1D time series of one specific channel.	
	def channel_plot(self, channel, max_value=None, figsize=None, show=True, file_name=None, where=None, **kwargs):
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
		#print(self.channels_num, type(self.channels_num))
		#index = np.where(self.channels_num == channel)[0]
		print(index)
		selected = self.data[:,index]
		t = self.times('matplotlib')
		
		plot.simple_plot(data=selected, t=t, channel=str(channel), units_y=self.units, max_value=max_value, spectrogram=False, show=show, figsize=figsize, title=self.start_time.isoformat()[:10], file_name=file_name, where=where, **kwargs)
		
	
	#Fast plotting function of the data as matrix. Maximum value can be adjusted to saturate the plot.	
	def plot(self, max_value=None, figsize=None, show=True, cmap='seismic', file_name=None, where=None, add_data=None, **kwargs):
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
	
		t = self.times(time_type='matplotlib')
		
		plot.gen_DAS_plot(data=self.data, t=t, channels=self.channels_num, units_y=self.units, max_value=max_value, figsize=figsize, show=show, title=self.start_time.isoformat()[:10], cmap=cmap, file_name=file_name, where=where, add_data=add_data, **kwargs)
		
		
	#Function for plotting single channel spectrogram
	def channel_spectrogram(self, channel, norm=False, trace=False, figsize=None, show=True, cmap='viridis', file_name=None, where=None, freq_lim=None, **kwargs):
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
		t = self.times(time_type='matplotlib')
		Sxx = np.flip(Sxx,axis=axis)
		Sxx = Sxx / Sxx.max(axis=axis) if norm == True else Sxx
		trace = spec if trace == True else None
		
		plot.simple_spectrogram(data=Sxx, freq=f, t=t, units_y=self.units, trace=trace, figsize=figsize, cmap=cmap, title=self.start_time.isoformat()[:10]+'  '+'Ch:'+str(channel), show=show, file_name=file_name, where=where, freq_lim=freq_lim, **kwargs)
		#simple_spectrogram(spec_matrix=selected, freqs=self.sampling_frequency, x=t, units_x='time', figsize=figsize, cmap=cmap, title=self.start_time.isoformat()[:10], show=show, file_name=file_name, where=where, **kwargs)
		
		
	def interactive_plot(self, channel=None, max_value=None, figsize=None, show=True, cmap='seismic', file_name=None, where=None, add_data=None, **kwargs):
		'''
		Co-authors: --
		Description: (UNDER CONSTRUCTION BUT WORKS SO FAR)
			Interactive plot based on plot() and channel_plot() functions. It plots the DAS Class matrix data as a colormap, and below it plots 
			a single-channel signal. The signal channel is shown in the matrix colormap plot as a yellow line, indicating the positon of the 
			current plotted channel. A box in the upper-right part of the channel plot allows the user to change the channel to visulize by
			enterirng the number of a new channel, and by hitting ENTER, the plot updates showing the new channel plot and indicating where on the
			colormap plot the channel is located.
		:Params:
			- channel(type:String or Int or Float): channel to plot initially.
			- max_value(type:Float; optional): maximum value of the colormap and y-axis. It will limit the plot in a range of -max_value to max_value.
			All values above this will look saturated with the color limits of the colormap. Default = None.
			- figsize(type:Tuple; optional): Tuple of 2 positions containing width and heigth of the figure. Default = None.
			- show(type:Boolean; optional): state if the plot must be shown. In case is False, the plot will not be shown, 
			but the figure instance would be open so the user can add further changes. Default = True.
			- cmap(type:String; optional): name of the matplotlib colormap to use for the data. Default = 'viridis'.
			- file_name(type:String; optional): in case the image want to be saved, this argument must be the name of the file, including the format 
			(f.e.: "example.png"). Default = None.
			- where(type:String; optional): path of the directory where the plot wants to be saved.
		:Return:
			- return1(type:--): --.  
		'''

		channel = 0 if channel == None else channel
		t = self.times(time_type='matplotlib')
		
		plot.DAS_interactive_plot(data=self.data, set_channel=channel, t=t, channels=self.channels_num, units_y=self.units, max_value=max_value, figsize=figsize, title=self.start_time.isoformat()[:10], cmap=cmap, show=show, file_name=file_name, where=where, add_data=add_data, **kwargs)
  

	#Function for plotting the DAS data as record structure/section.
	def plot_record_section(self, channels):
		'''
		Co-authors: --
		Description:
			Plot the DAS data as multiple seismograms in the same image (record section).
		:Params:
			- channels(type: list): List of the channels to plot into the record section. If not, all of them will be plotted as Default.
		:Return:
			- NA.
		'''

		ch0, chf = int(channels[0]), int(channels[1])
		ch0, chf = self.channels_num.index(ch0), self.channels_num.index(chf)
		das_data = self.data[:,ch0:chf+1]
		das_channels = self.channels_num[ch0:chf+1]
		t = self.times('matplotlib')
		date = self.times()[0].isoformat()[:10]

		plot.plot_record_section(signals=das_data, t=t, channels=das_channels, date=date)
