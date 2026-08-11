"""
Class "Fiber" for creating, storing and manipulating fiber optic sensing data.
Created on 2022-08-19 12:07:17

:authors:
	- Sergio Diaz-Meza (sergioad@gfz-potsdam.de)
	- Jonas Pätzel (jonas.patzel@ulb.be)

:contributors:
	- Christopher Wollin (wollin@gfz-potsdam.de)
:license:

"""

import sys
import copy
from warnings import warn

import numpy as np
from obspy.core import UTCDateTime as UTC

import pyqtgraph as pg
from pyqtgraph.Qt import QtWidgets

from .tools import file_io, utils, filters, signals, wavefield
from .plotting import plotting_mpl as plot
from .plotting import plotting_pyqt as plot_pyqt
from .plotting.pyqt_viewer import Viewer


class Fiber(object):
	'''
	NOTE: Most of the methods perform changes within the class permanently.
	Therefore it is useful to make a copy of the fiber instance with the method
	.copy() before performing any processing or changes.

	For most methods 'plot_mode' parameter can either be 'pyqt' or 'mpl',
	to not generate plot set to anything else, e.g. None
	'''

	def __init__(self, filepath, company="", range_ch=None, sensing="das",
	             load_data=True, show_progress=True, storage_opts=None):
		'''
		Initializes base class of Fobench, reading in data and metadata, data is manipulated mostly
		using numpy and scipy, tools are inspired by obspy
		needs 'filepath' to data file to read and manufacturer of interrogator, currently supported companies are:
		'silixa', 'febus', 'aragon', 'quantx', 'asn', 'terra15' and 'sintela'
		'sensing' indicates type of fiber optic sensing technology, defaulting to 'das'
		only part of data can be loaded when specifing the channel range; if 'load_data' is False,
		class is initialized containing only metadata
		'''
		if not company:
			raise ValueError(
				'\nNo company provided! Please choose one of:\n'
				' -"silixa"\n -"febus"\n -"bam"\n -"aragon"\n -"quantx"\n -"asn"\n'
				' -"terra15"\n -"sintela"'
			)

		self.__filepath__ = [filepath]

		self.company = company
		self.format = filepath.split('.')[-1]

		self.attributes = file_io.read_data(self.__filepath__[0], self.company, range_ch,
						self.format, load_data=load_data, show_progress=show_progress, storage_opts=storage_opts)

		self.__basefile__ = self.attributes['basefile'] # changed to the structure of the file
		self.fiber = self.attributes['fiber']
		self.properties = self.attributes['properties'] # all metadata of input file
		self.channels = self.attributes['chans'] # list of channels as array
		self.channels_num = self.attributes['chans_nums'] # list of channels as list
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
		self.distances = [(num + self.channel_offset) * self.spatial_interval for num in self.channels_num]

		self.ch_coord = None # coordinates of channels, requires more input ot be filled

	'''
	-----------------------------------------------------------------
	Internal methods
	-----------------------------------------------------------------
	'''


	def __str__(self):
		attributes = ['units', 'start_time', 'end_time', 'num_points', 'total_channels',
					'spatial_interval', 'sampling_frequency', 'gauge_length']

		return ('Instance of Fiber class\n'
				'recording parameters:\n'
				f'{"-" * 65}\n'+ "\n".join(f"{attr.ljust(25)} = {getattr(self, attr)}" for attr in attributes))

	__repr__ = __str__

	def __iadd__(self, other):
		if not isinstance(other, Fiber):
			raise TypeError('Object to add must be instance of Fiber class')

		return self.concatenate(other, fill_gaps=0)

	def __axis__(self, dim):
		'''
		translates string to numerical axis value for numpy arrays, dim can be 't' (time)
		or 'd' (distance)
		axis 0 corresponds to rows, axis 1 to columns
		'''
		axial = {'t':0, 'd':1}

		return axial[dim]

	'''
	-----------------------------------------------------------------
	Standard methods
	-----------------------------------------------------------------
	'''

	def metadata(self, meta_dict=False):
		'''
		print out metadata, optionally return all metadata as dictionary
		'''

		if meta_dict:
			metainfo = {key: value for key, value in vars(self).items() if not key.startswith('__')}
			return metainfo

		else:
			for prop, value in self.properties.items():
				print(f"{prop} = {value}")

	def copy(self):
		'''
		returns a deep copy the Fiber class object in its current state
		'''
		self.__dict__.pop('_viewer', None)
		return copy.deepcopy(self)

	def instr_correct(self, target='strain-rate', terra15_gl=None):
		'''
		performs instrument correction and data conversion for various instrument types,
		exact conversion depends on manufacturer and data format
		'''
		if not self.corrected:
			(self.data, self.units, self.channels, self.channels_num,
					self.total_channels, self.gauge_length, self.distances) = utils.instr_corr(self.data, vars(self),
									target=target, terra15_gl=terra15_gl)
			# self.distances = [(num + self.channel_offset) * self.spatial_interval for num in self.channels_num]
			self.corrected = True
			return self

		elif self.corrected:
			warn('Instrument correction has already been applied, doing nothing...')
			return self


	def trim(self, t0=None, tf=None):
		'''
		cuts data between given start and end times, t0 and tf can be UTC datetime
		or ISOformat style str
		'''

		data, start_time, end_time = utils.trim_time(t0=t0, tf=tf, data=self.data,
													 times=self.times(), start_time=self.start_time,
													 end_time=self.end_time)

		self.data, self.start_time, self.end_time, self.time_length, self.num_points = (data,
				start_time, end_time, end_time-start_time, data.shape[0])

		return self

	def restrict_channels(self, ch0, chf):
		'''
		trims data in space, between ch0 and chf, a single channel is returned
		when ch0 = chf, updates Fiber class attributes
		'''

		ch0, chf = int(min(ch0, chf)), int(max(ch0, chf)) # in case ch0 and chf not ordered
		ch0, chf = self.channels_num.index(ch0), self.channels_num.index(chf)
		self.data = self.data[:,ch0:chf+1]
		self.channels = self.channels[ch0:chf+1]
		self.channels_num = self.channels_num[ch0:chf+1]
		self.distances = self.distances[ch0:chf+1]
		self.total_channels = len(self.channels_num)

		return self

	def append_coord(self, n_ch, x_ch, y_ch, z_ch):
		'''
		attaches channel coordinates for later plotting
		takes 1D arrays of channel number (n_ch), longitude and latitude (x_ch and y_ch)
		and elevation in m (z_ch)
		'''
		coords = [n_ch,
				  np.zeros_like(n_ch) if x_ch is None else x_ch,
				  np.zeros_like(n_ch) if y_ch is None else y_ch,
				  np.zeros_like(n_ch) if z_ch is None else z_ch]

		self.ch_coord = np.column_stack(coords)

		return self

	def georeference(self, n_ch, x_ch, y_ch, z_ch, system='decimal', err=None):
		'''
		takes known channel locations, e.g. from tap tests and interpolates channel locations
		inbetween, attaches new coordinates
		takes 1D arrays of channel number (n_ch), longitude and latitude (x_ch and y_ch)
		and elevation in m (z_ch), coordinate system can be for lon and lat can be 'decimal' or 'utm'
		'err' is maximum accepted interpolation error between original metadata location and new interpolated
		location
		'''
		x_ch = np.zeros(n_ch.size) if x_ch is None else x_ch
		y_ch = np.zeros(n_ch.size) if y_ch is None else y_ch
		z_ch = np.zeros(n_ch.size) if z_ch is None else z_ch

		n_ch, x_ch, y_ch, z_ch = utils.interpolate_channels(n_ch, x_ch, y_ch,
													  z_ch, system, err, self.spatial_interval)
		self.append_coord(n_ch, x_ch, y_ch, z_ch)

		return self

	def get_data(self, channel=None):
		'''
		returns data similar to Fiber.data but has option to return only a specified channel
		'''
		if channel is not None:
			index = self.channels_num.index(int(channel))
			return self.data[:,index]

		return self.data

	def times(self, time_type='UTCDateTime'):
		'''
		returns array of sample times, can be 'UTCDateTime', 'isoformat', 'datetime64', 'matplotlib'
		or 'unix'
		'''
		return utils.return_times(self, time_type)

	def concatenate(self, input_das=None, fill_gaps=0):
		'''
		concats two Fiber class objects assuming they have the same sampling rate and channels
		order is determined automatically, gaps are filled with value of 'fill_gaps', overlapping times are taken
		from self, the rest filled with input_das
		updates class properties
		'''

		axis = self.__axis__('t')

		first, second = (self, input_das) if self.start_time <= input_das.start_time else (input_das, self)

		self.data = np.concatenate((first.data, second.data), axis=axis)
		self.start_time, self.end_time = first.start_time, second.end_time
		self.time_length = self.end_time - self.start_time
		self.num_points = self.data.shape[axis]
		self.__filepath__.extend(input_das.__filepath__)

		return self

	def to_traces(self, t_type='obspy'):
		'''
		returns channels as 'obspy' or 'pyrocko' streams, see
		fobench.tools.utils.to_traces for more details
		'''
		return utils.to_traces(self, t_type)

	def to_xarray(self, name=None, use_distance=True):
		"""Converts the Fiber class and data as xarray object.

		Returns:
			DataArray: DataArray of xarray containing the Fiber data and metadata.
						Can be useful for Xdas 
		"""
  
		return utils.to_xarray(self, name, use_distance)


	def write(self, save_path=None):
		"""
		Save the data of fiber in a new file with the original format from where it was read.

		Parameters
		----------
		save_path : str, optional
			Path to where to save the file including name of the file and the format.
		"""

		file_io.write_data(self, filepath=save_path, company=self.company)

	'''
	-----------------------------------------------------------------
	Signal Processing methods
	-----------------------------------------------------------------
	'''

	@utils._update_processing
	def spatial_resample(self, rs_type=None):
		'''
		modifies spatial sampling of the data by adding or removing channels,
		for 'upsampling' adds a channel between each channel pair by interpolating the values,
		for 'downsampling' removes every second channel, first channel is always unaltered
		see fobench.tools.utils.spatial_upsampling and .spatial_downsample for details
		'''
		if rs_type in ['upsampling', 'upsample']:
			self.data, self.channels_num = utils.spatial_upsampling(self)
			self.spatial_interval /= 2
		elif rs_type in ['downsampling', 'downsample']:
			self.data, self.channels_num = utils.spatial_downsampling(self)
			self.spatial_interval *= 2
		else:
			raise ValueError(f'\nInvalid resample type: "{rs_type}". Choose on of:\n'
					' -"upsampling"\n -"downsampling"')
		self.total_channels = len(self.channels_num)

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

	@utils._update_processing
	def decimate(self, new_freq=None, f_type='fir-remez'):
		'''
		decimates data to new sampling frequeny, target frequency should divide
		original sampling frequency evenly

		! careful when decimating using factors >= 13, it is then preferable to
		call decimation twice (see scipy.signal.decimate for more info) !

		options for filters are
			- 'fir-remez', an adaptative antialiasing filter (author: Marius Isken)
			- 'fir235' (author: Javier Quinteros)
			- 'None', scipy's default anti-aliasing order 8 Chebyshev Type I filter
		'''
		axis = self.__axis__('t')

		if new_freq is None:
			raise ValueError('new_freq must be provided as a positive number in Hz')
		if new_freq <= 0:
			raise ValueError(f'new_freq must be > 0 Hz, got {new_freq}')
		if new_freq > self.sampling_frequency:
			raise ValueError(f'new_freq ({new_freq} Hz) cannot exceed current sampling frequency ({self.sampling_frequency} Hz)')

		if self.sampling_frequency % new_freq != 0:
			warn(f'Decimation to {new_freq} Hz not possible! Decimating to {self.sampling_frequency / int(self.sampling_frequency / new_freq)} Hz instead')
		down_factor = int(self.sampling_frequency / new_freq)
		new_freq = self.sampling_frequency / down_factor

		self.data = filters.decimate(data=self.data, factor=down_factor, f_type=f_type, axis=axis)

		self.sampling_frequency  = new_freq
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

		self.data = signals.whiten_signal(data = self.data, freq_min=freq_min, freq_max=freq_max,
									sampling_frequency=self.sampling_frequency,
									total_channels=self.total_channels)

		return self

	@utils._update_processing
	def filter(self, f_type=None, freq=None, pre_process=True, alpha=0.05, order=1, sym=True,
			**options):
		'''
		filters data using specified filter, based on Obspy.signal.filter module
		for frequency filters, data is optionally pre-processed
		filter types are bandpass', 'bandstop', 'lowpass', 'highpass' and 'median''
		if 'bandpass' or 'bandstop', freq must be tuple(float, float)
		'''
		if pre_process and f_type != 'median':
			self.preprocess(alpha=alpha, sym=sym, order=order, axis=0)

		self.data = filters.point_filter(f_type=f_type, data=self.data,
								  df=self.sampling_frequency, freq=freq, **options)

		return self

	@utils._update_processing
	def preprocess(self, alpha=0.05, order=1, sym=True, axis=0, steps=(True, True, True)):
		'''
		performs demeaning, detrending and tapering, see fobench.tools.signals.filt_preprocess
		for more details
		'''
		self.data = signals.filt_preprocess(io_signal=self.data, order=order,
									  alpha=alpha, sym=sym, axis=axis)
		return self

	@utils._update_processing
	def fk_filter(self, bands=[{}], propagation='both', alpha=0.3, plot_mode=None,
				  verbose=False, results=False, mode='pass'):
		'''Applies frequency wavenumber filter to data'''

		out = filters.fk_filter(data=self.data, dt=self.dt, dx=self.spatial_interval,
								bands=bands, propagation=propagation, alpha=alpha,
								plot_mode=plot_mode, verbose=verbose, mode=mode)
		self.data = out[0] if len(out) == 3 else out
		if results:
			return (out[0], out[1], out[2]) if verbose else (out[0])
		return self



	@utils._update_processing
	def integrate(self, dim='t', taper=True):
		'''
		integrates data using cum-trapezoids methods, tapers data by default
		see fobench.fiber.stools.signals.integrate_signal for more details
		'''
		axis = self.__axis__(dim)
		dx = self.dt if axis == 0 else self.spatial_interval

		if taper: self.taper(dim=dim)
		self.data = signals.integrate_signal(data=self.data, dx=dx, axis=axis)

		return self

	@utils._update_processing
	def differentiate(self, method='gradient', dim='t'):
		'''
		differentiate data in space or time, methods are 'gradient' or  'diff'
		see fobench.fiber.tools.signals.differentiate_signal for more details
		'''
		if method not in ['gradient', 'diff']:
			raise ValueError(f'\nInvalid method: "{method}". Choose on of:\n'
					' -"gradient"\n -"diff"')

		axis = self.__axis__(dim)
		self.data = signals.differentiate_signal(self.data, method=method,
										   axis=axis, dt=self.dt)

		return self

	def SNR(self, dim='t', results=False, plot_mode='pyqt'):
		'''
		computes signal to noise ratio, defined as ratio between mean and standard
		deviation of signal
		'''
		axis = self.__axis__(dim)
		snr = self.data.mean(axis=axis) / self.data.std(axis=axis)
		if plot_mode == 'mpl':
			warn('⚠️ matplotlib plotting not implemented for this method, '
				 'plotting using pyqtgraph instead')
			plot_mode = 'pyqt'

		if plot_mode == 'pyqt':
			plot_pyqt.plot_distance(distances=self.distances, channels_num=self.channels_num,
									   data=snr, y_label='SNR [-]', title='SNR Profile')
		if results:
			return snr

	def rmsa(self, window=None, dim='t', plot_mode='pyqt', results=False,
		  vmin=None, vmax=None):
		'''
		computes root mean square amplitude for record, dependign on dimension,
		window is either seconds ('t') or number of channels ('d')
		see fobench.tools.wavefield.rmsa for more details
		'''
		axis = self.__axis__(dim)
		if window is not None and dim == 't': window =  window*self.sampling_frequency
		if plot_mode == 'mpl':
			warn('⚠️ matplotlib plotting not implemented for this method, '
				 'plotting using pyqtgraph instead')
			plot_mode = 'pyqt'
		rmsa = wavefield.rmsa(data=self.data, axis=axis, window=window, dim=dim,
							times=self.times('unix'), distances = self.distances,
							channels_num=self.channels_num, vmin=vmin, vmax=vmax,
							plot_mode=plot_mode)
		if results:
			return rmsa

	def p2p_amp(self, dim='t', results=False, plot_mode='pyqt'):
		'''
		computes peak-to-peak amplitude of data in time or space
		see fobench.fiber.tools.wavefield.peak_to_peak_amp for more details
		'''
		axis = self.__axis__(dim)
		p2p_amplitude, up_index, down_index = wavefield.peak_to_peak_amp(self.data,
											 self.sampling_frequency, axis=axis)
		if plot_mode=='pyqt' and dim=='t':
				plot_pyqt.plot_distance(distances=self.distances, channels_num=self.channels_num,
							data=p2p_amplitude, y_label='P2P Amplitude', x_label='Channel',
							  title='Peak-to-Peak Amplitude Profile')

		if plot_mode=='pyqt' and dim=='d':
				plot_pyqt.plot_timeseries(timestamps=self.times(time_type='unix'), data=p2p_amplitude,
									y_label='Amplitude',
									title='Peak-to-Peak Amplitude over time')

		if results:
			return p2p_amplitude, up_index, down_index

	'''
	-----------------------------------------------------------------
	Plotting methods
	-----------------------------------------------------------------
	'''

	def fx_plot(self, norm=False, vmin=None, vmax=None, order=1, nfft=None, figsize=None,
				 show=True, cmap='viridis', results=False, file_name=None,
				 where=None, plot_mode='pyqt', **kwargs):
		'''
		computes frequency-distance plot
		see fobench.tools.wavefield.frequency_content for more details
		fobench.plotting.plotting_mpl for matplotlib plotting options
		'''
		axis = self.__axis__('t')

		fx, freqs =  wavefield.frequency_content(data=self.data, fs=self.sampling_frequency,
										   order=order, nfft=nfft, norm=norm, axis=axis)
		p95 = np.percentile(fx, 95)
		if vmin is None: vmin = 0
		if vmax is None: vmax = p95
		if plot_mode == 'pyqt':
			plot_pyqt.plot_2d_distance(distances=self.distances, channels_num=np.array(self.channels_num),
							  y_ticks=freqs, data=np.flip(np.rot90(fx, k=1), axis=0),
							  cmap=cmap, vmin=vmin, vmax=vmax, y_label='Frequency [Hz]',
							  title='Frequency content', cbar_label=self.units)

		elif plot_mode == 'mpl':
			plot.mpl_fx_plot(spec_matrix=fx[::-1], freqs=freqs, x=self.channels_num,
					 units_y='Energy', figsize=figsize, title=str(self.start_time.date),
					 cmap=cmap, file_name=file_name, vmin=vmin, vmax=vmax, **kwargs)

		if results:
			return fx, freqs

	def spectrum(self, channel, plot_mode='pyqt', norm=False, pre_processing=True,
				 order=1, pad=0, nfft=None, mode='spectrum', figsize=None,
				 nperseg=None, file_name=None, legend=True, results=False, **kwargs):
		"""
		compute spectrum of channel(s), mode can be 'spectrum' or 'psd'
		see fobench.tools.signals.signal_spectrum for more details
		for mpl plotting options see fobench.plotting.plotting_mpl.simple_spectrum
		"""

		axis = self.__axis__('t')
		if isinstance(channel, np.ndarray):
			channel = sorted(channel)
		if isinstance(channel, tuple):
			channel = list(range(min(channel), max(channel) + 1))
		elif isinstance(channel, list):
			channel = sorted(channel)
		else:
			channel = [channel]
		ch_idx = np.array([self.channels_num.index(ch) for ch in channel])
		o_signal = np.take(self.data, indices=ch_idx, axis=self.__axis__('d'))

		f, spec = signals.signal_spectrum(o_signal=o_signal, fs=self.sampling_frequency, mode=mode,
				norm=norm, order=order, nfft=nfft, pre_processing=pre_processing, pad=pad, nperseg=nperseg,
				axis=axis)
		spec = spec[:,0] if spec.ndim > 1 else spec  # reduced dimensionality of single spectrum

		if plot_mode=='pyqt':
			units = self.units if mode == 'spectrum' else f'{self.units.split(" ")[-1]}²/Hz'
			plot_pyqt.plot_spectral(frequencies=f, amplitudes=spec, y_label =f'{units.title()}',
								title= f'{mode.title()}' if mode=='spectrum' else f'{mode.upper()}', labels=channel)
		elif plot_mode=='mpl':
			units = self.units if mode == 'spectrum' else f'{self.units.split(" ")[-1]}$^{{2}}$/Hz'
			plot.simple_spectrum(spectra=np.array([spec]), freqs=f, channels=[channel], y_units=units, legend=legend, figsize=figsize,
						title=str(self.start_time.date), file_name=file_name, **kwargs)

		if results:
			return f, spec

	def channel_plot(self, channel, max_value=None, figsize=None, file_name=None,
					plot_mode='pyqt', **kwargs):
		'''
		generates simple plot of channel data
		for mpl mode see fobench.plotting.plotting_mpl.simple_plot for details
		'''
		if isinstance(channel, np.ndarray):
			channel = sorted(channel)
		if isinstance(channel, tuple):
			channel = list(range(min(channel), max(channel) + 1))
		elif isinstance(channel, list):
			channel = sorted(channel)
		else:
			channel = [channel]

		ch_idx = np.array([self.channels_num.index(ch) for ch in channel])
		selected = np.take(self.data, indices=ch_idx, axis=self.__axis__('d'))

		if plot_mode=='pyqt':
			t = self.times(time_type='unix')
			plot_pyqt.plot_timeseries(data=selected, timestamps=t, y_label=self.units,
							 dt=self.dt, title='Channel Plot', labels=channel)

		elif plot_mode=='mpl':
			t = self.times('matplotlib')
			plot.simple_plot(data=selected.T, t=t, channel=channel, units_y=self.units,
					max_value=max_value, figsize=figsize, title=str(self.start_time.date),
					file_name=file_name, **kwargs)

	def plot(self, vmin=None, vmax=None, figsize=None, show=True, cmap='seismic',
		  file_name=None, where=None, add_data=None, plot_mode='pyqt', **kwargs):
		'''
		generates plot of data, for more details see fobench.plotting.plotting_pyqt
		and .plotting_mpl
		'''
		if plot_mode == 'pyqt':
			t = self.times(time_type='unix')
			p95 = np.percentile(self.data, 95)
			vmin = -p95 if vmin is None else vmin
			vmax = p95 if vmax is None else vmax
			plot_pyqt.plot_2d_timeseries(timestamps=t, y_ticks=np.array(self.channels_num),
						data=self.data, y_label='Channel', dt=self.dt,
						title='', vmin=vmin, vmax=vmax, cbar_label=self.units,
						distances=self.distances)

		elif plot_mode == 'mpl':
			t = self.times(time_type='matplotlib')
			plot.gen_DAS_plot(data=self.data, t=t, channels=self.channels_num,
					 units_y=self.units, figsize=figsize, title=str(self.start_time.date),
					 cmap=cmap, file_name=file_name, vmin=vmin, vmax=vmax, add_data=add_data,
					 **kwargs)

	def channel_spectrogram(self, channel, norm=False, trace=False, figsize=None,
						cmap='viridis', file_name=None, 	freq_lim=None,  results=False,
						plot_mode='pyqt', vmin=None, vmax=None, **kwargs):
		'''
		computes and plots spectrogram for a 'channel', is normalized to maximum value if norm is True
		if using 'mpl' plot mode, see fobench.plotting.plotting_mpl.simple_spectrogram
		for details on plotting parameters
		'''

		axis = self.__axis__('t')
		channel = int(channel)
		index = self.channels_num.index(channel)
		data = self.data[:, index]

		f, t, Sxx = signals.signal_spectrogram(data=data, sampling_frequency=self.sampling_frequency,
										 axis=axis, norm=norm)

		if plot_mode == 'pyqt':
			t = self.times(time_type='unix')
			if vmin is None: vmin = 0
			if vmax is None: vmax = np.percentile(Sxx, 95)
			plot_pyqt.plot_2d_timeseries(timestamps=t, y_ticks=f, dt=self.dt,
						data=np.rot90(Sxx, k=-1), y_label='Frequency [Hz]',
						title=f'Spectrogram channel {channel}', cmap='viridis',
						vmin=vmin, vmax=vmax, cbar_label=self.units)

		elif plot_mode == 'mpl':
			t = self.times(time_type='matplotlib')
			plot.simple_spectrogram(data=Sxx, freq=f, t=t, units_y=self.units,
						trace=data if trace == True else None, figsize=figsize, cmap=cmap,
						title=str(self.start_time.date)+'  '+'Ch:'+str(channel),
						file_name=file_name, freq_lim=freq_lim, **kwargs)

		if results:
			return Sxx, f, t

	def record_section(self, channels, plot_mode='pyqt'):
		'''
		plots record section of multiple channels, if channels is tuple the range
		between lower and upper limit will be plotted, if list only channels in list
		will be plotted
		'''

		if isinstance(channels, np.ndarray):
			channels = channels.tolist()

		if isinstance(channels, tuple):
			ch0, chf = sorted(map(int, channels))
			ch0, chf = self.channels_num.index(ch0), self.channels_num.index(chf)
			ch_idx = slice(ch0, chf + 1)
		elif isinstance(channels, list):
			channels = list(map(int, channels))
			ch_idx = sorted(self.channels_num.index(ch) for ch in channels)
		else:
			raise TypeError(f'Invalid type for channels: {type(channels).__name__}. Expected tuple, list, or np.ndarray.')

		das_data = self.data[:, ch_idx]
		das_channels = np.array(self.channels_num)[ch_idx]

		if plot_mode=='pyqt':
			plot_pyqt.plot_record_section(timestamps=self.times('unix'), data=das_data,
								 title='Record Section', numbers=das_channels, dt=self.dt,
								 y_label='Channel')
		elif plot_mode=='mpl':
			plot.plot_record_section(signals=das_data, t=self.times('matplotlib'),
							channels=das_channels, date=str(self.start_time.date))

	def acf_profile(self, max_lag, plot_mode='pyqt', deconvolve=False,
					window_size=None, results=False, vmin=None, vmax=None, **imshow_kwargs):
		'''
		computes autocorrelation profile, see fiber.tools.wavefield.autocorrelation_profile
		for more details
		'''
		axis = self.__axis__('t')
		max_shift = int(max_lag*self.sampling_frequency)
		if max_shift >= self.num_points:
			raise ValueError('Selected max_shift is too large')
		acf = wavefield.autocorrelation_profile(self.data, max_shift, axis, plot_mode,
												deconvolve, self.total_channels,
												self.distances, self.channels_num,
												self.sampling_frequency,
												window_size=window_size, vmin=vmin,
												vmax=vmax, **imshow_kwargs)

		if results:
			return acf

	def spatial_coherence(self, max_lag, results=False, plot_mode='pyqt', vmin=None,
					   vmax=None):
		'''
		computes sptial coherence matrix, see fiber.tools.wavefield.spatial_coherence_matrix
		for more details
		'''
		coh = wavefield.spatial_coherence_matrix(data=self.data.T, max_lag=max_lag,
										   distances=self.distances,
										   fs=self.sampling_frequency,
										   channel_nums=self.channels_num,
										   plot_mode=plot_mode, results=results,
										   vmin=vmin, vmax=vmax)
		if results:
			return coh

	def view(self):
		'''
		launches the Fobench Data Viewer
		'''
		print(f'{"-"*65}\nStarting Fobench Data Viewer')
		app = QtWidgets.QApplication.instance()
		if app is None:
			app = QtWidgets.QApplication(sys.argv)
		self._viewer = Viewer(self)
		self._viewer.show()
		pg.exec()
		print(f'{"-"*65}')