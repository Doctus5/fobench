import functools
import inspect
import numpy as np
import h5py as h5
from tqdm import trange
from pyrocko.orthodrome import distance_accurate50m_numpy as deg2m

from obspy.core.trace import Trace as oTrace
from obspy.core.stream import Stream

from obspy.core import UTCDateTime as UTC

from pyrocko.util import str_to_time
from pyrocko.trace import Trace as pTrace

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
TOOLS USED EXCLUSIVE FOR Fiber CLASS
'''
STRAIN_UNIT_MAP = {	  # mapping for strain(rate) data
-1: 'integrated strain',
0: 'strain',
1: 'strain-rate',
2: 'strain-acceleration',
3: 'strain-jerk'}

VEL_UNIT_MAP = {		 # mapping for velocity data, e.g. Terra15 output
-1: 'm',
0: 'm/s',
1: 'm/s^2',
2: 'm/s^3'}

TEMP_UNIT_MAP = {		# mapping for temperature data
-1: 'integrated temperature',
0: 'temperature',
1: 'temperature rate',
2: 'temperature acceleration'}

UNKNOW_UNIT_MAP = { # mapping for data with unknow unit
-1: 'integrated units',
0: 'units',
1: 'd/dt units',
2: 'd/dt^2 units'}



def _update_processing(func):
	'''
	decorator function that updates the Fiber.processing attribute after each
	processing step, in case of integration or differentiation of the data
	also updates Fiber.units
	'''

	@functools.wraps(func)
	def wrapper(*args, **kwargs):

		func_name = func.__name__

		# extract all arguments
		bound_arguments = inspect.signature(func).bind(*args, **kwargs)
		bound_arguments.apply_defaults()
		args_dict = bound_arguments.arguments
		args_dict.pop('self')

		# function call
		result = func(*args, **kwargs)

		# append info to fiber instance
		fiber = args[0]
		fiber.processing.append({func_name : args_dict})

		if func_name in ['integrate', 'differentiate'] and args_dict['dim']=='t':
			# determine unit map
			if (fiber.sensing == 'das' or fiber.sensing == 'dss'):
				unit_map = STRAIN_UNIT_MAP
				if fiber.company == 'terra15':
					if fiber.attributes['properties']['data_product'] == 'velocity': unit_map = VEL_UNIT_MAP
				if fiber.company == 'sintela': unit_map = UNKNOW_UNIT_MAP
			elif fiber.sensing == 'dts':
				unit_map = TEMP_UNIT_MAP

			# find current unit
			try: key = [i for i in unit_map if unit_map[i] == fiber.units][0]
			except: key = int(fiber.units[2])

			# depending on operation change key
			if func_name == 'integrate': key -= 1
			elif func_name == 'differentiate': key += 1

			# assign new unit
			try: fiber.units = unit_map[key]
			except: fiber.units = f'd^{key}/dt {unit_map[0]}'

		return result
	return wrapper

def instr_corr(data: np.ndarray = None, attributes: dict = None, target: str = 'strain-rate',
				terra15_gl: float = None) -> tuple[np.ndarray, str, np.ndarray, list, int, float]:
	'''
	performs instrument correction and data conversion for various instrument types

	Parameters
	----------
	data : np.ndarray
		data of Fiber class
	attributes : dict
		Fiber class attributes
	target : str, optional
		target unit. The default is 'strain-rate'.
	terra15_gl : float, optional
		gauge length for velocity to strain-rate conversion for Terra15 data.
		If not specified, original gauge length from Fiber.gauge_length is taken

	Returns
	-------
	data : np.ndarray
		corrected data
	target : str
		(new) unit of measurement
	attributes['channels'] : np.ndarray
		(new) channels
	attributes['channels_num'] : list
		(new) channel_numbers
	attributes['gauge_length'] : float
		(new) gauge length

		'''

	format, company, units = attributes['format'], attributes['company'], attributes['units']

	if format == 'tdms' and company == 'silixa': # Silixa TDMS

		if units == 'counts' and target == 'strain-rate':

			i_cst = 116E-9 # meters per radians.
			gauge_L = attributes['gauge_length'] # gauge lenght in meters.
			digital_N = 2**13 # magic number linked to the digitalization of the data. why not 2**16?
			fs = attributes['o_sampling_frequency'] # sampling frequency, which can be 1000 Hz for raw data.
			factor = i_cst*(fs/gauge_L)/digital_N # strain Rate per counts.
			data = np.multiply(data,factor)


	# if (format == 'h5' or format == 'hdf5') and company == 'febus': # FEBUS HDF5

	elif (format == 'h5' or format == 'hdf5') and company == 'silixa': # Silixa HDF5

		if units == 'counts' and target == 'strain-rate':

			i_cst = 116E-9 # meters per radians.
			gauge_L = attributes['gauge_length'] # gauge lenght in meters.
			digital_N = 2**13 # magic number linked to the digitalization of the data. why not 2**16?
			fs = attributes['o_sampling_frequency'] # sampling frequency, which can be 1000 Hz for raw data.
			factor = i_cst*(fs/gauge_L)/digital_N # strain Rate per counts.
			data = np.multiply(data,factor)

	# if format == 'npy' and company == 'bam': # .npy format for BAM. This might fail always since the unit is NON-COMMERCIAL!

	elif (format == 'h5' or format == 'hdf5') and company == 'terra15': # Terra15 HDF5
		if units == 'velocity' and target == 'strain-rate':
			gl = attributes['gauge_length'] if terra15_gl is None else terra15_gl
			gauge_samples = int(round(gl / attributes['spatial_interval']))
			gl = gauge_samples * attributes['spatial_interval']
			data = (data[:, gauge_samples:] - data[:, :-gauge_samples]) / gl
			n_decim = int(gauge_samples/2)
			attributes['channels'] = attributes['channels'][n_decim:-n_decim]
			attributes['channels_num'] = attributes['channels_num'][n_decim:-n_decim]
			attributes['total_channels'] = attributes['channels'].size
			attributes['gauge_length'] = gl

	elif (format == 'h5' or format == 'hdf5') and company == 'asn': # ASN OptoDAS HDF5 (It can be a bit more complex, so I'm trying to make it simple!)

		if units == 'rad/(strain*m)' and target == 'strain-rate':

			data = data / attributes['conv_factor'] # divide by sensitivities. It seems they already provide the conversion factor

	elif (format == 'h5' or format == 'hdf5') and company == 'quantx': # QuantX OptoaSense HDF5 (CHECK THIS!! WITH VERIFICATION OR CALIBRATION).

		if target == 'strain-rate':

			i_cst = 116E-9 # meters per radians.
			gauge_L = attributes['gauge_length'] # gauge lenght in meters.
			digital_N = int(attributes['units'][-4]) ** int(attributes['units'][-2:]) # magic number linked to the digitalization of the data. why not 2**ints
			fs = attributes['o_sampling_frequency'] # sampling frequency, which can be 1000 Hz for raw data.
			factor = i_cst*(fs/gauge_L)/digital_N # strain Rate per counts.
			data = np.multiply(data,factor)

	# ####################################################
	# CAUTION!! NON OFFICIAL / EXPERIMENTAL FORMATS, ONLY FOR SPECIAL CASES.
	# ####################################################

	elif format == 'npz' and company == 'bam': # .npy format for BAM. This might fail always since the unit is NON-COMMERCIAL!

		if units == 'counts' and target == 'strain':

			factor = 1E-6 / 18.4 # strain per count (WEIRD!)
			data = np.multiply(data,factor)

	elif (format == 'h5' or format == 'hdf5') and company == 'michelle': # Michelle HDF5 decimated from Silixa

		if units == 'counts' and target == 'strain-rate':

			i_cst = 116E-9 # meters per radians.
			gauge_L = attributes['gauge_length'] # gauge lenght in meters.
			digital_N = 2**13 # magic number linked to the digitalization of the data. why not 2**16?
			fs = attributes['o_sampling_frequency'] # sampling frequency, which can be 1000 Hz for raw data.
			factor = i_cst*(fs/gauge_L)/digital_N # strain Rate per counts.
			data = np.multiply(data,factor)

	return data, target, attributes['channels'], attributes['channels_num'], attributes['total_channels'], attributes['gauge_length']


def interpolate_channels(n_ch, x_ch, y_ch, z_ch, system='decimal', err=None, spacing=None):
	'''
	Co-authors: --
	Description:
		Do a linear interpolation between sections of georeferenced channels to georeference the non-located channels.
		Inputs must be the georeferences channels in ascending order.
	:Params:
		- n_ch(type:Numpy): 1D array of channel number.
		- x_ch(type:Numpy): 1D array of X (longitude) coordinates of the channels specified in "n_ch".
		- y_ch(type:Numpy): 1D array of Y (latitude) coordinates of the channels specified in "n_ch".
		- z_ch(type:Numpy): 1D array of Z (depth - meters) coordinates of the channels specified in "n_ch".
		- system(type:String): Defined the receiving coordinate systems for X and Y. It can be 'decimal' for decimal degrees
		or 'utm' for Universal Transverse Mercator. Default is 'decimal'.
		- err(type:Float - Optional): maximum accepted error of gauge length from interpolation in decimals 0.0 = 0% and 1.0 = 100%.
		- spacing(type:Int or Float - Optional): real channel spacing value in meters.
	:Return:
		- NA.
	'''

	new_ch, new_x, new_y, new_z =  [], [], [], []

	for i in range(n_ch.size-1):

		ch_start, ch_end = n_ch[i], n_ch[i+1]

		# deltas or differences between extremes.
		dch =  n_ch[i+1] - n_ch[i]
		dx = x_ch[i+1] - x_ch[i]
		dy = y_ch[i+1] - y_ch[i]
		dz = z_ch[i+1] - z_ch[i]

		# channels = np.linspace(ch_start, ch_end, ch_end - ch_start + 1)
		N = np.linspace(0, 1, ch_end - ch_start + 1) # number of channels in between, counting extremes.
		CH = n_ch[i] + dch * N
		X = x_ch[i] + dx * N
		Y = y_ch[i] + dy * N
		Z = z_ch[i] + dz * N

		ii = 0 if i == 0 else 1 # index for start aappending new georefereced channels.

		new_ch += list(CH[ii:].astype(int))
		new_x += list(X[ii:])
		new_y += list(Y[ii:])
		new_z += list(Z[ii:])

		if err != None:

			if system == 'decimal':

				r1 = deg2m(y_ch[i+1], x_ch[i+1], y_ch[i], x_ch[i], implementation='python')[0] # lat/lon distance in meters.
				r2 = np.sqrt(r1**2 + dz**2) # total calculated distance between known locations in 3D.

			if system == 'utm':

				r2 = np.sqrt(dx**2 + dy**2 + dz**2) # total calculated distance between known locations in 3D with UTM.

			calc_spacing = (1/(ch_end - ch_start)) * r2 # calculated channel spacing from georeferencing.
			dev = (calc_spacing - spacing) / spacing # calculated error between new channel spacing and real.
			print(f'Fiber section {i+1} -> channel spacing of {calc_spacing} m. Original spacing is {spacing} m.')

			if np.abs(dev) > err:

				# A warning is raised but program continues.
				print(f'WARNING: Error of calculated channels spacing in fiber section {i+1} is of {dev*100}%.\nCheck control points.')
				return None # Calculation is interrupted for not fulfilling standards.

	return np.array(new_ch).astype(int), np.array(new_x), np.array(new_y), np.array(new_z)


'''
####################################################
Signal Analysis functions below...
####################################################
'''

def auto_cascf(in_fs, out_fs):
	'''
	Gives the optimal cascadian decimation factors based on initial sampling frequency and final sampling frequency.
	:Params(type):
		- data(type:numpy): matrix data (2D) of the DAS Class.
		- in_fs(type:int or float): original (input) sampling frequency.
		- out_fs(type:int or float): desired (output) sampling frequency.
	:Return(type):
		- factors(type:list): list of factors to use to get from input to output frequencies.
	'''

	factors = []
	pref_factors = [2, 3, 4, 5] # prefered decimation factors in order of preference.
	n = in_fs / out_fs # overall decimation factor

	for d in pref_factors:
		while n % d == 0:
			factors.append(d)
			n //= d
	if n > 1:
		factors.append(n)  # Include any remaining prime factor

	return factors

def spatial_upsampling(das_class):
	'''
	Co-authors: --
	Description:
		Tool for upscaling spatialy the DDSS (DAS) data by double. Creates an interpolation between consecutive
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
	Description:
		Tool for downscaling spatialy the DDSS (DAS) data by half. Erase one channel between consecutive
		channels to simulate a decrease spatial resolution.
	:Params:
		- das_class(type:DAS): an initialized DAS Class with data.
	:Return:
		- new_data(type:numpy): 2D matrix containing the new spatial downsampled data.
		- new_channels_num(type: numpy): a list containing the new numbers of the channels, where the inermediate ones are eliminated.
	'''

	new_data = das_class.data[:,::2]
	new_channels_num = das_class.channels_num[::2] #only if the label of the channel wants to be fixed (0,2,4,6,...,N)
	#new_channels_num = [i for i in range(0,int(len(das_class.channels_num)/2))] #channel numbers change due to the downsampling (0,1,2,3,...,N/2)
	#new_channels_num = das_class.channels_num[:int(len(das_class.channels_num)/2)] #channel numbers change due to the downsampling (0,1,2,3,...,N/2)

	return new_data, new_channels_num

def to_traces(Fiber, t_type: str)-> Stream:
	'''
	Creates an obpsy/pyrocko Stream object and fill it with Traces in it. Each Trace
	would represent each channel of the DAS Class, including the metadata which
	are attributes of the Trace Class. This is mainly done so users can have access
	to obspy tools with this data. However, it can be slower and memory demanding.

	Parameters
	----------
	Fiber : fobench.Fiber object
		Fiber class object
	t_type : str
		type of stream to return, 'obspy' or 'pyrocko'

	Returns
	-------
	Stream
		Stream object

	'''
	stream = Stream() if t_type == 'obspy' else []

	for i in trange(Fiber.total_channels, desc='Creating Stream'):

		if t_type == 'obspy':
			trace = oTrace(data=Fiber.data[:,i])
			trace.stats.network = Fiber.fiber
			trace.stats.station = str(Fiber.channels_num[i]).zfill(5)
# 			trace.stats.npts = self.num_points #+ 1
			trace.stats.sampling_rate = Fiber.sampling_frequency
			trace.stats.delta = Fiber.dt
			trace.stats.starttime = Fiber.start_time
			trace.stats.calib = instr_corr(np.array(1), attributes=vars(Fiber))
			trace.stats.channel = 'H'
			#trace.stats.endtime = self.end_time
			stream.append(trace)
# 			print(stream)

		if t_type == 'pyrocko':
			trace = pTrace(ydata=Fiber.data[:,i])
			trace.network = Fiber.fiber
			trace.station = str(Fiber.channels_num[i]).zfill(5)
			trace.deltat = Fiber.dt
			trace.tmin = str_to_time(Fiber.start_time.isoformat().replace('T',' '))
			trace.tmax = str_to_time(Fiber.end_time.isoformat().replace('T',' '))
			stream.append(trace)

	return stream

def return_times(Fiber, time_type: str)-> np.ndarray:
	'''
	returns 1D array containing time-steps of data in the specified format

	Parameters
	----------
	Fiber : fobench.Fiber object
		Fiber class object
	time_type : str
		time format to return. options are 'UTCDateTime', 'isoformat', 'matplotlib'
		or 'unix´'

	Raises
	------
	ValueError
		unrecognized time format.

	Returns
	-------
	t : TYPE
		1D array containing time-steps of data in the specified format.

	'''
	converters = {
		'UTCDateTime': lambda t: t,
		'isoformat': lambda t: t.isoformat(),
		'matplotlib': lambda t: t.matplotlib_date,
		'unix': lambda t: t.timestamp,
	}

	if time_type not in converters:
		raise ValueError(
			f'Unrecognized time format "{time_type}"! Please choose one of:\n'
			  ' -"UTCDateTime"\n -"isoformat"\n -"matplotlib"\n -"unix"'
		)

	times = [Fiber.start_time + i * Fiber.dt for i in range(Fiber.data.shape[0])]
	return np.array([converters[time_type](t) for t in times])

def trim_time(t0: UTC | str, tf: UTC | str, data: np.ndarray, times: np.ndarray,
			  start_time: UTC, end_time: UTC) -> tuple[np.ndarray, UTC, UTC]:
	'''
	trims data in Fiber class in time dimension between given start and end times

	Parameters
	----------
	t0, tf : UTC datetime object or str
		new desired start and end timea of data
	data : np.ndarray
		original data to trim
	times : np.ndarray
		timestamps of Fiber class
	start_time, end_time : UTC datetime
		original start and end time of data

	Raises
	------
	ValueError
		t0 > tf

	Returns
	-------
	data : np.ndarray
		trimmed data
	start_time, end_time : UTC datetime
		new start and end time of data
	'''
	t0, tf = UTC(t0), UTC(tf)

	t0 = max(t0, start_time)
	tf = min(tf, end_time)

	if tf < t0: raise ValueError("End time (tf) must be after start time (t0).")
	t0_pos = max(0, np.searchsorted(times, t0, side='right') - 1)
	tf_pos = max(0, np.searchsorted(times, tf, side='right') - 1)

	data = data[t0_pos:tf_pos, :]
	start_time = times[t0_pos]
	end_time = times[tf_pos]

	return data, start_time, end_time