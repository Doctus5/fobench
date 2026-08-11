import functools
import inspect
import numpy as np
import h5py as h5
from tqdm import trange
from pyrocko.orthodrome import distance_accurate50m_numpy as deg2m
from pyrocko.util import str_to_time
from pyrocko.trace import Trace as pTrace
from obspy.core.trace import Trace as oTrace
from obspy.core.stream import Stream
from obspy.core import UTCDateTime as UTC

def scan_hdf5(path, recursive=True, tab_step=2):
	"""Function to sscan hdf5 file hierarchycally."""
	def scan_node(g, tabs=0):
		print(" " * tabs, g.name)
		for k, v in g.items():
			if isinstance(v, h5.Dataset):
				print(" " * tabs + " " * tab_step + " -", v.name)
			elif isinstance(v, h5.Group) and recursive:
				scan_node(v, tabs=tabs + tab_step)
	with h5.File(path, "r") as f:
		scan_node(f)


"""TOOLS USED EXCLUSIVE FOR Fiber CLASS"""

STRAIN_UNIT_MAP = {	  # mapping for strain(rate) data
                   -1: "integrated strain",
                   0: "strain",
                   1: "strain-rate",
                   2: "strain-acceleration",
                   3: "strain-jerk"}

VEL_UNIT_MAP = {		 # mapping for velocity data, e.g. Terra15 output
                -1: "m",
                0: "m/s",
                1: "m/s^2",
                2: "m/s^3"}

TEMP_UNIT_MAP = {		# mapping for temperature data
                 -1: "integrated temperature",
                 0: "temperature",
                 1: "temperature rate",
                 2: "temperature acceleration"}

UNKNOW_UNIT_MAP = {      # mapping for data with unknow unit
                   -1: "integrated units",
                   0: "units",
                   1: "d/dt units",
                   2: "d/dt^2 units"}


def _update_processing(func):
	"""Decorator function that updates the Fiber.processing attribute after each
	processing step, in case of integration or differentiation of the data
	also updates Fiber.units
	"""

	@functools.wraps(func)
	def wrapper(*args, **kwargs):

		func_name = func.__name__

		# extract all arguments
		bound_arguments = inspect.signature(func).bind(*args, **kwargs)
		bound_arguments.apply_defaults()
		args_dict = bound_arguments.arguments
		args_dict.pop("self")

		result = func(*args, **kwargs) # function call

		# append info to fiber instance
		fiber = args[0]
		fiber.processing.append({func_name : args_dict})

		if func_name in ["integrate", "differentiate"] and args_dict["dim"]=="t":
			# determine unit map
			if (fiber.sensing == "das" or fiber.sensing == "dss"):
				unit_map = STRAIN_UNIT_MAP
				if fiber.company == "terra15" and "strain" not in fiber.units:
					if fiber.attributes["properties"]["data_product"] == "velocity": unit_map = VEL_UNIT_MAP
				if fiber.company == "sintela": unit_map = UNKNOW_UNIT_MAP
			elif fiber.sensing == "dts":
				unit_map = TEMP_UNIT_MAP

			# find current unit
			try: key = [i for i in unit_map if unit_map[i] == fiber.units][0]
			except: key = int(fiber.units[2])

			# depending on operation change key
			if func_name == "integrate": key -= 1
			elif func_name == "differentiate": key += 1

			# assign new unit
			try: fiber.units = unit_map[key]
			except: fiber.units = f"d^{key}/dt {unit_map[0]}"

		return result
	return wrapper

def instr_corr(data: np.ndarray = None, attributes: dict = None, target: str = "strain-rate",
				terra15_gl: float = None) -> tuple[np.ndarray, str, np.ndarray, list, int, float]:
	"""Performs instrument correction and data conversion for various instrument types

	Parameters
	----------
	data : np.ndarray
		Data of Fiber class.
	attributes : dict
		Fiber class attributes.
	target : str, optional
		Target unit. The default is "strain-rate".
	terra15_gl : float, optional
		Gauge length for velocity to strain-rate conversion for Terra15 data.
		If not specified, original gauge length from Fiber.gauge_length is taken.

	Returns
	-------
	data : np.ndarray
		Corrected data.
	target : str
		(New) unit of measurement.
	attributes["channels"] : np.ndarray
		(New) channels.
	attributes["channels_num"] : list
		(New) channel_numbers.
	attributes["gauge_length"] : float
		(New) gauge length.

	"""

	format, company, units = attributes['format'], attributes['company'], attributes['units']

	if format in ("tdms", "h5", "hdf5") and company == "silixa": # Silixa TDMS/HDF5
		if units == "counts" and target == "strain-rate":
			i_cst = 116E-9 # meters per radians.
			gauge_L = attributes["gauge_length"] # gauge length in meters.
			digital_N = 2**13 # magic number linked to the digitalization of the data. why not 2**16?
			fs = attributes["o_sampling_frequency"] # sampling frequency, which can be 1000 Hz for raw data.
			factor = i_cst*(fs/gauge_L)/digital_N # strain-rate per counts.
			data = np.multiply(data, factor)

	elif (format == "h5" or format == "hdf5") and company == "terra15": # Terra15 HDF5
		if units == "m/s" and target == "strain-rate":
			gl = attributes["gauge_length"] if terra15_gl is None else terra15_gl
			gauge_samples = int(round(gl / attributes["spatial_interval"]))
			gauge_samples = 1 if gauge_samples < 1 else gauge_samples
			gl = gauge_samples * attributes["spatial_interval"]
			print(f"\n⚠️ Applying the nearest possible gauge length: {gl}m")
			data = (data[:, gauge_samples:] - data[:, :-gauge_samples]) / gl
			n_left = gauge_samples // 2
			n_right = gauge_samples - n_left
			if gauge_samples != 0:
				attributes["channels"] = attributes["channels"][n_left:-n_right]
				attributes["channels_num"] = attributes["channels_num"][n_left:-n_right]
				attributes["distances"] = attributes["distances"][n_left:-n_right]
				attributes["total_channels"] = attributes["channels"].size
			attributes["gauge_length"] = gl
		else:
			print("\n⚠️ Data not in velocity units, doing nothing ...")

	elif (format == "h5"	 or format == "hdf5"	) and company == "asn": # ASN OptoDAS HDF5 (It can be a bit more complex, so I"	m trying to make it simple!)
		if units == "rad/(strain*m)" and target == "strain-rate":
			data = data / attributes["conv_factor"] # divide by sensitivities. It seems they already provide the conversion factor

	elif (format == "h5" or format == "hdf5") and company == "quantx": # QuantX OptoaSense HDF5 (CHECK THIS!! WITH VERIFICATION OR CALIBRATION).
		if target == "strain-rate":
			i_cst = 116E-9 # meters per radians.
			gauge_L = attributes["gauge_length"] # gauge lenght in meters.
			digital_N = int(attributes["units"][-4]) ** int(attributes["units"][-2:]) # magic number linked to the digitalization of the data. why not 2**ints
			fs = attributes["o_sampling_frequency"] # sampling frequency, which can be 1000 Hz for raw data.
			factor = i_cst*(fs/gauge_L)/digital_N # strain Rate per counts.
			data = np.multiply(data,factor)

	# ####################################################
	# CAUTION!! NON OFFICIAL / EXPERIMENTAL FORMATS, ONLY FOR SPECIAL CASES.
	# ####################################################

	elif format == "npz" and company == "bam": # .npy format for BAM. This might fail always since the unit is NON-COMMERCIAL!
		if units == "counts" and target == "strain":
			factor = 1E-6 / 18.4 # strain per count (WEIRD!)
			data = np.multiply(data,factor)

	elif (format == "h5" or format == "hdf5") and company == "michelle": # Michelle HDF5 decimated from Silixa
		if units == "counts" and target == "strain-rate":
			i_cst = 116E-9 # meters per radians.
			gauge_L = attributes["gauge_length"] # gauge lenght in meters.
			digital_N = 2**13 # magic number linked to the digitalization of the data. why not 2**16?
			fs = attributes["o_sampling_frequency"] # sampling frequency, which can be 1000 Hz for raw data.
			factor = i_cst*(fs/gauge_L)/digital_N # strain Rate per counts.
			data = np.multiply(data,factor)

	return data, target, attributes['channels'], attributes['channels_num'], attributes['total_channels'], attributes['gauge_length'], attributes['distances']


def interpolate_channels(n_ch: np.ndarray, x_ch: np.ndarray, y_ch: np.ndarray,
                         z_ch: np.ndarray, system: str = "decimal", err: float = None,
                         spacing: int | float = None)-> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
	"""Do a linear interpolation between sections of georeferenced channels to
	georeference the non-located channels. Inputs are the georeferenced channels
	in ascending order.

	Parameters
	----------
	n_ch : np.ndarray
		1D array of channel number.
	x_ch : np.ndarray
		1D array of X (longitude) coordinates of the channels specified in ``"n_ch"``.
	y_ch : np.ndarray
		1D array of Y (latitude) coordinates of the channels specified in ``"n_ch"``.
	z_ch : np.ndarray
		1D array of Z (depth - meters) coordinates of the channels specified in ``"n_ch"``.
	system : str
		Defines the output coordinate system for X and Y. It can be ``'decimal'``
		for decimal degrees or ``'utm'`` for Universal Transverse Mercator.
	err: float
		The maximum accepted error of gauge length from interpolation in decimals 0.0 = 0% and 1.0 = 100%.
	spacing : int, float
		Real channel spacing value in meters.

	Returns
	-------
    new_ch : np.ndarray
        -
    new_x : np.ndarray
        X coordinates.
    new_y : np.ndarray
        Y coordinates.
    new_z : np.ndarray
        Z coordinates.

	"""

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
			if system == "decimal":
				r1 = deg2m(y_ch[i+1], x_ch[i+1], y_ch[i], x_ch[i], implementation="python")[0] # lat/lon distance in meters.
				r2 = np.sqrt(r1**2 + dz**2) # total calculated distance between known locations in 3D.
			elif system == "utm":
				r2 = np.sqrt(dx**2 + dy**2 + dz**2) # total calculated distance between known locations in 3D with UTM.

			calc_spacing = (1/(ch_end - ch_start)) * r2 # calculated channel spacing from georeferencing.
			dev = (calc_spacing - spacing) / spacing # calculated error between new channel spacing and real.
			print(f"Fiber section {i+1} -> channel spacing of {calc_spacing} m. Original spacing is {spacing} m.")

			if np.abs(dev) > err:
				print(f"\n⚠️ Error of calculated channels spacing in fiber section "
					f"{i+1} is of {dev*100}%.\nCheck control points.")
				return None # Calculation is interrupted for not fulfilling standards.

	return np.array(new_ch).astype(int), np.array(new_x), np.array(new_y), np.array(new_z)


"""Signal Analysis functions below"""

def auto_cascf(in_fs: int | float, out_fs: int | float) -> list:
	"""Gives the optimal cascadian decimation factors based on initial sampling
	frequency and desired final sampling frequency.

	Parameters
	----------
	in_fs : int, float
		Original (input) sampling frequency.
	out_fs: int, float:
		Desired (output) sampling frequency.

	Returns
	-------
	factors: list
		List of factors to use to get from input to output frequency.
	"""

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

def spatial_upsampling(Fiber)-> tuple[np.ndarray, list]:
	"""Function for spatially upsampling the DAS, doubling the number of channels.
	Creates an interpolation between consecutive channels to simulate an increase
	in spatial resolution.

	Parameters
	----------
	Fiber : fobench.Fiber object
		Fiber class object.

	Returns
	-------
	new_data : np.ndarray
		2D matrix containing the new spatial upsampled data.
	new_channels_num : list
		List containing the new numbers of the channels, with newly created intermediary ones.
	"""

	new_channels_num = [Fiber.channels_num[0]]
	shape = (len(Fiber.data[:,0]),1)
	new_data = Fiber.data[:,0].reshape(shape) #reshaping to not affect original dimensionality.

	for i in range(Fiber.total_channels-1):
		first, second = new_data[:,-1].reshape(shape), Fiber.data[:,i+1].reshape(shape)
		inter = (first + second) / 2
		new_data = np.concatenate((new_data, inter), axis=1)
		new_data = np.concatenate((new_data, second), axis=1)
		inter_num = (new_channels_num[-1] + Fiber.channels_num[i+1]) / 2
		new_channels_num += ([inter_num, Fiber.channels_num[i+1]])

	return new_data, new_channels_num


def spatial_downsampling(Fiber)-> tuple[np.ndarray, list]:
	"""Function for spatially downsampling DAS data by half. Erases one channel between consecutive
		channels to simulate a decrease of the spatial resolution.

	Parameters
	----------
	Fiber : fobench.Fiber object
		Fiber class object.

	Returns
	-------
	new_data : np.ndarray
		2D matrix containing the new spatial downsampled data.
	new_channels_num : list
		List containing the new numbers of the channels, with the inermediate ones eliminated.
	"""

	new_data = Fiber.data[:,::2]
	new_channels_num = Fiber.channels_num[::2] #only if the label of the channel wants to be fixed (0,2,4,6,...,N)
	#new_channels_num = [i for i in range(0,int(len(Fiber.channels_num)/2))] #channel numbers change due to the downsampling (0,1,2,3,...,N/2)
	#new_channels_num = Fiber.channels_num[:int(len(Fiber.channels_num)/2)] #channel numbers change due to the downsampling (0,1,2,3,...,N/2)

	return new_data, new_channels_num

def to_traces(Fiber, t_type: str)-> Stream:
	"""Creates an obpsy or pyrocko Stream object and fills it with Traces. Each Trace
	represents a channel of the DAS Class, including the metadata which
	from the attributes of the Trace Class. This is mainly done so users can have access
	to obspy tools with their data. However, it can be slow and memory demanding.

	Parameters
	----------
	Fiber : fobench.Fiber object
		Fiber class object.
	t_type : str
		type of stream to return, ``'obspy'`` or ``'pyrocko'``.

	Returns
	-------
	stream
		Stream object.

	"""

	stream = Stream() if t_type == 'obspy' else []

	for i in trange(Fiber.total_channels, desc="Creating Stream"):
		if t_type == "obspy":
			trace = oTrace(data=Fiber.data[:,i])
			trace.stats.network = Fiber.fiber
			trace.stats.station = str(Fiber.channels_num[i]).zfill(5)
# 			trace.stats.npts = self.num_points #+ 1
			trace.stats.sampling_rate = Fiber.sampling_frequency
			trace.stats.delta = Fiber.dt
			trace.stats.starttime = Fiber.start_time
			trace.stats.calib = instr_corr(np.array(1), attributes=vars(Fiber))
			trace.stats.channel = "H"
			#trace.stats.endtime = self.end_time
			stream.append(trace)
# 			print(stream)

		elif t_type == "pyrocko":
			trace = pTrace(ydata=Fiber.data[:,i])
			trace.network = Fiber.fiber
			trace.station = str(Fiber.channels_num[i]).zfill(5)
			trace.deltat = Fiber.dt
			trace.tmin = str_to_time(Fiber.start_time.isoformat().replace("T"," "))
			trace.tmax = str_to_time(Fiber.end_time.isoformat().replace("T"," "))
			stream.append(trace)

	return stream


def to_xarray(Fiber, name=None, use_distance=False):
    
    try:
        import xarray as xr
    except ImportError as exc:
        raise ImportError("xarray package is required, but no module is found.")
    
    # create the different labels and coordinates
    data_label = name or Fiber.units or "data"
    times = Fiber.times("datetime64")
    channels = np.asarray(Fiber.channels_num)
    distances = np.asarray(Fiber.distances, dtype=float)
    attrs = clean_metadata(Fiber)
    
    if use_distance:
        dims = ("time", "distance")
        coords = {"time": ("time",times),
                  "distance": ("distance",distances),
                  "channel": ("distance",channels)}
    else:
        dims = ("time", "channel")
        coords = {"time": ("time",times),
                "channel": ("channel",channels),
				"distance": ("channel",distances)}
    
    return xr.DataArray(data=Fiber.data, dims=dims, coords=coords, name=data_label, attrs=attrs)
        
    
def clean_metadata(Fiber) -> dict:
	"""Returns a clean metadata of Fiber

	Args:
		Fiber : _description_

	Returns:
		dict: clean metadata dicitonary of Fiber class
	"""
    
	attrs = {"fiber": Fiber.fiber,
			"company": Fiber.company,
			"format": Fiber.format,
			"units": Fiber.units,
			"sampling_frequency": Fiber.sampling_frequency,
			"o_sampling_frequency": Fiber.o_sampling_frequency,
			"dt": Fiber.dt,
			"spatial_interval": Fiber.spatial_interval,
			"gauge_length": Fiber.gauge_length,
			"channel_offset": Fiber.channel_offset,
			"start_time": Fiber.start_time.isoformat(),
			"end_time": Fiber.end_time.isoformat(),
			"time_length": Fiber.time_length,
			"total_channels": Fiber.total_channels,
			"conv_factor": Fiber.conv_factor,
			"sensing": Fiber.sensing,
			# "basefile": getattr(Fiber, "__basefile__", None),
			"source_files": getattr(Fiber, "__filepath__", None),
			# "processing": Fiber.processing,
			# "properties": Fiber.properties
   			}
    
	return attrs


def return_times(Fiber, time_type: str)-> np.ndarray:
	"""Returns a 1D array containing time-steps of data in the specified format.

	Parameters
	----------
	Fiber : fobench.Fiber object
		``Fiber`` class object.
	time_type : str
		time format to return. options are ``'UTCDateTime'``, ``'isoformat'``,
		``datetime64``, ``'matplotlib'`` or ``'unix'``

	Raises
	------
	ValueError
		Unrecognized time format.

	Returns
	-------
	t : np.array
		1D array containing time-steps of data in the specified format.

	"""
	converters = {"UTCDateTime": lambda t: t,
				  "isoformat": lambda t: t.isoformat(),
				  "datetime64": lambda t: np.datetime64(t.datetime, "ns"),
				  "matplotlib": lambda t: t.matplotlib_date,
				  "unix": lambda t: t.timestamp}

	if time_type not in converters:
		raise ValueError(
			f"\n⚠️ Unrecognized time format '{time_type}'! Please choose one of:\n"
			" -'UTCDateTime'\n -'isoformat'\n -'matplotlib'\n -'unix'")

	times = [Fiber.start_time + i * Fiber.dt for i in range(Fiber.data.shape[0])]

	return np.array([converters[time_type](t) for t in times])


def _to_seconds(value) -> float:
    
    if isinstance(value, (int, float, np.integer, np.floating)):
        
        return float(value)
    
    return float(UTC(value))


def _interp_nan_along_axis0_inplace(data: np.ndarray) -> np.ndarray:
    
    x = np.arange(data.shape[0])
    
    for j in range(data.shape[1]):
        
        y = data[:, j]
        valid = np.isfinite(y)
        
        if np.any(valid) and not np.all(valid):
            
            y[~valid] = np.interp(x[~valid], x[valid], y[valid])
    
    return data


def time_concatenation(data1: np.ndarray, data2: np.ndarray, start1, end1, dt1: float, start2, end2, 
                      dt2: float, axis: int = 0, overlap: str = 'data2', gap: str = 'nan', tolerance: float = 1.0,
                      out: np.ndarray = None) -> np.ndarray:
	"""2D concatenaiton of matrices with time series data.

	Parameters
	----------
	data1 : np.ndarray
		_description_
	data2 : np.ndarray
		_description_
	start1 : _type_
		_description_
	end1 : _type_
		_description_
	dt1 : float
		_description_
	start2 : _type_
		_description_
	end2 : _type_
		_description_
	dt2 : float
		_description_
	axis : int, optional
		_description_, by default 0
	overlap : str, optional
		_description_, by default 'data2'
	gap : str, optional
		_description_, by default 'nan'
	tolerance : float, optional
		_description_, by default 1.0
	out : np.ndarray, optional
		_description_, by default None

	Returns
	-------
	np.ndarray
		_description_

	Raises
	------
	ValueError
		_description_
	ValueError
		_description_
	ValueError
		_description_
	ValueError
		_description_
	ValueError
		_description_
	ValueError
		_description_
	ValueError
		_description_
	ValueError
		_description_
	ValueError
		_description_
	ValueError
		_description_
	ValueError
		_description_
	"""

	if data1.ndim != 2 or data2.ndim != 2:
		raise ValueError('Both inputs must be 2D matrices.')
	if axis not in (0, 1):
		raise ValueError('axis must be 0 or 1.')
	if overlap not in ('data1', 'data2', 'mean'):
		raise ValueError('overlap must be "data1", "data2", or "mean".')
	if gap not in ('nan', 'zero', 'linear'):
		raise ValueError('gap must be "nan", "zero", or "linear".')

	a = np.moveaxis(data1, axis, 0)
	b = np.moveaxis(data2, axis, 0)
	if a.shape[1] != b.shape[1]:
		raise ValueError('Input shapes are incompatible on non-concatenation axis.')

	if not np.isclose(dt1, dt2):
		raise ValueError(f'dt mismatch: dt1={dt1}, dt2={dt2}.')
	dt = float(dt1)
	if dt <= 0:
		raise ValueError('dt must be positive.')

	s1, e1 = _to_seconds(start1), _to_seconds(end1)
	s2, e2 = _to_seconds(start2), _to_seconds(end2)
	n1_exp = int(round((e1 - s1) / dt)) + 1
	n2_exp = int(round((e2 - s2) / dt)) + 1
	if abs(n1_exp - a.shape[0]) > 1 or abs(n2_exp - b.shape[0]) > 1:
		raise ValueError('Start/end/dt are not consistent with input matrix length.')

	# Snap data2 to data1 grid.
	off = (s2 - s1) / dt
	off_r = int(round(off))
	adjust_sec = abs(off - off_r) * dt
	if adjust_sec > (tolerance * dt):
		raise ValueError(f'Time adjustment too large ({adjust_sec:.6f}s > {tolerance * dt:.6f}s).')

	# Scalar arithmetic only (no large index arrays).
	n1, n2 = a.shape[0], b.shape[0]
	i0 = min(0, off_r)
	a0 = -i0
	b0 = off_r - i0
	n_out = max(a0 + n1, b0 + n2)
	sa = slice(a0, a0 + n1)
	sb = slice(b0, b0 + n2)

	out_dtype = float if (gap in ('nan', 'linear') or overlap == 'mean') else np.result_type(a.dtype, b.dtype)
	shape_out = (n_out, a.shape[1])

	if out is not None:
		if out.shape != shape_out:
			raise ValueError(f'Provided out has wrong shape {out.shape}, expected {shape_out}.')
		if not np.can_cast(out_dtype, out.dtype, casting='safe'):
			raise ValueError(f'Provided out dtype {out.dtype} is incompatible with required {out_dtype}.')
		out_arr = out
		out_arr[:] = np.nan if gap in ('nan', 'linear') else 0
	else:
		if gap in ('nan', 'linear'):
			out_arr = np.full(shape_out, np.nan, dtype=out_dtype)
		else:
			out_arr = np.zeros(shape_out, dtype=out_dtype)

	if overlap == 'data2':
		out_arr[sa, :] = a
		out_arr[sb, :] = b
	elif overlap == 'data1':
		out_arr[sb, :] = b
		out_arr[sa, :] = a
	else:  # overlap == 'mean'
		out_arr[sa, :] = a
		out_arr[sb, :] = b
		ov0 = max(a0, b0)
		ov1 = min(a0 + n1, b0 + n2)
		if ov1 > ov0:
			ia = ov0 - a0
			ib = ov0 - b0
			n_ov = ov1 - ov0
			out_arr[ov0:ov1, :] = 0.5 * (a[ia:ia+n_ov, :] + b[ib:ib+n_ov, :])

	if gap == 'linear':
		_interp_nan_along_axis0_inplace(out_arr)

	return np.moveaxis(out_arr, 0, axis)



def _to_seconds(value) -> float:
    
    if isinstance(value, (int, float, np.integer, np.floating)):
        
        return float(value)
    
    return float(UTC(value))


def _interp_nan_along_axis0_inplace(data: np.ndarray) -> np.ndarray:
    
    x = np.arange(data.shape[0])
    
    for j in range(data.shape[1]):
        
        y = data[:, j]
        valid = np.isfinite(y)
        
        if np.any(valid) and not np.all(valid):
            
            y[~valid] = np.interp(x[~valid], x[valid], y[valid])
    
    return data


def time_concatenation(data1: np.ndarray, data2: np.ndarray, start1, end1, dt1: float, start2, end2, 
                      dt2: float, axis: int = 0, overlap: str = 'data2', gap: str = 'nan', tolerance: float = 1.0,
                      out: np.ndarray = None) -> np.ndarray:
	"""2D concatenaiton of matrices with time series data.

	Parameters
	----------
	data1 : np.ndarray
		_description_
	data2 : np.ndarray
		_description_
	start1 : _type_
		_description_
	end1 : _type_
		_description_
	dt1 : float
		_description_
	start2 : _type_
		_description_
	end2 : _type_
		_description_
	dt2 : float
		_description_
	axis : int, optional
		_description_, by default 0
	overlap : str, optional
		_description_, by default 'data2'
	gap : str, optional
		_description_, by default 'nan'
	tolerance : float, optional
		_description_, by default 1.0
	out : np.ndarray, optional
		_description_, by default None

	Returns
	-------
	np.ndarray
		_description_

	Raises
	------
	ValueError
		_description_
	ValueError
		_description_
	ValueError
		_description_
	ValueError
		_description_
	ValueError
		_description_
	ValueError
		_description_
	ValueError
		_description_
	ValueError
		_description_
	ValueError
		_description_
	ValueError
		_description_
	ValueError
		_description_
	"""

	if data1.ndim != 2 or data2.ndim != 2:
		raise ValueError('Both inputs must be 2D matrices.')
	if axis not in (0, 1):
		raise ValueError('axis must be 0 or 1.')
	if overlap not in ('data1', 'data2', 'mean'):
		raise ValueError('overlap must be "data1", "data2", or "mean".')
	if gap not in ('nan', 'zero', 'linear'):
		raise ValueError('gap must be "nan", "zero", or "linear".')

	a = np.moveaxis(data1, axis, 0)
	b = np.moveaxis(data2, axis, 0)
	if a.shape[1] != b.shape[1]:
		raise ValueError('Input shapes are incompatible on non-concatenation axis.')

	if not np.isclose(dt1, dt2):
		raise ValueError(f'dt mismatch: dt1={dt1}, dt2={dt2}.')
	dt = float(dt1)
	if dt <= 0:
		raise ValueError('dt must be positive.')

	s1, e1 = _to_seconds(start1), _to_seconds(end1)
	s2, e2 = _to_seconds(start2), _to_seconds(end2)
	n1_exp = int(round((e1 - s1) / dt)) + 1
	n2_exp = int(round((e2 - s2) / dt)) + 1
	if abs(n1_exp - a.shape[0]) > 1 or abs(n2_exp - b.shape[0]) > 1:
		raise ValueError('Start/end/dt are not consistent with input matrix length.')

	# Snap data2 to data1 grid.
	off = (s2 - s1) / dt
	off_r = int(round(off))
	adjust_sec = abs(off - off_r) * dt
	if adjust_sec > (tolerance * dt):
		raise ValueError(f'Time adjustment too large ({adjust_sec:.6f}s > {tolerance * dt:.6f}s).')

	# Scalar arithmetic only (no large index arrays).
	n1, n2 = a.shape[0], b.shape[0]
	i0 = min(0, off_r)
	a0 = -i0
	b0 = off_r - i0
	n_out = max(a0 + n1, b0 + n2)
	sa = slice(a0, a0 + n1)
	sb = slice(b0, b0 + n2)

	out_dtype = float if (gap in ('nan', 'linear') or overlap == 'mean') else np.result_type(a.dtype, b.dtype)
	shape_out = (n_out, a.shape[1])

	if out is not None:
		if out.shape != shape_out:
			raise ValueError(f'Provided out has wrong shape {out.shape}, expected {shape_out}.')
		if not np.can_cast(out_dtype, out.dtype, casting='safe'):
			raise ValueError(f'Provided out dtype {out.dtype} is incompatible with required {out_dtype}.')
		out_arr = out
		out_arr[:] = np.nan if gap in ('nan', 'linear') else 0
	else:
		if gap in ('nan', 'linear'):
			out_arr = np.full(shape_out, np.nan, dtype=out_dtype)
		else:
			out_arr = np.zeros(shape_out, dtype=out_dtype)

	if overlap == 'data2':
		out_arr[sa, :] = a
		out_arr[sb, :] = b
	elif overlap == 'data1':
		out_arr[sb, :] = b
		out_arr[sa, :] = a
	else:  # overlap == 'mean'
		out_arr[sa, :] = a
		out_arr[sb, :] = b
		ov0 = max(a0, b0)
		ov1 = min(a0 + n1, b0 + n2)
		if ov1 > ov0:
			ia = ov0 - a0
			ib = ov0 - b0
			n_ov = ov1 - ov0
			out_arr[ov0:ov1, :] = 0.5 * (a[ia:ia+n_ov, :] + b[ib:ib+n_ov, :])

	if gap == 'linear':
		_interp_nan_along_axis0_inplace(out_arr)

	return np.moveaxis(out_arr, 0, axis)


def trim_time(t0: UTC | str, tf: UTC | str, data: np.ndarray, times: np.ndarray,
			  start_time: UTC, end_time: UTC) -> tuple[np.ndarray, UTC, UTC]:
	"""Trims data in Fiber class in time dimension between given start and end times

	Parameters
	----------
	t0, tf : UTC datetime object or str
		New desired start and end timea of data.
	data : np.ndarray
		Original data to trim.
	times : np.ndarray
		Timestamps of Fiber class.
	start_time, end_time : UTC datetime
		Original start and end time of data.

	Raises
	------
	ValueError
		``t0 > tf``.

	Returns
	-------
	data : np.ndarray
		Trimmed data.
	start_time, end_time : UTC datetime
		New start and end time of data.

	"""
	t0, tf = UTC(t0), UTC(tf)

	t0 = max(t0, start_time)
	tf = min(tf, end_time)

	if tf < t0: raise ValueError("\n⚠️ End time (tf) must be after start time (t0).")
	t0_pos = max(0, np.searchsorted(times, t0, side="right") - 1)
	tf_pos = max(0, np.searchsorted(times, tf, side="right") - 1)

	data = data[t0_pos:tf_pos, :]
	start_time = times[t0_pos]
	end_time = times[tf_pos]

	return data, start_time, end_time
