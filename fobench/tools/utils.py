
import functools
import inspect
import numpy as np
import h5py as h5
from pyrocko.orthodrome import distance_accurate50m_numpy as deg2m

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
STRAIN_UNIT_MAP = {      # mapping for strain(rate) data
-1: 'integrated strain',
0: 'strain',
1: 'strain-rate',
2: 'strain-acceleration',
3: 'strain-jerk'}

VEL_UNIT_MAP = {         # mapping for velocity data, e.g. Terra15 output
-1: 'm',
0: 'm/s',
1: 'm/s^2',
2: 'm/s^3'}

TEMP_UNIT_MAP = {        # mapping for temperature data
-1: 'integrated temperature',
0: 'temperature',
1: 'temperature rate',
2: 'temperature acceleration'}


def _update_processing(func):
	'''
	Co-authors: Jonas Pätzel
	Description: 
		internal decorator function that updates the Fiber.processing attribute after each processing function
		in case of integration or differentiation of the data also updates Fiber.units
	:Params:
		-func(type: function): preprocessing function that will be logged
	:Return:
		- NA
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
			fs = attributes['o_sampling_frequency'] # sampling frequency, which can be 1000 Hz for raw data.
			factor = i_cst*(fs/gauge_L)/digital_N # strain Rate per counts.
			data = np.multiply(data,factor)
     
            
    # if (format == 'h5' or format == 'hdf5') and company == 'febus': # FEBUS HDF5

	if (format == 'h5' or format == 'hdf5') and company == 'silixa': # Silixa HDF5

		if units == 'counts' and target == 'strain-rate':

			i_cst = 116E-9 # meters per radians.
			gauge_L = attributes['gauge_length'] # gauge lenght in meters.
			digital_N = 2**13 # magic number linked to the digitalization of the data. why not 2**16?
			fs = attributes['o_sampling_frequency'] # sampling frequency, which can be 1000 Hz for raw data.
			factor = i_cst*(fs/gauge_L)/digital_N # strain Rate per counts.
			data = np.multiply(data,factor)
        
    # if format == 'npy' and company == 'bam': # .npy format for BAM. This might fail always since the unit is NON-COMMERCIAL!
    
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
			fs = attributes['o_sampling_frequency'] # sampling frequency, which can be 1000 Hz for raw data.
			factor = i_cst*(fs/gauge_L)/digital_N # strain Rate per counts.
			data = np.multiply(data,factor)

	# ####################################################
	# CAUTION!! NON OFFICIAL / EXPERIMENTAL FORMATS, ONLY FOR SPECIAL CASES.
	# ####################################################
   
	if format == 'npz' and company == 'bam': # .npy format for BAM. This might fail always since the unit is NON-COMMERCIAL!
    
		if units == 'counts' and target == 'strain':
        
			factor = 1E-6 / 18.4 # strain per count (WEIRD!)
			data = np.multiply(data,factor)

	if (format == 'h5' or format == 'hdf5') and company == 'michelle': # Michelle HDF5 decimated from Silixa

		if units == 'counts' and target == 'strain-rate':

			i_cst = 116E-9 # meters per radians.
			gauge_L = attributes['gauge_length'] # gauge lenght in meters.
			digital_N = 2**13 # magic number linked to the digitalization of the data. why not 2**16?
			fs = attributes['o_sampling_frequency'] # sampling frequency, which can be 1000 Hz for raw data.
			factor = i_cst*(fs/gauge_L)/digital_N # strain Rate per counts.
			data = np.multiply(data,factor)
		
	return data, target, attributes['channels'], attributes['channels_num'], attributes['total_channels']


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
