"""
File managing the different fiber optic sensing formats.

Created on 28-01-2024 12:07:17
Last modification on 2023-09-14 14:51:00

:author:
	- Sergio Diaz (sergioad@gfz-potsdam.de)
:contributors:
	- Christopher Wollin (wollin@gfz-potsdam.de)
:license:

"""

# Necessary packes to read the file formats
import nptdms as tdms
import logging
from nptdms.log import log_manager
log_manager.set_level(logging.ERROR)
import h5py as h5
import numpy as np

from collections import namedtuple
from obspy.core import UTCDateTime as UTC



def read_data(filepath=None, company=None, range_ch=None, format=None, load_data=True):
    '''
    Co-authors: --
    Description:
        Manages the import of the data in the best way it can be done according to the format.
        Format and company availables (format, company)...
        - 'tmds':'silixa'
        - 'hdf5':'silixa'
        - 'hdf5':'febus'
        - 'hfd5':'tera15'
        - 'hdf5':'asn'
        - 'npy':'bam' (non-commercial. Errors can be present during reading).
        - 'hdf5':'quantx'
    :Params:
        - filepath(type:String): compelte path fot he file to be read.
        - company(type:String): manufacturer or the instrument that generates the data.
        - range_ch(type:Int or List): channel number(s) to load only in data. Method to avoid loading all the data (Not teste for other than Silixa).
    :Return:
        - variables(type:tuple): tuple of the variable that are attributed of class Fiber.  
    '''

    # modify range_ch variable
    if isinstance(range_ch, int): # check if it is single value
        
        range_ch = [range_ch]
        
    if isinstance(range_ch, np.ndarray): # check if is array
        
        range_ch = list(range_ch)
    
    
    # Start finding the attributes...
    if format == 'tdms' and company == 'silixa': # Silixa TDMS

        print('Reading TDMS file (Silixa Format)...')
        file_file = tdms.TdmsFile.read(filepath)
        properties = file_file.properties
        dataset = None
        chans = file_file['Measurement'].channels() if range_ch == None else [file_file['Measurement'].channels()[i] for i in range_ch]
        chans_nums = [int(chan.name) for chan in chans]
        # loading of the data conditioned.
        data = __data__(file_file['Measurement'], format, company, chans_nums) if load_data == True else None
        fiber = properties['name'].split('_')[0]
        sampling_frequency = properties['SamplingFrequency[Hz]']
        o_sampling_frequency = sampling_frequency
        dt = 1/sampling_frequency
        start_time = UTC(properties['ISO8601 Timestamp'])
        end_time = UTC(start_time + (len(chans[0])-1)*dt)
        spatial_interval = properties['SpatialResolution[m]']
        time_length = end_time - start_time
        num_points = len(chans[0])
        gauge_length = properties['GaugeLength']
        channel_offset = properties['OffsetLength']
        units = 'counts'
        conv_factor = None # conversion factor if given explicitly
        
    elif (format == 'h5' or format == 'hdf5') and company == 'silixa': # Silixa HDF5
    
        print('Reading H5 file (Silixa Format)...')
        file_file = h5.File(filepath,'r')
        properties = {}
        instrument = file_file.keys()
        dataset = file_file['Acquisition']['Raw[0]']['RawData']

        # load of all hidden properties
        for vv in list(file_file['Acquisition'].attrs): properties[vv] = file_file['Acquisition'].attrs[vv]  #mandatory!!
        for vv in list(file_file['Acquisition']['Custom'].attrs): properties[vv] = file_file['Acquisition']['Custom'].attrs[vv]  #optional
        for vv in list(file_file['Acquisition']['Custom']['AdvancedUserSettings'].attrs): properties[vv] = file_file['Acquisition']['Custom']['AdvancedUserSettings'].attrs[vv]  #mandatory!!
        for vv in list(file_file['Acquisition']['Custom']['SystemInformation']['Chassis'].attrs): properties[vv] = file_file['Acquisition']['Custom']['SystemInformation']['Chassis'].attrs[vv]  #optional
        for vv in list(file_file['Acquisition']['Custom']['SystemInformation']['Devices0'].attrs): properties[vv] = file_file['Acquisition']['Custom']['SystemInformation']['Devices0'].attrs[vv]  #optional
        for vv in list(file_file['Acquisition']['Custom']['SystemInformation']['Devices1'].attrs): properties[vv] = file_file['Acquisition']['Custom']['SystemInformation']['Devices1'].attrs[vv]  #optional
        for vv in list(file_file['Acquisition']['Custom']['SystemInformation']['Devices2'].attrs): properties[vv] = file_file['Acquisition']['Custom']['SystemInformation']['Devices2'].attrs[vv]  #optional
        for vv in list(file_file['Acquisition']['Custom']['SystemInformation']['GPS'].attrs): properties[vv] = file_file['Acquisition']['Custom']['SystemInformation']['GPS'].attrs[vv]  #mandatory
        for vv in list(file_file['Acquisition']['Custom']['SystemInformation']['OSVersion'].attrs): properties[vv] = file_file['Acquisition']['Custom']['SystemInformation']['OSVersion'].attrs[vv]  #optional/mandatory
        for vv in list(file_file['Acquisition']['Custom']['SystemInformation']['ProcessingUnit'].attrs): properties[vv] = file_file['Acquisition']['Custom']['SystemInformation']['ProcessingUnit'].attrs[vv]  #optional
        for vv in list(file_file['Acquisition']['Custom']['SystemSettings'].attrs): properties[vv] = file_file['Acquisition']['Custom']['SystemSettings'].attrs[vv]  #mandatory!!
        for vv in list(file_file['Acquisition']['Custom']['UserSettings'].attrs): properties[vv] = file_file['Acquisition']['Custom']['UserSettings'].attrs[vv]  #mandatory!!
        for vv in list(file_file['Acquisition']['Raw[0]'].attrs): properties[vv] = file_file['Acquisition']['Raw[0]'].attrs[vv]  #mandatory!!
        for vv in list(file_file['Acquisition']['Raw[0]']['RawData'].attrs): properties[vv] = file_file['Acquisition']['Raw[0]']['RawData'].attrs[vv]  #optional
        for vv in list(file_file['Acquisition']['Raw[0]']['RawDataTime'].attrs): properties[vv] = file_file['Acquisition']['Raw[0]']['RawDataTime'].attrs[vv]  #optional

        chans = [i for i in range(properties['NumberOfLoci'])] if not range_ch else range_ch
        chans_nums = np.array(chans)
        # loading the data conditioned
        data = __data__(dataset, format, company, list(chans_nums[range_ch])) if load_data == True else None
        fiber = properties['FibreType']
        sampling_frequency = properties['OutputDataRate']
        o_sampling_frequency = sampling_frequency
        dt = 1/sampling_frequency
        start_time = UTC(properties['PartStartTime'])
        end_time = UTC(start_time + (properties['Count']-1)*dt)
        spatial_interval = properties['SpatialResolution']
        time_length = end_time - start_time
        num_points = properties['Count']
        gauge_length = properties['GaugeLength']
        channel_offset = abs(properties['PreTriggerSamples'])
        units = 'counts'
        conv_factor = None # conversion factor if given explicitly
    
    # Required checking! Contact providers! TEST!
    elif (format == 'h5' or format == 'hdf5') and company == 'febus': # FEBUS HDF5
    
        print('Reading H5 file (Febus Format)...')
        LAG = 201 #important parameter! Sometimes the data is repeated in batches. This number indicates the position in the minibatch where the data begins to be repeated.
        file_file = h5.File(filepath,'r')
        instrument = list(file_file.keys())[0]
        properties = file_file[instrument]['Source1']['Zone1'].attrs
        fiber = 'febus' # VARIABLE TO CHANGE FOR REAL NAME EXTRACTION
        meassure_type = list(file_file[instrument]['Source1']['Zone1'].keys())[0]
        dataset = file_file[instrument]['Source1']['Zone1'][meassure_type]
        chans_nums = [i for i in range(dataset.shape[2])]
        chans = np.array(chans_nums)
        # loading the data conditioned
        data = __data__(dataset, format, company, list(chans_nums[range_ch]), LAG) if load_data == True else None
        sampling_frequency = 1/(properties['Spacing'][1]*1e-3)
        o_sampling_frequency = sampling_frequency
        dt = 1/sampling_frequency
        start_time = UTC(file_file[instrument]['Source1']['time'][0]) #Time is retaken to correct for LAG.
        end_time = UTC(file_file[instrument]['Source1']['time'][-1]) + (dt * LAG)
        spatial_interval = properties['Spacing'][0]
        time_length = end_time - start_time
        num_points = int(time_length/dt) # check!
        gauge_length = properties['GaugeLength'] # CHECK THIS!!
        channel_offset = 0 # FIX THIS!!
        units = 'counts'
        conv_factor = None # conversion factor if given explicitly
    
    #file_file.close()
    
    elif (format == 'h5' or format == 'hdf5') and company == 'terra15': # Terra15 HDF5
        
        print('Reading HDF5 file (Terra15 Format)...')
        file_file = h5.File(filepath,'r')
        properties = file_file.attrs
        dataset = file_file['data_product']
        chans_nums = [i for i in range(properties['nx'])] if range_ch == None else range_ch
        chans = np.array(chans_nums)
        # loading the data conditioned
        data = __data__(dataset['data'], format, company, list(chans_nums)) if load_data == True else None
        fiber = 'standard'
        dt = float(properties['dt_computer'])
        sampling_frequency = 1 / dt
        o_sampling_frequency = sampling_frequency
        num_points = int(properties['nt'])
        print(properties['nt'])
        start_time = UTC(properties['file_start_gps_time']) if properties['file_start_gps_time'] else UTC(properties['file_start_computer_time'])
        end_time = UTC(start_time + num_points * dt)
        spatial_interval = float(properties['dx'])
        time_length = end_time - start_time
        gauge_length = float(properties['gauge_length'])
        channel_offset = int(properties['sensing_range_start'] / spatial_interval)
        units = properties['data_product_units']
        conv_factor = None # conversion factor if given explicitly
        
    elif (format == 'h5' or format == 'hdf5') and company == 'asn': # ASN OptoDAS HDF5 (It can be a bit more complex, so I'm trying to make it simple!)
        
        print('Reading HDF5 file (ASN Format)...')
        file_file = h5.File(filepath,'r')
        properties = file_file['acqSpec']
        dataset = file_file['data']
        chans_nums = [i for i in range(int(file_file['header']['dimensionRanges']['dimension1']['size'][()]))] if range_ch == None else range_ch
        chans = np.array(chans_nums)
        # loading the data conditioned
        data = __data__(dataset, format, company, list(chans_nums)) if load_data == True else None
        fiber = 'standard'
        dt = float(file_file['header']['dt'][()])
        sampling_frequency = 1 / dt
        o_sampling_frequency = sampling_frequency
        num_points = int(file_file['header']['dimensionRanges']['dimension0']['size'][()])
        start_time = UTC(float(file_file['header']['time'][()]))
        end_time = UTC(start_time + num_points * dt)
        original_channels = np.array(file_file['header']['channels'])
        spatial_interval = float(file_file['header']['dx'][()]) * (original_channels[1] - original_channels[0])
        time_length = end_time - start_time
        gauge_length = float(file_file['header']['gaugeLength'][()])
        channel_offset = int(original_channels[0])
        units = str(file_file['header']['sensitivityUnits'][()])[3:-2]
        conv_factor = file_file['header']['sensitivities'][0,0]
        
    elif (format == 'h5' or format == 'hdf5') and company == 'quantx': # QuantX OptoaSense HDF5
        
        print('Reading HDF5 file (QuantX Format)...')
        file_file = h5.File(filepath,'r')
        properties = file_file['Acquisition'].attrs
        dataset = file_file['Acquisition']['Raw[0]']
        chans_nums = [i for i in range(int(file_file['Acquisition']['Raw[0]'].attrs['NumberOfLoci']))] if range_ch == None else range_ch
        chans = np.array(chans_nums)
        # loading the data conditioned
        data = __data__(dataset['RawData'], format, company, list(chans_nums)) if load_data == True else None
        fiber = 'standard'
        sampling_frequency = float(dataset.attrs['OutputDataRate'])
        o_sampling_frequency = sampling_frequency
        dt = 1 / sampling_frequency
        num_points = int(file_file['Acquisition']['Raw[0]']['RawDataTime'].attrs['Count'])
        start_time = UTC(str(properties['MeasurementStartTime'])[2:-1])
        end_time = UTC(start_time + num_points * dt)
        spatial_interval = float(properties['SpatialSamplingInterval'])
        time_length = end_time - start_time
        gauge_length = float(properties['GaugeLength'])
        channel_offset = int(properties['StartLocusIndex'])
        units = str(dataset.attrs['RawDataUnit'])[2:-1]
        conv_factor = None # conversion factor if given explicitly
    
    # UNDER CONSTRUCTION   
    elif (format == 'h5' or format == 'hdf5') and company == 'aragon': # Aragon Photonics HDAS HDF5
        
        print('Reading HDF5 file (Aragon Photonics Format)...')
        file_file = h5.File(filepath,'r')
        properties = file_file.attrs
        dataset = file_file['strain']
        chans_nums = [i for i in range(file_file['position'].size)] if range_ch == None else range_ch
        chans = np.array(chans_nums)
        # loading the data conditioned
        data = __data__(dataset, format, company, list(chans_nums)) if load_data == True else None
        fiber = 'standard'
        sampling_frequency = float(properties['trigger_frequency'][0])
        o_sampling_frequency = sampling_frequency
        dt = 1 / sampling_frequency
        num_points = int(file_file['time'].size)
        start_time = UTC(properties['Global_FileSaveIO_TimeStamp_s'][0])
        end_time = UTC(start_time + num_points * dt)
        spatial_interval = float(properties['spatial_sampling'][0])
        time_length = end_time - start_time
        gauge_length = float(properties['Global_RAM_User_SET_Pulse_Width_(meter)'][0])
        channel_offset = int(properties['fiber_position_offset'][0]/spatial_interval)
        units = [key for key in file_file.keys()][1]
        conv_factor = None # conversion factor if given explicitly
        
	# ####################################################
	# CAUTION!! NON OFFICIAL / EXPERIMENTAL FORMATS, ONLY FOR SPECIAL CASES.
	# ####################################################
    
    elif format == 'npy' and company == 'bam': # .npy format for BAM. This might fail always since the unit is NON-COMMERCIAL!
    
        print('File format is a Numpy Class. It contains only the unitsdata, and so the metadata must be filled automatically in the code.')
        file_file = None
        properties = None
        chans = None
        dataset = np.load(filepath)
        chans_nums = [i for i in range(dataset.shape[1])]
        chans = np.array(chans_nums)
        # loading the data conditioned
        data = __data__(filepath, format, company, list(chans_nums[range_ch])) if load_data == True else None
        fiber = 'La Chida'
        sampling_frequency = 100000
        o_sampling_frequency = sampling_frequency
        dt = 1/sampling_frequency
        start_time = UTC('2023-03-01T00:00:00')
        spatial_interval = 0.4
        num_points = int(dataset.shape[0])
        end_time = UTC(start_time + num_points*dt)
        time_length = end_time - start_time
        gauge_length = 2
        channel_offset = None
        units = None
        conv_factor = None # conversion factor if given explicitly

    elif format == 'npz' and company == 'bam': # .npy format for BAM. This might fail always since the unit is NON-COMMERCIAL!
    
        print('File format is a Numpy Zip Class. No Gaueg Length specified. Do not attempt to convert to Stran-Rate.')
        file_file = None
        properties = None
        chans = None
        dataset = np.load(filepath)
        chans_nums = [i for i in range(len(dataset['distance']))]
        chans = np.array(chans_nums)
        # loading the data conditioned
        data = __data__(filepath, format, company, list(chans_nums[range_ch])) if load_data == True else None
        fiber = 'La Chida'
        sampling_frequency = dataset['freq']
        o_sampling_frequency = sampling_frequency
        dt = 1/sampling_frequency
        start_time = UTC('2023-03-01T00:00:00')
        spatial_interval = dataset['distance'][1] - dataset['distance'][0]
        num_points = int(dataset['time'][:-2].shape[0])
        end_time = UTC(start_time + num_points*dt)
        time_length = end_time - start_time
        gauge_length = 2
        channel_offset = None
        units = None
        conv_factor = None # conversion factor if given explicitly
        
    elif (format == 'h5' or format == 'hdf5') and company == 'michelle': # .h5 format for Michelle INGV's decimated files from Silixa.

        print('Reading HDF5 file (Michelle INGV decimated Format)...')
        file_file = h5.File(filepath,'r')
        properties = file_file.attrs
        dataset = file_file['Fiber']
        chans_nums = list(file_file['ChannelMap']) if range_ch == None else range_ch
        chans = np.array(chans_nums)
        # loading the data conditioned
        data = __data__(dataset, format, company, list(chans_nums)) if load_data == True else None
        fiber = properties['Fibre Type'].decode('UTF-8')
        dt = float(properties['Sampletime'][0])
        sampling_frequency = 1 / dt
        o_sampling_frequency = properties['SamplingFrequency[Hz]'][0]
        num_points = int(dataset.shape[0])
        start_time = UTC(properties['StartTime_txt'].decode('UTF-8')) if properties['StartTime_txt'] else UTC(properties['CPUTimeStamp_txt'].decode('UTF-8'))
        end_time = UTC(start_time + (num_points-1) * dt)
        spatial_interval = float(properties['SpatialResolution[m]'][0])
        time_length = end_time - start_time
        gauge_length = float(properties['GaugeLength'][0])
        channel_offset = int(properties['OffsetLength'])
        units = 'counts'
        conv_factor = None # conversion factor if given explicitly 

    else:
        
        # Terminate if file format can not be handled.
        raise ValueError('"'+format+'" is not a recognized file format.')
    
    # Free space by erasing old content.
    database = None
    file_file = None
    
    # Attributed for the Fiber class.
    result_tuple = namedtuple('attributes',[
        'file',
        'fiber',
        'dataset',
        'properties',
        'chans',
        'chans_nums',
        'list_chans_num',
        'sampling_frequency',
        'o_sampling_frequency',
        'dt',
        'start_time',
        'end_time',
        'spatial_interval',
        'num_points',
        'time_length',
        'gauge_length',
        'channel_offset',
        'data',
        'units',
        'conv_factor'
    ])
    
    attributes = result_tuple(
                file_file, 
                fiber,
                dataset, 
                h5_to_dict(properties), 
                chans, 
                chans_nums, 
                len(chans_nums),
                sampling_frequency,
                o_sampling_frequency,
                dt,
                start_time,
                end_time,
                spatial_interval,
                num_points,
                time_length,
                gauge_length,
                channel_offset,
                data,
                units,
                conv_factor
                )

    return attributes


# Recurive method to convert all h5py Objects into dictionaries.
def h5_to_dict(h5_obj):
    
    result = {key: h5_to_dict(item) if isinstance(item, h5.Group) else item[()] if isinstance(item, h5.Dataset) else item for key, item in h5_obj.items()}
    
    return result


#Loads the data of the tdms file into a numpy array. Axis 0 is the time, and axis 1 are the channels.
def __data__(extract_point, format, company, range_ch, LAG=None):
    '''
    Co-authors: --
    Description: 
        Extracts the data depending on the file type. This is done automatically during initialization of the class.
    :Params:
        - NA.
    :Return:
        - values(type:Numpy): 2D numpy matrix with values in time per channel. Axis 0 (rows) is time and axis 1 (columns) are the channels.  
    '''

    values = np.asarray([])
    
    if format == 'tdms' and company == 'silixa':
        
        values = extract_point.as_dataframe().to_numpy()[:,range_ch] # New way to load data. Cuts time by half.
        
    if (format == 'h5' or format == 'hdf5') and company == 'febus':
    
        dims = extract_point.shape
        values = extract_point[:,:LAG,:].reshape(int(dims[0]*LAG),dims[2])[:,range_ch]

    if (format == 'h5' or format == 'hdf5') and company == 'silixa':

        values = np.array(extract_point[:,range_ch])

    if (format == 'h5' or format == 'hdf5') and company == 'terra15':

        values = np.array(extract_point[:,range_ch])

    if (format == 'h5' or format == 'hdf5') and company == 'asn':

        values = np.array(extract_point[:,range_ch])

    if (format == 'h5' or format == 'hdf5') and company == 'quantx':

        values = np.array(extract_point[:,range_ch])
    
    if (format == 'h5' or format == 'hdf5') and company == 'aragon':

        values = np.array(extract_point[:,range_ch])
        values *= (10**-9) # from nanostrain to strain.
        
	# ####################################################
	# CAUTION!! NON OFFICIAL / EXPERIMENTAL FORMATS, ONLY FOR SPECIAL CASES.
	# ####################################################
        
    if format == 'npy' and company == 'bam':
    
        values = np.load(extract_point)[:,range_ch]

    if format == 'npz' and company == 'bam':
    
        values = np.load(extract_point)['ph'][:,range_ch]
        
    if (format == 'h5' or format == 'hdf5') and company == 'michelle':

        values = np.array(extract_point[:,range_ch])

    return values.astype('float')