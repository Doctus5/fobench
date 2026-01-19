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
from tqdm import tqdm

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
        - 'hfd5':'terra15'
        - 'hdf5':'asn'
        - 'npy':'bam' (non-commercial. Errors can be present during reading).
        - 'hdf5':'quantx'
        - 'h5' : 'aragon'
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

    elif isinstance(range_ch, np.ndarray): # check if is array
        range_ch = list(range_ch)

    if format == 'tdms' and company == 'silixa': # Silixa TDMS

        pbar = tqdm(total=1, leave=True, desc='Reading Silixa TDMS file')
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

        pbar = tqdm(total=1, leave=True, desc='Reading Silixa HDF5 file')
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
        data = __data__(dataset, format, company, list(chans_nums)) if load_data == True else None
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

        pbar = tqdm(total=1, leave=True, desc='Reading Febus HDF5 file')
        file_file = h5.File(filepath,'r')
        instrument = list(file_file.keys())[0]
        properties = file_file[instrument]['Source1']['Zone1'].attrs
        fiber = 'febus' # VARIABLE TO CHANGE FOR REAL NAME EXTRACTION
        measure_type = list(file_file[instrument]['Source1']['Zone1'].keys())[0]
        dataset = file_file[instrument]['Source1']['Zone1'][measure_type]
        chans_nums = [i for i in range(dataset.shape[2])] if not range_ch else range_ch
        chans = np.array(chans_nums)
        # loading the data conditioned
        #LAG = properties['BlockOverlap'][0] # 201 is standard. How much the data is repeated in batches.
        sampling_frequency = 1/(properties['Spacing'][1]*1e-3)
        o_sampling_frequency = sampling_frequency
        dt = 1/sampling_frequency
        SEMILAG = properties['BlockRate'][0]*1e-3 if getattr(properties['BlockRate'], "size", 0) == 1 else properties['BlockRate']*1e-3 # unpacking value sif is inside list.
        SEMILAG = int(np.round((1/SEMILAG) / dt))
        data = __data__(dataset, format, company, list(chans_nums), SEMILAG) if load_data == True else None
        start_time = UTC(file_file[instrument]['Source1']['time'][0]) #Time is retaken to correct for LAG.
        end_time = UTC(file_file[instrument]['Source1']['time'][-1]) + (dt * SEMILAG)
        spatial_interval = properties['Spacing'][0]
        time_length = end_time - start_time
        num_points = int(time_length/dt) # check!
        gauge_length = properties['GaugeLength'][0] # CHECK THIS!!
        channel_offset = 0 # FIX THIS!!
        units = 'counts'
        conv_factor = None # conversion factor if given explicitly

    elif (format == 'h5' or format == 'hdf5') and company == 'terra15': # Terra15 HDF5

        pbar = tqdm(total=1, leave=True, desc='Reading Terra15 HDF5 file')
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

        pbar = tqdm(total=1, leave=True, desc='Reading ASN HDF5 file')
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

    elif (format == 'h5' or format == 'hdf5') and company == 'quantx': # QuantX OptaSense HDF5

        pbar = tqdm(total=1, leave=True, desc='Reading QuantX HDF5 file')
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

    elif (format == 'h5' or format == 'hdf5') and company == 'aragon': # Aragon Photonics HDAS HDF5

        pbar = tqdm(total=1, leave=True, desc='Reading Aragon Photonics HDF5 file')
        file_file = h5.File(filepath,'r')
        properties = file_file.attrs
        dataset = file_file['strain']
        chans_nums = list(range(file_file['position'].size)) if range_ch is None else list(range(range_ch[0], range_ch[1] + 1))
        chans = np.array(chans_nums)
        # loading the data conditioned
        data = __data__(dataset, format, company, list(chans_nums)) if load_data == True else None
        pbar.set_description('Extracting Attributes')
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
        units = [key for key in file_file.keys()][2]
        conv_factor = None # conversion factor if given explicitly

    elif (format == 'h5' or format == 'hdf5') and company == 'sintela': # Sintela Onyx HDF%
        print('Data units (strain, strain-rate...) can not be extracted from Sintela files\n'
             'you can set it manually by editing Fiber.units')
        pbar = tqdm(total=1, leave=True, desc='Reading Sintela HDF5 file')
        file_file = h5.File(filepath, 'r')
        properties = file_file['Acquisition'].attrs
        dataset = file_file['Acquisition/Raw[0]/RawData']
        chans_nums = chans_nums = list(range(properties['NumberOfLoci'])) if range_ch == None else list(range(range_ch[0], range_ch[1] + 1))
        chans = np.array(chans_nums)
        data = __data__(dataset, format, company, list(chans_nums)) if load_data == True else None
        fiber = properties['FiberID']
        sampling_frequency = properties['PulseRate']
        o_sampling_frequency = sampling_frequency
        num_points = dataset.shape[0]
        start_time = UTC(properties['MeasurementStartTime'].decode('utf-8'))
        dt = 1 / sampling_frequency
        end_time = UTC(start_time + num_points * dt)
        time_length = end_time - start_time
        spatial_interval = float(properties['SpatialSamplingInterval'])
        channel_offset = properties['StartLocusIndex']*spatial_interval
        gauge_length = float(properties['GaugeLength'])
        units = 'units' # unit is not coded in metadata for Sintela .h5
        conv_factor = None

	# ####################################################
	# CAUTION!! NON OFFICIAL / EXPERIMENTAL FORMATS, ONLY FOR SPECIAL CASES.
	# ####################################################

    elif format == 'npy' and company == 'bam': # .npy format for BAM. This might fail always since the unit is NON-COMMERCIAL!

        print('File format is a Numpy Class. It contains only the unitsdata, and so the metadata must be filled automatically in the code.')
        file_file = None
        properties = None
        chans = None
        dataset = np.load(filepath)
        chans_nums = [i for i in range(dataset.shape[1])] if not range_ch else range_ch
        chans = np.array(chans_nums)
        # loading the data conditioned
        data = __data__(filepath, format, company, list(chans_nums)) if load_data == True else None
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

        print('File format is a Numpy Zip Class. No Gauge Length specified. Do not attempt to convert to Strain-Rate.')
        file_file = None
        properties = None
        chans = None
        dataset = np.load(filepath)
        chans_nums = [i for i in range(len(dataset['distance']))] if not range_ch else range_ch
        chans = np.array(chans_nums)
        # loading the data conditioned
        data = __data__(filepath, format, company, list(chans_nums)) if load_data == True else None
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
        raise ValueError(f'{format} is not a recognized file format!')
    del dataset

    # Attributed for the Fiber class.
    attr_keys = [
        'basefile',
        'format',
        'company',
        'fiber',
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
        ]

    attributes = [filepath,
                format,
                company,
                fiber,
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
                ]

    attributes = dict(zip(attr_keys,attributes))
    pbar.update(1)
    pbar.set_description('Read File ✓')
    # pbar.close()
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
        values = extract_point.as_dataframe().to_numpy()[:, range_ch]

    elif format in ('h5', 'hdf5'):
        if company == 'febus':
            dims = extract_point.shape
            values = extract_point[:, :LAG, :].reshape(dims[0] * LAG, dims[2])[:, range_ch[0]:range_ch[1]]
        elif company == 'terra15':
            values = np.array(extract_point[:, range_ch[0]:range_ch[1]])
        elif company in ('silixa', 'asn', 'quantx', 'aragon', 'sintela'):
            values = np.array(extract_point[:, range_ch])
            if company == 'aragon':
                values *= 1e-9  # Convert from nanostrain to strain

    # ####################################################
    # CAUTION!! NON OFFICIAL / EXPERIMENTAL FORMATS, ONLY FOR SPECIAL CASES.
    # ####################################################

    elif format == 'npy' and company == 'bam':
        values = np.load(extract_point)[:, range_ch]

    elif format == 'npz' and company == 'bam':
        values = np.load(extract_point)['ph'][:, range_ch]

    elif (format == 'h5' or format == 'hdf5') and company == 'michelle':
        values = np.array(extract_point[:, range_ch])

    return values.astype('float')