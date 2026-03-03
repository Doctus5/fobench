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
import copy
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
    
    template = None

    if format == 'tdms' and company == 'silixa': # Silixa TDMS

        pbar = tqdm(total=1, leave=True, desc='Reading Silixa TDMS file')
        file_file = tdms.TdmsFile.read(filepath)
        template = __scan_template__(core_file=file_file, company=company, format=format)
        properties = file_file.properties
        dataset = None
        measurement = file_file['Measurement']
        chans = measurement.channels() if range_ch == None else [measurement.channels()[i] for i in range_ch]
        chans_nums = [int(chan.name) for chan in chans]
        # loading of the data conditioned.
        data = __data__(measurement, format, company, chans_nums) if load_data else None
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
        
        with h5.File(filepath,'r') as file_file:
            
            template = __scan_template__(core_file=file_file, company=company, format=format)
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
            data = __data__(dataset, format, company, list(chans_nums)) if load_data else None
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
        
        with h5.File(filepath,'r') as file_file:
            
            instrument = list(file_file.keys())[0]
            properties = dict(file_file[instrument]['Source1']['Zone1'].attrs)
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
            data = __data__(dataset, format, company, list(chans_nums), SEMILAG) if load_data else None
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
        
        with h5.File(filepath,'r') as file_file:
            
            properties = dict(file_file.attrs)
            dataset = file_file['data_product']
            chans_nums = [i for i in range(properties['nx'])] if range_ch == None else range_ch
            chans = np.array(chans_nums)
            # loading the data conditioned
            data = __data__(dataset['data'], format, company, list(chans_nums)) if load_data else None
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
        
        with h5.File(filepath,'r') as file_file:
            
            properties = h5_to_dict(file_file['acqSpec'])
            dataset = file_file['data']
            chans_nums = [i for i in range(int(file_file['header']['dimensionRanges']['dimension1']['size'][()]))] if range_ch == None else range_ch
            chans = np.array(chans_nums)
            # loading the data conditioned
            data = __data__(dataset, format, company, list(chans_nums)) if load_data else None
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
        
        with h5.File(filepath,'r') as file_file:
            
            properties = dict(file_file['Acquisition'].attrs)
            dataset = file_file['Acquisition']['Raw[0]']
            chans_nums = [i for i in range(int(file_file['Acquisition']['Raw[0]'].attrs['NumberOfLoci']))] if range_ch == None else range_ch
            chans = np.array(chans_nums)
            # loading the data conditioned
            data = __data__(dataset['RawData'], format, company, list(chans_nums)) if load_data else None
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
        
        with h5.File(filepath,'r') as file_file:
            
            template = __scan_template__(core_file=file_file, company=company, format=format)
            properties = dict(file_file.attrs)
            dataset = file_file['strain']
            chans_nums = list(range(file_file['position'].size)) if range_ch is None else list(range(range_ch[0], range_ch[1] + 1))
            chans = np.array(chans_nums)
            # loading the data conditioned
            data = __data__(dataset, format, company, list(chans_nums)) if load_data else None
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

    elif (format == 'h5' or format == 'hdf5') and company == 'sintela': # Sintela Onyx HDF5
        print('Data units (strain, strain-rate...) can not be extracted from Sintela files\n'
             'you can set it manually by editing Fiber.units')
        pbar = tqdm(total=1, leave=True, desc='Reading Sintela HDF5 file')
        
        with h5.File(filepath, 'r') as file_file:
            
            template = __scan_template__(core_file=file_file, company=company, format=format)
            properties = dict(file_file['Acquisition'].attrs)
            dataset = file_file['/Acquisition/Raw[0]/RawData']
            chans_nums = chans_nums = list(range(properties['NumberOfLoci'])) if range_ch == None else list(range(range_ch[0], range_ch[1] + 1))
            chans = np.array(chans_nums)
            data = __data__(dataset, format, company, list(chans_nums)) if load_data else None
            fiber = properties['FiberID']
            sampling_frequency = properties['PulseRate']
            o_sampling_frequency = sampling_frequency
            num_points = dataset.shape[0]
            start_time = UTC(properties['MeasurementStartTime'].decode('utf-8')) if isinstance(properties['MeasurementStartTime'],bytes) else UTC(properties['MeasurementStartTime'])
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
        data = __data__(filepath, format, company, list(chans_nums)) if load_data else None
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
        data = __data__(filepath, format, company, list(chans_nums)) if load_data else None
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

        pbar = tqdm(total=1, leave=True, desc='Reading Michelle INGV decimated HDF5 file')
        
        with h5.File(filepath, 'r') as file_file:
            
            properties = file_file.attrs
            dataset = file_file['Fiber']
            chans_nums = list(file_file['ChannelMap']) if range_ch == None else range_ch
            chans = np.array(chans_nums)
            # loading the data conditioned
            data = __data__(dataset, format, company, list(chans_nums)) if load_data else None
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

    attributes = [template if template is not None else filepath,
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



    



# Save data in original format
def write_data(Fiber, filepath=None, company=None):
    """
    Manages the saving of the data in the original format as it was written from the manufacturer.
    Check format and company availables (format, company) from read_data() function.
    "We can possibly try to give options to save in another format but it will require to known previosuly the formats somehow"

    Parameters
    ----------
    Fiber : fiber class
        _description_
    filepath : str, optional
        compelte path including name to be saved and the format after the dot.
        "path/where/to/save/file.h5", by default None
    company : str, optional
        manufacturer or the instrument that generates the data, by default None

    Raises
    ------
    ValueError
        filepath is required
    """
    
    # warnings
    if filepath is None:
        raise ValueError("filepath is required")
    
    format = filepath.split(".")[-1] # accessing the desired format from the string of the file path.
    
    if format == 'tdms' and company == 'silixa': # Silixa TDMS

        pbar = tqdm(total=1, leave=True, desc='Saving in Silixa TDMS file')
        template = __clone_template__(Fiber.__basefile__)
        root_properties = dict(template["root_properties"])

        # Resolve TDMS group/channels from template structure
        groups = template.get("groups", template.get("group_properties", {}))
        gname = template.get("group_name", "Measurement")
        if gname not in groups and len(groups) > 0:
            gname = next(iter(groups))
        group_info = groups[gname]
        group_properties = dict(group_info.get("properties", {}))
        channels_template = group_info.get("channels", template.get("channels", []))

        # map Fiber attrs back to TDMS root properties
        root_properties["SamplingFrequency[Hz]"] = float(Fiber.sampling_frequency)
        root_properties["ISO8601 Timestamp"] = Fiber.start_time.isoformat()
        root_properties["SpatialResolution[m]"] = float(Fiber.spatial_interval)
        root_properties["GaugeLength"] = float(Fiber.gauge_length)
        root_properties["OffsetLength"] = float(Fiber.channel_offset)

        # keep only current channels (supports crop)
        ch_map = {str(ch["name"]): ch for ch in channels_template}
        selected = [ch_map[str(c)] for c in Fiber.channels_num]

        objects = [tdms.RootObject(properties=root_properties), tdms.GroupObject(gname, properties=group_properties)]
        
        for i, ch in enumerate(selected):
            objects.append(
                tdms.ChannelObject(
                    gname,
                    str(ch["name"]),
                    np.asarray(Fiber.data[:, i]),
                    # np.rint(np.clip(Fiber.data[:, i], np.iinfo(np.int16).min, np.iinfo(np.int16).max)).astype(np.int16), # scales back to int16, buuuut with precision errors if dta gets distorted by processing.
                    properties=dict(ch.get("properties", {})),
                )
            )

        with tdms.TdmsWriter(filepath, mode="w") as w:
            w.write_segment(objects)
            
    if (format == 'h5' or format == 'hdf5') and company == 'sintela': # Sitela H5 files.
        
        pbar = tqdm(total=1, leave=True, desc='Saving in Sintela HDF5 file')
        dataset_path = "/Acquisition/Raw[0]/RawData"
        acquisition_path = "/Acquisition"
        raw_group_path = "/Acquisition/Raw[0]"
        template = __clone_template__(Fiber.__basefile__)

        dataset_node = __find_h5_node_by_path__(template, dataset_path)
        acquisition_node = __find_h5_node_by_path__(template, acquisition_path)
        raw_group_node = __find_h5_node_by_path__(template, raw_group_path)

        # Preserve original type as much as possible for compatibility.
        target_dtype = np.dtype(dataset_node["dtype"])
        output_data = np.asarray(Fiber.data)
        if np.issubdtype(target_dtype, np.integer):
            info = np.iinfo(target_dtype)
            output_data = np.rint(np.clip(output_data, info.min, info.max)).astype(target_dtype)
        else:
            output_data = output_data.astype(target_dtype, copy=False)

        dataset_node["data"] = output_data
        dataset_node["shape"] = output_data.shape

        # Keep chunking only if dimensions still fit new shape.
        chunks = dataset_node.get("chunks")
        if chunks is not None:
            if len(chunks) != len(output_data.shape) or any(ch > sh for ch, sh in zip(chunks, output_data.shape)):
                dataset_node["chunks"] = None

        # Update Sintela metadata from Fiber state.
        acquisition_node.setdefault("attrs", {})
        acquisition_node["attrs"].update({
            "FiberID": Fiber.fiber,
            "PulseRate": Fiber.sampling_frequency,
            "MeasurementStartTime": (Fiber.start_time.isoformat()).encode("utf-8"),
            "NumberOfLoci": int(Fiber.data.shape[1]),
            "SpatialSamplingInterval": Fiber.spatial_interval,
            "StartLocusIndex": int(round(Fiber.channel_offset / Fiber.spatial_interval)),
            "GaugeLength": Fiber.gauge_length
        })

        # If Raw[0] carries channel count, keep it coherent.
        raw_group_node.setdefault("attrs", {})
        if "NumberOfLoci" in raw_group_node["attrs"]:
            raw_group_node["attrs"]["NumberOfLoci"] = int(Fiber.data.shape[1])

        with h5.File(filepath, "w") as f:
            __h5_writter__(template, f)

    if (format == 'h5' or format == 'hdf5') and company == 'aragon': # Aragon H5 files.

        pbar = tqdm(total=1, leave=True, desc='Saving in Aragon HDF5 file')
        template = __clone_template__(Fiber.__basefile__)

        strain_path = "/strain"
        time_path = "/time"
        position_path = "/position"

        strain_node = __find_h5_node_by_path__(template, strain_path)
        time_node = __find_h5_node_by_path__(template, time_path)
        position_node = __find_h5_node_by_path__(template, position_path)

        # Aragon reader converts nanostrain -> strain (x1e-9), so invert on write.
        strain_out = np.asarray(Fiber.data) * 1e9
        strain_dtype = np.dtype(strain_node["dtype"])
        if np.issubdtype(strain_dtype, np.integer):
            info = np.iinfo(strain_dtype)
            strain_out = np.rint(np.clip(strain_out, info.min, info.max)).astype(strain_dtype)
        else:
            strain_out = strain_out.astype(strain_dtype, copy=False)
        strain_node["data"] = strain_out
        strain_node["shape"] = strain_out.shape

        # Keep chunking only if dimensions still fit new shape.
        chunks = strain_node.get("chunks")
        if chunks is not None:
            if len(chunks) != len(strain_out.shape) or any(ch > sh for ch, sh in zip(chunks, strain_out.shape)):
                strain_node["chunks"] = None

        # Rebuild time axis in file units (seconds from start).
        t_dtype = np.dtype(time_node["dtype"])
        t_vals = np.arange(Fiber.data.shape[0], dtype=np.float64) * (1.0 / Fiber.sampling_frequency)
        t_vals = t_vals.astype(t_dtype, copy=False)
        time_node["data"] = t_vals
        time_node["shape"] = t_vals.shape

        # Rebuild position axis in meters from current channel selection.
        if Fiber.channels_num is not None and len(Fiber.channels_num) == Fiber.data.shape[1]:
            p_vals = np.asarray(Fiber.channels_num, dtype=np.float64) * float(Fiber.spatial_interval)
        else:
            p_vals = np.arange(Fiber.data.shape[1], dtype=np.float64) * float(Fiber.spatial_interval)
        p_dtype = np.dtype(position_node["dtype"])
        p_vals = p_vals.astype(p_dtype, copy=False)
        position_node["data"] = p_vals
        position_node["shape"] = p_vals.shape

        # Update known Aragon root attrs while preserving original dtype/shape style.
        root_attrs = template.setdefault("attrs", {})
        __update_h5_attr_like__(root_attrs, "trigger_frequency", float(Fiber.sampling_frequency))
        __update_h5_attr_like__(root_attrs, "spatial_sampling", float(Fiber.spatial_interval))
        __update_h5_attr_like__(root_attrs, "Global_RAM_User_SET_Pulse_Width_(meter)", float(Fiber.gauge_length))
        __update_h5_attr_like__(root_attrs, "fiber_position_offset", float(Fiber.channel_offset) * float(Fiber.spatial_interval))
        __update_h5_attr_like__(root_attrs, "Global_FileSaveIO_TimeStamp_s", float(Fiber.start_time))

        with h5.File(filepath, "w") as f:
            __h5_writter__(template, f)

    pbar.update(1)
    pbar.set_description('File Saved ✓')
    
#########################################################
### H5/HDF5 explorers, templates, internals and writters
#########################################################

def __h5_writter__(template, file):
    """Handles the writting of the H5/HDF5 files.

    Parameters
    ----------
    template : dict
        lightweight template of the h5 structure file.
    file : h5 obj
        opened H5 file to where to save the data.

    Raises
    ------
    ValueError
        _description_
    """

    node_type = template.get("type") or template.get("kind")  # support both keys if needed

    # Determine name in parent
    if template["path"] == "/":  # root
        current_obj = file
    else:
        name = template["path"].split("/")[-1]
        if node_type == "group":
            current_obj = file.create_group(name)
        elif node_type == "dataset":
            data = template.get("data")
            kwargs = {
                "dtype": template["dtype"],
                "chunks": template.get("chunks"),
                "compression": template.get("compression"),
                "maxshape": template.get("maxshape"),
            }
            kwargs = {k: v for k, v in kwargs.items() if v is not None}

            if data is None:
                current_obj = file.create_dataset(
                    name,
                    shape=template["shape"],
                    **kwargs,
                )
            else:
                current_obj = file.create_dataset(
                    name,
                    data=data,
                    **kwargs,
                )
        else:
            raise ValueError(f"Unknown node type at {template['path']}")

    # Set attributes
    for k, v in template.get("attrs", {}).items():
        try:
            current_obj.attrs[k] = v
        except (TypeError, ValueError):
            # Some attrs like cross-file references are not transferable to a new file.
            continue

    # Recurse into children
    for child_name, child_node in template.get("children", {}).items():
        __h5_writter__(child_node, current_obj)
        
        
# scanning the files as saving them as lightweight versions for later file writting.
def __scan_template__(core_file, company=None, format=None):
    """
    Scans the structure of the file and saves it as a lightweight version dictionary.

    Parameters
    ----------
    core_file : file object
        file object once it is opened.
    company : str,
        manufacturer or the instrument that generates the data, by default None
    format : str, optional
        format of the file where the data is coming from, by default None

    Returns
    -------
    dict
        lightweight version of the architecture of the file. Used for later recreate the file and save information in original format.
    """
    
    native_template = None

    if format == 'tdms' and company == 'silixa': # Silixa TDMS
        
        native_template = {
            "root_properties": dict(core_file.properties),  # unchanged
            "group_properties": {
                "Measurement": {
                    "properties": dict(core_file["Measurement"].properties),  # unchanged
                    "channels": [
                        {
                            "name": ch.name,
                            "properties": dict(ch.properties),  # unchanged
                            "dtype": ch.dtype,                  # for writer
                        }
                        for ch in core_file["Measurement"].channels()
                    ],
                }
            },
        }    
    
    elif format in ('h5', 'hdf5'): # evaluate for HDF5 files. Company specification MIGHT not be needed.
        
        native_template = __scan_hdf5__(core_file)

    return native_template


# scans the hdf5 file. It needs to be separated from __scan_template__ because we need recurrence for deper layers.
def __scan_hdf5__(obj):
    """
    Scans a H5/HDF5 paths and structure to a lightweight version. Scan recurrence included for deeper layers.

    Parameters
    ----------
    obj : H5 object
        H5/HDF5 object, usually given after reading it.

    Returns
    -------
    dict
        return a dictionary containing the structure of the scanned H5 file. No dataset.
    """
    
    node = {
        "path": obj.name,
        "type": "group" if isinstance(obj, h5.Group) else "dataset",
        "attrs": dict(obj.attrs)
    }

    if isinstance(obj, h5.Dataset):
        
        node.update({
            "shape": obj.shape,
            "dtype": obj.dtype.str,
            "chunks": obj.chunks,
            "compression": obj.compression,
            "maxshape": obj.maxshape,
        })
        
    else:
        
        node["children"] = {
            key: __scan_hdf5__(obj[key])
            for key in obj.keys()
        }

    return node


def __find_h5_node_by_path__(template, target_path):
    """
    Return a node from a scanned HDF5 template tree by absolute path.

    Parameters
    ----------
    template : dict
        Root node of the lightweight HDF5 template created by __scan_hdf5__.
    target_path : str
        Absolute HDF5-style path to retrieve (for example: /Acquisition/Raw[0]/RawData).

    Returns
    -------
    dict
        Template node corresponding to ``target_path``.

    Raises
    ------
    KeyError
        If any path component is missing in the template tree.
    """
    if target_path == "/":
        return template

    parts = [p for p in target_path.strip("/").split("/") if p]
    node = template

    for part in parts:
        children = node.get("children", {})
        if part not in children:
            raise KeyError(f'HDF5 template path not found: "{target_path}" (missing "{part}")')
        node = children[part]

    return node


def __update_h5_attr_like__(attrs, key, value):
    """
    Update an HDF5 attribute while preserving the original storage style.

    Parameters
    ----------
    attrs : dict
        Attribute dictionary of a template node.
    key : str
        Attribute name to update.
    value : object
        New value for the attribute.

    """
    if key not in attrs:
        attrs[key] = value
        return

    current = attrs[key]
    if isinstance(current, np.ndarray):
        out = np.array(current, copy=True)
        if out.size == 0:
            attrs[key] = np.array([value], dtype=out.dtype)
        else:
            out.flat[0] = value
            attrs[key] = out
    elif isinstance(current, np.generic):
        attrs[key] = np.array(value, dtype=current.dtype).item()
    else:
        attrs[key] = value


def __clone_template__(obj):
    """
    Clone a template object recursively with safe fallbacks for HDF5 objects.

    Parameters
    ----------
    obj : object
        Template object to clone.

    Returns
    -------
    object
        Cloned object. Non-copyable objects (for example some HDF5 references)
        are returned as-is.

    """
    
    if isinstance(obj, dict):
        return {k: __clone_template__(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [__clone_template__(v) for v in obj]
    if isinstance(obj, tuple):
        return tuple(__clone_template__(v) for v in obj)
    if isinstance(obj, np.ndarray):
        return np.array(obj, copy=True)
    if isinstance(obj, (str, bytes, int, float, bool, type(None), np.generic)):
        return obj

    try:
        return copy.deepcopy(obj)
    except Exception:
        # Keep object as-is if it cannot be deep-copied (happens often with h5py references). These H5 are waaay more complex and more diverse than what I thought.
        return obj
