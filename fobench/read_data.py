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



def read_data(filepath=None, company=None, range_ch=None, format=None):
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

    if format == 'tdms' and company == 'silixa': # Silixa TDMS

        print('Reading TDMS file (Silixa Format)...')
        file_file = tdms.TdmsFile.read(filepath)
        properties = file_file.properties
        dataset = None
        chans = file_file['Measurement'].channels() if not range_ch else [file_file['Measurement'].channels()[i] for i in range_ch]
        chans_nums = [int(chan.name) for chan in chans]
        fiber = properties['name'].split('_')[0]
        sampling_frequency = properties['SamplingFrequency[Hz]']
        dt = 1/sampling_frequency
        start_time = UTC(properties['ISO8601 Timestamp'])
        end_time = UTC(start_time + len(chans[0])*dt)
        spatial_interval = properties['SpatialResolution[m]']
        time_length = end_time - start_time
        num_points = int(time_length/dt)
        gauge_length = properties['GaugeLength']
        channel_offset = properties['OffsetLength']
        units = 'counts'
        conv_factor = None # conversion factor if given explicitly
            
    if (format == 'h5' or format == 'hdf5') and company == 'febus': # FEBUS HDF5
    
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
        sampling_frequency = 1/(properties['Spacing'][1]*1e-3)
        dt = 1/sampling_frequency
        start_time = UTC(file_file[instrument]['Source1']['time'][0]) #Time is retaken to correct for LAG.
        end_time = UTC(file_file[instrument]['Source1']['time'][-1]) + (dt * LAG)
        spatial_interval = properties['Spacing'][0]
        time_length = end_time - start_time
        num_points = int(time_length/dt)
        gauge_length = properties['GaugeLength'] # CHECK THIS!!
        channel_offset = 0 # FIX THIS!!
        units = 'counts'
        conv_factor = None # conversion factor if given explicitly

    if (format == 'h5' or format == 'hdf5') and company == 'silixa': # Silixa HDF5
    
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
        fiber = properties['FibreType']
        sampling_frequency = properties['OutputDataRate']
        dt = 1/sampling_frequency
        start_time = UTC(properties['PartStartTime'])
        end_time = UTC(start_time + properties['Count']*dt)
        spatial_interval = properties['SpatialResolution']
        time_length = end_time - start_time
        num_points = properties['Count']
        gauge_length = properties['GaugeLength']
        channel_offset = abs(properties['PreTriggerSamples'])
        units = 'counts'
        conv_factor = None # conversion factor if given explicitly
        
    if format == 'npy' and company == 'bam': # .npy format for BAM. This might fail always since the unit is NON-COMMERCIAL!
    
        print('File format is a Numpy Class. It contains only the unitsdata, and so the metadata must be filled automatically in the code.')
        file_file = None
        properties = None
        chans = None
        dataset = np.load(filepath)
        chans_nums = [i for i in range(dataset.shape[1])]
        chans = np.array(chans_nums)
        fiber = 'La Chida'
        sampling_frequency = 100000
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

    if format == 'npz' and company == 'bam': # .npy format for BAM. This might fail always since the unit is NON-COMMERCIAL!
    
        print('File format is a Numpy Zip Class. No Gaueg Length specified. Do not attempt to convert to Stran-Rate.')
        file_file = None
        properties = None
        chans = None
        dataset = np.load(filepath)
        chans_nums = [i for i in range(len(dataset['distance']))]
        chans = np.array(chans_nums)
        fiber = 'La Chida'
        sampling_frequency = dataset['freq']
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
    
    #file_file.close()
    
    if (format == 'h5' or format == 'hdf5') and company == 'terra15': # Terra15 HDF5
        
        print('Reading HDF5 file (Terra15 Format)...')
        file_file = h5.File(filepath,'r')
        properties = file_file.attrs
        dataset = file_file['data_product']
        chans_nums = [i for i in range(properties['nx'])]
        chans = np.array(chans_nums)
        fiber = 'standard'
        dt = float(properties['dt_computer'])
        sampling_frequency = 1 / dt
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
        
    if (format == 'h5' or format == 'hdf5') and company == 'asn': # ASN OptoDAS HDF5 (It can be a bit more complex, so I'm trying to make it simple!)
        
        print('Reading HDF5 file (ASN Format)...')
        file_file = h5.File(filepath,'r')
        properties = file_file['acqSpec']
        dataset = file_file['data']
        chans_nums = [i for i in range(int(file_file['header']['dimensionRanges']['dimension1']['size'][()]))]
        chans = np.array(chans_nums)
        fiber = 'standard'
        dt = float(file_file['header']['dt'][()])
        sampling_frequency = 1 / dt
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
        
    if (format == 'h5' or format == 'hdf5') and company == 'quantx': # QuantX OptoaSense HDF5
        
        print('Reading HDF5 file (QuantX Format)...')
        file_file = h5.File(filepath,'r')
        properties = file_file['Acquisition'].attrs
        dataset = file_file['Acquisition']['Raw[0]']
        chans_nums = [i for i in range(int(file_file['Acquisition']['Raw[0]'].attrs['NumberOfLoci']))]
        chans = np.array(chans_nums)
        fiber = 'standard'
        sampling_frequency = float(dataset.attrs['OutputDataRate'])
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
        'dt',
        'start_time',
        'end_time',
        'spatial_interval',
        'num_points',
        'time_length',
        'gauge_length',
        'channel_offset',
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
                dt,
                start_time,
                end_time,
                spatial_interval,
                num_points,
                time_length,
                gauge_length,
                channel_offset,
                units,
                conv_factor
                )

    return attributes


# Recurive method to convert all h5py Objects into dictionaries.
def h5_to_dict(h5_obj):
    
    result = {key: h5_to_dict(item) if isinstance(item, h5.Group) else item[()] if isinstance(item, h5.Dataset) else item for key, item in h5_obj.items()}
    
    return result
