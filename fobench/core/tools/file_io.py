"""File managing the different fiber optic sensing formats.

:Authors:
    - Sergio Diaz-Meza
    - Jonas Pätzel

:Contributors:
    - Christopher Wollin

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

def read_data(filepath: str = None, company: str = None, range_ch: int|list = None,
              format: str = None, load_data: bool = True,
              show_progress: bool = True, storage_opts = None):
    """
    Function to read out the data of various data formats.

    Parameters
    ----------
    filepath : str
        Path to file to read.
    company : str
        Interrogator manufacturer. One of ``"silixa"``, ``"febus"``, ``"aragon"``,
        ``"quantx"``, ``"asn"``, ``"terra15"`` ``bam`` and ``"sintela"``
    range_ch : tuple(int, int)
            Range of channels to load.
    format : str, optional
        Format specifier.
    load_data : bool
        If ``False``, only metadata is loaded.
    show_progress : True
        Show progress bar when loading data.
    storage_opts : TYPE, optional
        DESCRIPTION.

    Returns
    -------
    attributes : dict
        Dictionary holding all data and metadata

    Raises
    ------
    ValueError
        Unrecognized file format.
    """

    # modify range_ch variable
    if isinstance(range_ch, int): # check if it is single value
        range_ch = [range_ch]

    elif isinstance(range_ch, np.ndarray): # check if is array
        range_ch = list(range_ch)

    # option for reading S3 stored files. The check of instance is necessary because the remote file not necessary starts with...
    if isinstance(filepath, str) and filepath.startswith("s3://"):

        return s3_file(filepath, company, range_ch=range_ch, format=format, load_data=load_data,
                        show_progress=show_progress, storage_opts=storage_opts)

    template = None

    if format == "tdms" and company == "silixa": # Silixa TDMS

        pbar = tqdm(total=1, leave=True, desc="Reading Silixa TDMS file", disable=not show_progress)
        file_file = tdms.TdmsFile.read_metadata(filepath) # before it was tdms.TdmsFile.read()
        template = None
        properties = file_file.properties
        dataset = None
        measurement = file_file["Measurement"]
        tdms_chann = measurement.channels() if range_ch is None else [measurement.channels()[i] for i in range_ch]
        chans = np.asarray([int(chan.name) for chan in tdms_chann], dtype=int)

        # loading of the data conditioned. Old way
        # data_range = None if range_ch is None else list(chans_nums)
        # data = __data__(measurement, format, company, data_range) if load_data else None

        # new way
        data = __tdms_biReader__(filepath, file_file, range_ch) if load_data else None

        fiber = properties["name"].split("_")[0]
        sampling_frequency = properties["SamplingFrequency[Hz]"]
        o_sampling_frequency = sampling_frequency
        dt = 1/sampling_frequency
        start_time = UTC(properties["ISO8601 Timestamp"])
        end_time = UTC(start_time + (len(tdms_chann[0])-1)*dt)
        spatial_interval = properties["SpatialResolution[m]"]
        time_length = end_time - start_time
        num_points = len(tdms_chann[0])
        gauge_length = properties["GaugeLength"]
        channel_offset = properties["OffsetLength"]
        units = "counts"
        conv_factor = None # conversion factor if given explicitly

    elif (format == "h5" or format == "hdf5") and company == "silixa": # Silixa HDF5

        pbar = tqdm(total=1, leave=True, desc="Reading Silixa HDF5 file", disable=not show_progress)

        with h5.File(filepath,"r") as file_file:

            template = None
            properties = {}
            instrument = file_file.keys()
            dataset = file_file["Acquisition"]["Raw[0]"]["RawData"]

            # load of all hidden properties
            for vv in list(file_file["Acquisition"].attrs): properties[vv] = file_file["Acquisition"].attrs[vv]  #mandatory!!
            for vv in list(file_file["Acquisition"]["Custom"].attrs): properties[vv] = file_file["Acquisition"]["Custom"].attrs[vv]  #optional
            for vv in list(file_file["Acquisition"]["Custom"]["AdvancedUserSettings"].attrs): properties[vv] = file_file["Acquisition"]["Custom"]["AdvancedUserSettings"].attrs[vv]  #mandatory!!
            for vv in list(file_file["Acquisition"]["Custom"]["SystemInformation"]["Chassis"].attrs): properties[vv] = file_file["Acquisition"]["Custom"]["SystemInformation"]["Chassis"].attrs[vv]  #optional
            for vv in list(file_file["Acquisition"]["Custom"]["SystemInformation"]["Devices0"].attrs): properties[vv] = file_file["Acquisition"]["Custom"]["SystemInformation"]["Devices0"].attrs[vv]  #optional
            for vv in list(file_file["Acquisition"]["Custom"]["SystemInformation"]["Devices1"].attrs): properties[vv] = file_file["Acquisition"]["Custom"]["SystemInformation"]["Devices1"].attrs[vv]  #optional
            for vv in list(file_file["Acquisition"]["Custom"]["SystemInformation"]["Devices2"].attrs): properties[vv] = file_file["Acquisition"]["Custom"]["SystemInformation"]["Devices2"].attrs[vv]  #optional
            for vv in list(file_file["Acquisition"]["Custom"]["SystemInformation"]["GPS"].attrs): properties[vv] = file_file["Acquisition"]["Custom"]["SystemInformation"]["GPS"].attrs[vv]  #mandatory
            for vv in list(file_file["Acquisition"]["Custom"]["SystemInformation"]["OSVersion"].attrs): properties[vv] = file_file["Acquisition"]["Custom"]["SystemInformation"]["OSVersion"].attrs[vv]  #optional/mandatory
            for vv in list(file_file["Acquisition"]["Custom"]["SystemInformation"]["ProcessingUnit"].attrs): properties[vv] = file_file["Acquisition"]["Custom"]["SystemInformation"]["ProcessingUnit"].attrs[vv]  #optional
            for vv in list(file_file["Acquisition"]["Custom"]["SystemSettings"].attrs): properties[vv] = file_file["Acquisition"]["Custom"]["SystemSettings"].attrs[vv]  #mandatory!!
            for vv in list(file_file["Acquisition"]["Custom"]["UserSettings"].attrs): properties[vv] = file_file["Acquisition"]["Custom"]["UserSettings"].attrs[vv]  #mandatory!!
            for vv in list(file_file["Acquisition"]["Raw[0]"].attrs): properties[vv] = file_file["Acquisition"]["Raw[0]"].attrs[vv]  #mandatory!!
            for vv in list(file_file["Acquisition"]["Raw[0]"]["RawData"].attrs): properties[vv] = file_file["Acquisition"]["Raw[0]"]["RawData"].attrs[vv]  #optional
            for vv in list(file_file["Acquisition"]["Raw[0]"]["RawDataTime"].attrs): properties[vv] = file_file["Acquisition"]["Raw[0]"]["RawDataTime"].attrs[vv]  #optional

            chans = np.arange(properties["NumberOfLoci"], dtype=int) if range_ch is None else np.asarray(range_ch, dtype=int)
            # loading the data conditioned
            data_range = None if range_ch is None else chans.tolist()

            # data loading
            # data = None
            # if load_data:
            #     data = __h5_biReader__(filepath, dataset, data_range, copy_data=True)
            #     if data is None:
            #         data = __data__(dataset, format, company, data_range)
            data = __data__(dataset, format, company, data_range) if load_data else None
            fiber = properties["FibreType"]
            sampling_frequency = properties["OutputDataRate"]
            o_sampling_frequency = sampling_frequency
            dt = 1/sampling_frequency
            start_time = UTC(properties["PartStartTime"])
            end_time = UTC(start_time + (properties["Count"]-1)*dt)
            spatial_interval = properties["SpatialResolution"]
            time_length = end_time - start_time
            num_points = properties["Count"]
            gauge_length = properties["GaugeLength"]
            channel_offset = abs(properties["PreTriggerSamples"])
            units = "counts"
            conv_factor = None # conversion factor if given explicitly

    # Required checking! Contact providers! TEST!
    elif (format == "h5" or format == "hdf5") and company == "febus": # FEBUS HDF5

        pbar = tqdm(total=1, leave=True, desc="Reading Febus HDF5 file", disable=not show_progress)

        with h5.File(filepath,"r") as file_file:

            instrument = list(file_file.keys())[0]
            template = None
            properties = dict(file_file[instrument]["Source1"]["Zone1"].attrs)
            fiber = "febus" # VARIABLE TO CHANGE FOR REAL NAME EXTRACTION
            measure_type = list(file_file[instrument]["Source1"]["Zone1"].keys())[0]
            dataset = file_file[instrument]["Source1"]["Zone1"][measure_type]
            chans = np.arange(dataset.shape[2], dtype=int) if range_ch is None else np.asarray(range_ch, dtype=int)
            # loading the data conditioned
            #LAG = properties["BlockOverlap"][0] # 201 is standard. How much the data is repeated in batches.
            sampling_frequency = 1/(properties["Spacing"][1]*1e-3)
            o_sampling_frequency = sampling_frequency
            dt = 1/sampling_frequency
            SEMILAG = properties["BlockRate"][0]*1e-3 if getattr(properties["BlockRate"], "size", 0) == 1 else properties["BlockRate"]*1e-3 # unpacking value sif is inside list.
            SEMILAG = int(np.round((1/SEMILAG) / dt))
            data_range = None if range_ch is None else chans.tolist()

            # data loading
            # data = None
            # if load_data:
            #     data = __h5_biReader__(filepath, dataset, data_range, copy_data=True)
            #     if data is None:
            #         data = __data__(dataset, format, company, data_range)
            data = __data__(dataset, format, company, data_range, SEMILAG) if load_data else None
            start_time = UTC(file_file[instrument]["Source1"]["time"][0]) #Time is retaken to correct for LAG.
            end_time = UTC(file_file[instrument]["Source1"]["time"][-1]) + (dt * SEMILAG)
            spatial_interval = properties["Spacing"][0]
            time_length = end_time - start_time
            num_points = int(time_length/dt) # check!
            gauge_length = properties["GaugeLength"][0] # CHECK THIS!!
            channel_offset = 0 # FIX THIS!!
            units = "counts"
            conv_factor = None # conversion factor if given explicitly

    elif (format == "h5" or format == "hdf5") and company == "terra15": # Terra15 HDF5

        pbar = tqdm(total=1, leave=True, desc="Reading Terra15 HDF5 file", disable=not show_progress)

        with h5.File(filepath,"r") as file_file:

            properties = dict(file_file.attrs)
            template = None
            dataset = file_file["data_product"]
            chans = np.arange(properties["nx"], dtype=int) if range_ch is None else np.asarray(range_ch, dtype=int)
            # loading the data conditioned
            data_range = None if range_ch is None else chans.tolist()

            # data loading
            # data = None
            # if load_data:
            #     data = __h5_biReader__(filepath, dataset["data"], data_range, copy_data=True)
            #     if data is None:
            #         data = __data__(dataset["data"], format, company, data_range)
            data = __data__(dataset["data"], format, company, data_range) if load_data else None
            fiber = "standard"
            dt = float(properties["dt_computer"])
            sampling_frequency = 1 / dt
            o_sampling_frequency = sampling_frequency
            num_points = int(properties["nt"])
            start_time = UTC(properties["file_start_gps_time"]) if properties["file_start_gps_time"] else UTC(properties["file_start_computer_time"])
            end_time = UTC(start_time + num_points * dt)
            spatial_interval = float(properties["dx"])
            time_length = end_time - start_time
            gauge_length = float(properties["gauge_length"])
            channel_offset = int(properties["sensing_range_start"] / spatial_interval)
            units = 'strain-rate' if properties['data_product'] == 'strainrate' else properties['data_product']
            units = "m/s" if "velocity" in units else units

            conv_factor = None # conversion factor if given explicitly

    elif (format == "h5" or format == "hdf5") and company == "asn": # ASN OptoDAS HDF5 (It can be a bit more complex, so I"m trying to make it simple!)

        pbar = tqdm(total=1, leave=True, desc="Reading ASN HDF5 file", disable=not show_progress)

        with h5.File(filepath,"r") as file_file:

            template = None
            properties = h5_to_dict(file_file["acqSpec"])
            dataset = file_file["data"]
            n_channels = int(np.asarray(file_file["header"]["dimensionRanges"]["dimension1"]["size"][()]).squeeze())
            chans = np.arange(n_channels, dtype=int) if range_ch is None else np.asarray(range_ch, dtype=int)
            # loading the data conditioned
            data_range = None if range_ch is None else chans.tolist()

            # data loading
            # data = None
            # if load_data:
            #     data = __h5_biReader__(filepath, dataset, data_range, copy_data=True)
            #     if data is None:
            #         data = __data__(dataset, format, company, data_range)
            data = __data__(dataset, format, company, data_range) if load_data else None
            fiber = "standard"
            dt = float(file_file["header"]["dt"][()])
            sampling_frequency = 1 / dt
            o_sampling_frequency = sampling_frequency
            num_points = int(file_file["header"]["dimensionRanges"]["dimension0"]["size"][()].squeeze())
            start_time = UTC(float(file_file["header"]["time"][()]))
            end_time = UTC(start_time + num_points * dt)
            original_channels = file_file["header"]["channels"]
            spatial_interval = float(file_file["header"]["dx"][()]) * (int(original_channels[1]) - int(original_channels[0]))
            time_length = end_time - start_time
            gauge_length = float(file_file["header"]["gaugeLength"][()])
            channel_offset = int(original_channels[0])
            units_value = file_file["header"]["sensitivityUnits"][()].squeeze().item()
            units = units_value.decode("utf-8") if isinstance(units_value, bytes) else str(units_value)
            conv_factor = float(file_file["header"]["sensitivities"][()].squeeze())
            # units = str(file_file["header"]["sensitivityUnits"][()])[3:-2]
            # conv_factor = file_file["header"]["sensitivities"][0,0]

    elif (format == "h5" or format == "hdf5") and company == "optasense": # OptaSense HDF5

        pbar = tqdm(total=1, leave=True, desc="Reading OptaSense HDF5 file", disable=not show_progress)

        with h5.File(filepath,"r") as file_file:

            properties = dict(file_file["Acquisition"].attrs)
            template = None
            dataset = file_file["Acquisition"]["Raw[0]"]
            chans = np.arange(int(file_file["Acquisition"]["Raw[0]"].attrs["NumberOfLoci"]), dtype=int) if range_ch is None else np.asarray(range_ch, dtype=int)

            # loading the data conditioned
            data_range = None if range_ch is None else chans.tolist()

            # data loading
            # data = None
            # if load_data:
            #     data = __h5_biReader__(filepath, dataset["RawData"], data_range, copy_data=True)
            #     if data is None:
            #         data = __data__(dataset["RawData"], format, company, data_range)
            data = __data__(dataset["RawData"], format, company, data_range) if load_data else None
            fiber = "standard"
            sampling_frequency = float(dataset.attrs["OutputDataRate"])
            o_sampling_frequency = sampling_frequency
            dt = 1 / sampling_frequency
            num_points = int(file_file["Acquisition"]["Raw[0]"]["RawDataTime"].attrs["Count"])
            start_time = UTC(str(properties["MeasurementStartTime"])[2:-1])
            end_time = UTC(start_time + num_points * dt)
            spatial_interval = float(properties["SpatialSamplingInterval"])
            time_length = end_time - start_time
            gauge_length = float(properties["GaugeLength"])
            channel_offset = int(properties["StartLocusIndex"])
            units = str(dataset.attrs["RawDataUnit"])[2:-1]
            conv_factor = None # conversion factor if given explicitly

    elif (format == "h5" or format == "hdf5") and company == "aragon": # Aragon Photonics HDAS HDF5

        pbar = tqdm(total=1, leave=True, desc="Reading Aragon Photonics HDF5 file", disable=not show_progress)

        with h5.File(filepath,"r") as file_file:

            template = None
            properties = dict(file_file.attrs)
            dataset = file_file["strain"]
            if range_ch is None:
                chans = np.arange(file_file["position"].size, dtype=int)
                data_range = None
            else:
                chans = np.arange(range_ch, dtype=int)
                data_range = chans.tolist()
            # loading the data conditioned
            # data loading
            # data = None
            # if load_data:
            #     data = __h5_biReader__(filepath, dataset, data_range, copy_data=True)
            #     if data is None:
            #         data = __data__(dataset, format, company, data_range)
            #     else:
            #         data *= 1e-9
            data = __data__(dataset, format, company, data_range) if load_data else None
            pbar.set_description("Extracting Attributes")
            fiber = "standard"
            sampling_frequency = float(properties["trigger_frequency"][0])
            o_sampling_frequency = sampling_frequency
            dt = 1 / sampling_frequency
            num_points = int(file_file["time"].size)
            start_time = UTC(properties["Global_FileSaveIO_TimeStamp_s"][0])
            end_time = UTC(start_time + num_points * dt)
            spatial_interval = float(properties["spatial_sampling"][0])
            time_length = end_time - start_time
            gauge_length = float(properties["Global_RAM_User_SET_Pulse_Width_(meter)"][0])
            channel_offset = int(properties["fiber_position_offset"][0]/spatial_interval)
            units = list(file_file.keys())[2]
            conv_factor = None # conversion factor if given explicitly

    elif (format == "h5" or format == "hdf5") and company == "sintela": # Sintela Onyx HDF5
        print("Data units (strain, strain-rate...) can not be extracted from Sintela files\n"
             "you can set it manually by editing Fiber.units")
        pbar = tqdm(total=1, leave=True, desc="Reading Sintela HDF5 file", disable=not show_progress)

        with h5.File(filepath, "r") as file_file:

            template = None
            properties = dict(file_file["Acquisition"].attrs)
            dataset = file_file["/Acquisition/Raw[0]/RawData"]
            if range_ch is None:
                chans = np.arange(properties["NumberOfLoci"], dtype=int)
                data_range = None
            else:
                chans = np.asarray(range_ch, dtype=int)
                data_range = chans.tolist()

            # data loading
            # data = None
            # if load_data:
            #     data = __h5_biReader__(filepath, dataset, data_range, copy_data=True)
            #     if data is None:
            #         data = __data__(dataset, format, company, data_range)
            data = __data__(dataset, format, company, data_range) if load_data else None
            fiber = properties["FiberID"]
            sampling_frequency = properties["PulseRate"]
            o_sampling_frequency = sampling_frequency
            num_points = dataset.shape[0]
            start_time = UTC(properties["MeasurementStartTime"].decode("utf-8")) if isinstance(properties["MeasurementStartTime"],bytes) else UTC(properties["MeasurementStartTime"])
            dt = 1 / sampling_frequency
            end_time = UTC(start_time + num_points * dt)
            time_length = end_time - start_time
            spatial_interval = float(properties["SpatialSamplingInterval"])
            channel_offset = properties["StartLocusIndex"]*spatial_interval
            gauge_length = float(properties["GaugeLength"])
            units = "units" # unit is not coded in metadata for Sintela .h5
            conv_factor = None

	# ####################################################
	# CAUTION!! NON OFFICIAL / EXPERIMENTAL FORMATS, ONLY FOR SPECIAL CASES.
	# ####################################################

    elif format == "npy" and company == "bam": # .npy format for BAM. This might fail always since the unit is NON-COMMERCIAL!

        print("File format is a Numpy Class. It contains only the unitsdata, and so the metadata must be filled automatically in the code.")
        file_file = None
        properties = {}
        chans = None
        dataset = np.load(filepath)
        chans = np.arange(dataset.shape[1], dtype=int) if range_ch is None else np.asarray(range_ch, dtype=int)
        # loading the data conditioned
        data_range = None if range_ch is None else chans.tolist()
        data = __data__(filepath, format, company, data_range) if load_data else None
        fiber = "La Chida"
        sampling_frequency = 100000
        o_sampling_frequency = sampling_frequency
        dt = 1/sampling_frequency
        start_time = UTC("2023-03-01T00:00:00")
        spatial_interval = 0.4
        num_points = int(dataset.shape[0])
        end_time = UTC(start_time + num_points*dt)
        time_length = end_time - start_time
        gauge_length = 2
        channel_offset = None
        units = None
        conv_factor = None # conversion factor if given explicitly

    elif format == "npz" and company == "bam": # .npy format for BAM. This might fail always since the unit is NON-COMMERCIAL!

        print("File format is a Numpy Zip Class. No Gauge Length specified. Do not attempt to convert to Strain-Rate.")
        file_file = None
        properties = {}
        chans = None
        dataset = np.load(filepath)
        chans = np.arange(len(dataset["distance"]), dtype=int) if range_ch is None else np.asarray(range_ch, dtype=int)
        # loading the data conditioned
        data_range = None if range_ch is None else chans.tolist()
        data = __data__(filepath, format, company, data_range) if load_data else None
        fiber = "La Chida"
        sampling_frequency = dataset["freq"]
        o_sampling_frequency = sampling_frequency
        dt = 1/sampling_frequency
        start_time = UTC("2023-03-01T00:00:00")
        spatial_interval = dataset["distance"][1] - dataset["distance"][0]
        num_points = int(dataset["time"][:-2].shape[0])
        end_time = UTC(start_time + num_points*dt)
        time_length = end_time - start_time
        gauge_length = 2
        channel_offset = None
        units = None
        conv_factor = None # conversion factor if given explicitly

    elif (format == "h5" or format == "hdf5") and (company == "michele"): # .h5 format for Michele INGV"s decimated files from Silixa.

        pbar = tqdm(total=1, leave=True, desc="Reading Michele INGV decimated HDF5 file", disable=not show_progress)

        with h5.File(filepath, "r") as file_file:

            properties = dict(file_file.attrs)
            template = None
            dataset = file_file["Fiber"]
            chans = np.asarray(file_file["ChannelMap"]) if range_ch is None else np.asarray(range_ch, dtype=int)
            # loading the data conditioned
            data_range = None if range_ch is None else chans.tolist()

            # data loading
            # data = None
            # if load_data:
            #     data = __h5_biReader__(filepath, dataset, data_range, copy_data=True)
            #     if data is None:
            #         data = __data__(dataset, format, company, data_range)
            data = __data__(dataset, format, company, data_range) if load_data else None
            fiber = properties["Fibre Type"].decode("UTF-8")
            dt = float(properties["Sampletime"][0])
            sampling_frequency = 1 / dt
            o_sampling_frequency = properties["SamplingFrequency[Hz]"][0]
            num_points = int(dataset.shape[0])
            start_time = UTC(properties["StartTime_txt"].decode("UTF-8")) if properties["StartTime_txt"] else UTC(properties["CPUTimeStamp_txt"].decode("UTF-8"))
            end_time = UTC(start_time + (num_points-1) * dt)
            spatial_interval = float(properties["SpatialResolution[m]"][0])
            time_length = end_time - start_time
            gauge_length = float(properties["GaugeLength"][0])
            channel_offset = int(properties["OffsetLength"])
            units = "counts"
            conv_factor = None # conversion factor if given explicitly

    else:
        raise ValueError(f"{format} is not a recognized file format!")
    del dataset

    # Attributed for the Fiber class.
    attr_keys = [
        "basefile",
        "format",
        "company",
        "fiber",
        "properties",
        "chans",
        "total_channels",
        "sampling_frequency",
        "o_sampling_frequency",
        "dt",
        "start_time",
        "end_time",
        "spatial_interval",
        "num_points",
        "time_length",
        "gauge_length",
        "channel_offset",
        "data",
        "units",
        "conv_factor"
        ]

    # coonvert the type of the data to floating, ready for processing
    # if not np.issubdtype(data.dtype, np.floating) and load_data is True:
    #     data = data.astype(float, copy=False)

    attributes = [template if template is not None else filepath,
                format,
                company,
                fiber,
                h5_to_dict(properties),
                np.array(chans),
                len(chans),
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
    pbar.set_description("Read File ✓")
    # pbar.close()
    return attributes



"""Core reading functions"""

def h5_to_dict(h5_obj):
    """Recursive method to convert all h5py Objects into dictionaries."""
    result = {key: h5_to_dict(item) if isinstance(item, h5.Group) else item[()] if isinstance(item, h5.Dataset) else item for key, item in h5_obj.items()}
    return result


def __data__(extract_point, format, company, range_ch, LAG=None):
    """Extracts the data depending on the file type"""

    # it would be nice if manufacturer would save their data as (distance, time) shape initially.
    # it would allow for fast processing.
    if format in ("h5", "hdf5"):

        if company in ("silixa", "asn", "optasense", "aragon", "sintela", "terra15", "michele"):
            if range_ch is None:
                values = extract_point[:,:]
            else:
                values = extract_point[:, range_ch]
            if company == "aragon":
                values *= 1e-9  # Convert from nanostrain to strain
            # return np.ascontiguousarray(values.T)
            return values

        elif company == "febus":
            dims = extract_point.shape
            values = extract_point[:, :LAG, :].reshape(dims[0] * LAG, dims[2])
            if range_ch is not None:
                values = values[:, range_ch]
            # return np.ascontiguousarray(values.T)
            return values

    elif format == "tdms" and company == "silixa":

        values = extract_point.as_dataframe().to_numpy()
        if range_ch is not None:
            values = values[:, range_ch]
        # return np.ascontiguousarray(values.T)
        return values

    # ####################################################
    # CAUTION!! NON OFFICIAL / EXPERIMENTAL FORMATS, ONLY FOR SPECIAL CASES.
    # ####################################################

    elif format == "npy" and company == "bam":
        values = np.load(extract_point)
        if range_ch is not None:
            values = values[:, range_ch]
        # return np.ascontiguousarray(values.T)
        return values

    elif format == "npz" and company == "bam":
        values = np.load(extract_point)["ph"]
        if range_ch is not None:
            values = values[:, range_ch]
        # return np.ascontiguousarray(values.T)
        return values


def s3_file(filepath, company, range_ch=None, format=None, load_data=True, show_progress=True, storage_opts=None):
    """s3_file Reads a fibre optic sensing file stored in a S3-compatible filesystem.

    Parameters
    ----------
    filepath : str
        Complete S3 URI of the file
    company : str
        Manufacturer of the fibre optic sensing unit that produce the data (check read_data())
    range_ch : int or list, optional
        Channels range to load
    format : str, optional
        Format of the file to read (suffix)
    load_data : bool, optional
        Load the data when variable is ``True``. else it's only the metadata
    show_progress : bool, optional
        Shows the file-reading progress bar
    storage_opts : dict, optional
        Options to be passed to the S3 filesystem. Needed as credentials and temporal file handling, by default None

    Returns
    -------
    dict
        Fiber data and metadata attributes

    """


    import s3fs
    s3 = s3fs.S3FileSystem(**(storage_opts or {}))

    with s3.open(filepath,"rb") as remote_file:

        attributes = read_data(filepath=remote_file, company=company, range_ch=range_ch,
                                format=format, load_data=load_data, show_progress=show_progress)
    attributes["basefile"] = filepath

    return attributes


def __tdms_biReader__(filepath, tdms_metadata, range_ch, copy_data=True):
    """Fast simple reader of TDMS as binary, superior to nptdms package,
    but only works for a certain specific layout. Lets hope that tdms stays like that.
    DAMN I HATE TDMS, so unceessary, so inconvenient.

    Parameters
    ----------
    filepath : str
        global or partial path where the data is located.
    tdms_metadata : _type_
        _description_
    range_ch : array or list
        range fo channels to slice the data
    copy_data : bool, optional
        load the data into memory, by default True

    Returns
    -------
    np.array()
        full data matrix of the sensing

    """

    segments = tdms_metadata._reader._segments
    segment = segments[0]

    content = tdms_metadata["Measurement"]
    channels = content.channels()
    n_channels, n_samples, dtype = len(channels), len(channels[0]), np.dtype(channels[0].dtype)
    n_bytes = segment.next_segment_pos -segment.data_position

    values = np.memmap(filepath, mode="r", dtype=dtype, offset=segment.data_position, shape=(n_samples, n_channels), order="C")

    if range_ch is not None:
        values = values[:,range_ch]
    if copy_data:
        # values = np.array(values.T, copy=True, order="C")
        values = np.array(values, copy=True)

    return values


def __h5_biReader__(filepath, dataset, range_ch=None, copy_data=True):
    """Fast simple reader of H5 as binary, superior to h5py package,
    but only works for a certain specific layout.

    Parameters
    ----------
    filepath : str
        global or partial path where the data is located.
    tdms_metadata : _type_
        _description_
    range_ch : array or list
        range fo channels to slice the data
    copy_data : bool, optional
        load the data into memory, by default True

    Returns
    -------
    np.array()
        full data matrix of the sensing

    """

    # if the dataset has the following complications, then the fast reader would not work, and so we must fall back to the old reader.
    # still the old one is fast, but this one is muuuuch faster :)
    if not isinstance(dataset, h5.Dataset):
        return None

    if not isinstance(filepath, str):
        return None

    if dataset.chunks is not None:
        return None

    if dataset.compression is not None:
        return None

    if dataset.shuffle or dataset.fletcher32:
        return None

    offset = dataset.id.get_offset()

    if offset is None:
        return None

    values = np.memmap(filepath, mode="r", dtype=dataset.dtype, offset=offset, shape=dataset.shape, order="C")

    if range_ch is not None:
        values = values[:, range_ch]

    if copy_data:
        values = np.array(values, copy=True)

    return values


"""Methods for writing in original format"""

def write_data(Fiber, filepath=None, company=None):
    """Manages saving of data in the original format as it was written from the manufacturer.
    Check format and company available (format, company) from read_data() function.
    "We can possibly try to give options to save in another format but it will
    require to known previosuly the formats somehow"

    Parameters
    ----------
    Fiber : Fiber class
        Fiber class instance.
    filepath : str, optional
        Complete path including name to be saved and the format after the dot.
        "path/where/to/save/file.h5".
    company : str, optional
        Manufacturer or the instrument that generates the data.

    Raises
    ------
    ValueError
        Filepath is required.

    """

    if filepath is None:
        raise ValueError("filepath is required")

    format = filepath.split(".")[-1] # accessing the desired format from the string of the file path.

    if format == "tdms" and company == "silixa": # Silixa TDMS

        pbar = tqdm(total=1, leave=True, desc="Saving in Silixa TDMS file")
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
        selected = [ch_map[str(c)] for c in Fiber.channels]

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

        pbar = tqdm(total=1, leave=True, desc="Saving in Sintela HDF5 file")
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
            __h5_writer__(template, f)

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
        if Fiber.channels is not None and len(Fiber.channels) == Fiber.data.shape[1]:
            p_vals = np.asarray(Fiber.channels, dtype=np.float64) * float(Fiber.spatial_interval)
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
            __h5_writer__(template, f)

    if (format == 'h5' or format == 'hdf5') and company == 'asn': # ASN H5 files.

        pbar = tqdm(total=1, leave=True, desc='Saving in ASN HDF5 file')
        template = __clone_template__(Fiber.__basefile__)

        # Core data and mandatory header paths.
        data_node = __find_h5_node_by_path__(template, "/data")
        dt_node = __find_h5_node_by_path__(template, "/header/dt")
        t0_node = __find_h5_node_by_path__(template, "/header/time")
        n0_node = __find_h5_node_by_path__(template, "/header/dimensionRanges/dimension0/size")
        n1_node = __find_h5_node_by_path__(template, "/header/dimensionRanges/dimension1/size")
        ch_node = __find_h5_node_by_path__(template, "/header/channels")
        dx_node = __find_h5_node_by_path__(template, "/header/dx")
        gl_node = __find_h5_node_by_path__(template, "/header/gaugeLength")

        # Preserve on-disk dtype as much as possible for compatibility/size.
        out = np.asarray(Fiber.data)
        out_dtype = np.dtype(data_node["dtype"])
        if np.issubdtype(out_dtype, np.integer):
            info = np.iinfo(out_dtype)
            out = np.rint(np.clip(out, info.min, info.max)).astype(out_dtype)
        else:
            out = out.astype(out_dtype, copy=False)
        data_node["data"] = out
        data_node["shape"] = out.shape
        chunks = data_node.get("chunks")
        if chunks is not None:
            if len(chunks) != len(out.shape) or any(ch > sh for ch, sh in zip(chunks, out.shape)):
                data_node["chunks"] = None

        # Scalar datasets in header.
        dt_val = np.array(1.0 / float(Fiber.sampling_frequency), dtype=np.dtype(dt_node["dtype"]))
        dt_node["data"] = dt_val
        dt_node["shape"] = dt_val.shape

        t0_val = np.array(float(Fiber.start_time), dtype=np.dtype(t0_node["dtype"]))
        t0_node["data"] = t0_val
        t0_node["shape"] = t0_val.shape

        n0_val = np.array(int(Fiber.data.shape[0]), dtype=np.dtype(n0_node["dtype"]))
        n1_val = np.array(int(Fiber.data.shape[1]), dtype=np.dtype(n1_node["dtype"]))
        n0_node["data"], n0_node["shape"] = n0_val, n0_val.shape
        n1_node["data"], n1_node["shape"] = n1_val, n1_val.shape

        # Channels and spacing.
        if Fiber.channels is not None and len(Fiber.channels) == Fiber.data.shape[1]:
            channels = np.asarray(Fiber.channels, dtype=np.float64)
        else:
            channels = np.arange(Fiber.data.shape[1], dtype=np.float64) + float(Fiber.channel_offset)
        ch_dtype = np.dtype(ch_node["dtype"])
        if np.issubdtype(ch_dtype, np.integer):
            info = np.iinfo(ch_dtype)
            channels = np.rint(np.clip(channels, info.min, info.max)).astype(ch_dtype)
        else:
            channels = channels.astype(ch_dtype, copy=False)
        ch_node["data"] = channels
        ch_node["shape"] = channels.shape

        chan_step = 1.0
        if channels.size > 1:
            step = float(channels[1]) - float(channels[0])
            if step != 0:
                chan_step = step
        dx_val = np.array(float(Fiber.spatial_interval) / chan_step, dtype=np.dtype(dx_node["dtype"]))
        dx_node["data"] = dx_val
        dx_node["shape"] = dx_val.shape

        gl_val = np.array(float(Fiber.gauge_length), dtype=np.dtype(gl_node["dtype"]))
        gl_node["data"] = gl_val
        gl_node["shape"] = gl_val.shape

        # Optional unit update if present in file structure.
        try:
            u_node = __find_h5_node_by_path__(template, "/header/sensitivityUnits")
            u_dtype = np.dtype(u_node["dtype"])
            if u_dtype.kind in ("S", "a"):
                u_val = np.array(str(Fiber.units).encode("utf-8"), dtype=u_dtype)
            elif u_dtype.kind == "U":
                u_val = np.array(str(Fiber.units), dtype=u_dtype)
            else:
                u_val = np.array(str(Fiber.units))
            u_node["data"] = u_val
            u_node["shape"] = u_val.shape
        except KeyError:
            pass

        with h5.File(filepath, "w") as f:
            __h5_writer__(template, f)

    pbar.update(1)
    pbar.set_description("File Saved ✓")

"""H5/HDF5 explorers, templates, internals and writers"""

def __h5_writer__(template, file):
    """Handles writing of H5/HDF5 files.

    Parameters
    ----------
    template : dict
        Lightweight template of the h5 structure file.
    file : h5 obj
        Opened H5 file to where to save the data.

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
            dtype_out, data = __coerce_h5_dtype_and_data__(template.get("dtype"), data)
            shape_out = template["shape"] if data is None else np.asarray(data).shape
            kwargs = {
                "dtype": dtype_out,
                "chunks": template.get("chunks"),
                "compression": template.get("compression"),
                "maxshape": template.get("maxshape"),
            }
            kwargs = {k: v for k, v in kwargs.items() if v is not None}

            # h5 scalar datasets (shape=() or with empty shapes) cannot be extendable/chunked/compressed.
            if shape_out == ():
                kwargs.pop("maxshape", None)
                kwargs.pop("chunks", None)
                kwargs.pop("compression", None)

            if data is None:
                current_obj = file.create_dataset(
                    name,
                    shape=shape_out,
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
        __h5_writer__(child_node, current_obj)


def __coerce_h5_dtype_and_data__(dtype_in, data):
    """Converts template dtype/data to a HDF5-compatible representation.

    Parameters
    ----------
    dtype_in : str or numpy.dtype or None
        dtype extracted from template metadata.
    data : object
        Dataset payload. Can be ``None`` for placeholder datasets.

    Returns
    -------
    tuple
        ``(dtype_out, data_out)``.

    """

    try:
        np_dtype = np.dtype(dtype_in) if dtype_in is not None else None
    except Exception:
        np_dtype = None

    # Most dtypes are already valid. Unicode/object need normalization below.
    if np_dtype is not None and np_dtype.kind not in ("O", "U", "S"):
        return np_dtype, data

    # Object dtypes are not directly writable by h5py. Map to safe alternatives.
    if data is None:
        return h5.string_dtype(encoding="utf-8"), None

    arr = np.asarray(data)
    if arr.dtype.kind in ("U", "S"):
        def _to_text(x):
            if isinstance(x, (bytes, np.bytes_)):
                return x.decode("utf-8", errors="replace")
            return str(x)
        if arr.shape == ():
            return h5.string_dtype(encoding="utf-8"), _to_text(arr.item())
        out = np.vectorize(_to_text, otypes=[object])(arr)
        return h5.string_dtype(encoding="utf-8"), out

    if arr.dtype.kind != "O":
        return arr.dtype, arr

    # Object payload: infer from first non-None value.
    sample = None
    for vv in arr.ravel():
        if vv is not None:
            sample = vv
            break

    if isinstance(sample, (str, bytes, np.str_, np.bytes_)) or sample is None:
        def _to_text(x):
            if x is None:
                return ""
            if isinstance(x, (bytes, np.bytes_)):
                return x.decode("utf-8", errors="replace")
            return str(x)
        if arr.shape == ():
            return h5.string_dtype(encoding="utf-8"), _to_text(arr.item())
        out = np.vectorize(_to_text, otypes=[object])(arr)
        return h5.string_dtype(encoding="utf-8"), out

    # Fallback: try numeric cast, else stringify.
    try:
        out = arr.astype(np.float64)
        return out.dtype, out
    except Exception:
        if arr.shape == ():
            return h5.string_dtype(encoding="utf-8"), "" if arr.item() is None else str(arr.item())
        out = np.vectorize(lambda x: "" if x is None else str(x), otypes=[object])(arr)
        return h5.string_dtype(encoding="utf-8"), out


def scan_template(filepath, company="", format="", storage_opts=None):

    if isinstance(filepath, str) and filepath.startswith("s3://"):
        import s3fs
        s3 = s3fs.S3FileSystem(**(storage_opts or {}))
        opener = s3.open(filepath, "rb")
    else:
        opener = open(filepath, "rb")

    with opener as raw_file:

        if format == "tdms":
            opened_file = tdms.TdmsFile.read(raw_file)

            return __scan_template__(opened_file, company=company, format=format)

        if format in ("hdf5","h5"):
            with h5.File(raw_file, "r") as opened_file:

                return __scan_template__(opened_file, company=company, format=format)


def __scan_template__(core_file, company="", format=""):
    """Scans the structure of the file and saves it as a lightweight version dictionary.

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
        Lightweight version of the architecture of the file.
        Used for later recreation of the file to save information in original format.

    """

    native_template = None

    if format == "tdms" and company == "silixa": # Silixa TDMS

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

def __scan_hdf5__(obj):
    """Scans a H5/HDF5 paths and structure to a lightweight version.
    Scan recurrence included for deeper layers.

    Parameters
    ----------
    obj : H5 object
        H5/HDF5 object, usually given after reading it.

    Returns
    -------
    dict
        Dictionary containing the structure of the scanned H5 file. No dataset.

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
    """Returns a node from a scanned HDF5 template tree by absolute path.

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
    """Update an HDF5 attribute while preserving the original storage style.

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
    """Clones a template object recursively with safe fallbacks for HDF5 objects.

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