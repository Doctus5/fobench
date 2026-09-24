"""Functions for managing, sorting and manipulating files and ``Databases`` in
folders.

:Authors:
    - Sergio Diaz-Meza
    - Jonas Pätzel

:Contributors:
    - Christopher Wollin

"""

import os
import glob
import warnings
import json
import pandas as pd
import numpy as np

# Import Fiber for both install layouts:
# - editable install (`from fobench.fiber import Fiber`)
# - repo namespace usage (`from fobench.fobench.fiber import Fiber`)
try:
    from fobench.core.fiber import Fiber
except ModuleNotFoundError:
    from ..fobench.core.fiber import Fiber


def scan_folder(folder_path, format=None, storage_opts=None, skip_files=None):
    """Function that scans all within a folder to find all files inside.

    Parameters
    ----------
    folder_path : str
        Global folder path where files will be searched.
    format : str, optional
        If a format extension is specified (without the point), only the files with such an extension will be returned.
    storage_opts : dict, optional
        Options passed to the S3 filesystems (profile, endpoint, anonymous access settings, etc.).
    skip_files : list or str, optional
        File paths or file names to not take into account.

    Returns
    -------
    list of str
        Sorted local paths or S3 UIRs.

    """

    # Old method of scanning. No recursivity though.
    # format = '*' if format == None else '*.' + format
    # files = glob.glob(folder_path + format)

    format = "." + format if not format.startswith(".") else format
    storage_opts = {} if storage_opts is None else dict(storage_opts)

    if str(folder_path).startswith("s3://"): # s3 storage fyle systems
        try:
            import s3fs
        except ImportError:
            raise ImportError("S3 scanning library needed (s3fs).")

        s3 = s3fs.S3FileSystem(**storage_opts) # pass the credentials if you have them localy
        prefix = folder_path.removeprefix("s3://").rstrip("/")
        files = s3.glob(f"{prefix}/**/*{format}")
        files = ["s3://"+file for file in files]

    else: # for the local files
        # Use glob to recursively match files ending with the given extension
        search_pattern = os.path.join(folder_path, '**', f'*{format}')
        files = glob.glob(search_pattern, recursive=True)
        
    skip_files = {str(file) for file in (skip_files or [])}
    files = [file for file in files if file not in skip_files and os.path.basename(file) not in skip_files]

    if not files:
        warnings.warn(f"⚠️ No files with the specified format {format} or files at all were found.",
                      category=UserWarning)

    return files


def files2database(files, company, storage_opts=None):
    """Function reada all the FOS files with Fiber Class and reports a big
    Database from them with all its important metadata.

    Parameters
    ----------
    files : str
        List of complete paths from the files to be handled.
    companry : str
        Name of manufacturer.

    Returns
    -------
    database : pandas.Dataframe
        Dataframe indicating paths and essential metadata from each file.

	"""

    N = len(files) # number of data files.
    filtered_keys = ["start_time", "end_time", "dt", "sampling_rate", "n_channels",
                     "spatial_interval", "gauge_length", "channel_offset", "units", "scale_factor"] # attributes of interest for holding consistency.
    #filtered_keys = ["start_time", "end_time", "dt"] # attributes of interest for holding consistency.

    database = []

    # scan the files
    for i in range(N):

        file = files[i] # select file
        # by checking the start and end times in the files
        print("File: " + file)
        d_file = Fiber(file, company=company, load_data=False, storage_opts=storage_opts)

        info = d_file.metadata(meta_dict=True) # getting public relevant attributes.
        info["scale_factor"] = d_file.instr_correct(return_factor=True)
        info = {key: info[key] for key in filtered_keys if key in info}  # Filter
        info["file"] = d_file.__filepath__ # we also take the filepath.

        # correct format for start and end times into strings (UTC also can do the thing).
        info["start_time"] = pd.to_datetime(info["start_time"].isoformat())
        info["end_time"] = pd.to_datetime(info["end_time"].isoformat())

        database.append(info)
        del d_file # free memory.

    # handling it as a DataFrame
    database = pd.DataFrame.from_dict(database)

    return database


def init_dataframe(dictionary):
    """Initializes a DataFrame from dictionary and converts the dates in
    isoformat to datetime for management.

    Parameters
    ----------
    dictionary : dict
        Database from metadata as dictionary.
    Returns
    -------
    database : pandas.Datframe
        DataFrame version of the database.

    """

    database = pd.DataFrame(dictionary)
    database = metadates_2_isoformat(database, reverse=True)

    return database

def open_metadatafile(json_file):
    """Opens a metadata JSON file.

    Parameters
    ----------
    json_file : str
        Dataframe indicating paths and essential metadata from each file.

    Returns
    -------
    metadata : dict
        Metadata as dictionary.

    """

    with open(json_file, "r") as file:
        metadata = json.load(file)

    return metadata


def chrono_order(database):
    """Aranges all fiber optic sensing files in a chronological way.
    Requires reading the file without loading the data.

    Parameters
    ----------
    database : pandas.Dataframe
        Dataframe indicating paths and essential metadata from each file.

    Returns
    -------
    database : pandas.Dataframe
        Dataframe indicating paths and essential metadata from each file.

	"""

    database = database.sort_values(by=["start_time", "end_time"])

    return database


def database_discontinuities(df, split=True):
    """Checks if there are time and parameters discontinuities along a reported Dataset.
    Necessary to separate chuncks of data. Optionally returns the database split.

    Parameters
    ----------
    df : pandas.Dataframe
        Database indicating paths and essential metadata from each file.

    Returns
    -------
    database_chunks : list of pandas.DataFrame
        Database divided into continuous acquisitions when ``split=True``.
    continuity : pandas.Series
        Boolean continuity indicator when ``split=False``. A ``False`` value
        indicates the beginning of a new acquisition.
	"""

    # attributes of metadata to evaluate continuity.
    var_conditions = ["sampling_rate", "n_channels", "spatial_interval", "gauge_length", "channel_offset", "units", "scale_factor"]
    
    # if database is empty
    if df.empty:
        if split:
            return []
        return pd.Series(dtype=bool, index=df.index)

    # Initialize continuity check by checking time.
    # End time of file i + dt must match the start time of file i+1
    df["start_time"] = pd.to_datetime(df["start_time"])
    df["end_time"] = pd.to_datetime(df["end_time"])

    # check chrono continuity
    expected_start = df["end_time"] + pd.to_timedelta(df["dt"], unit="s")
    continuity = expected_start.shift(1) == df["start_time"]
    continuity.iloc[0] = True
    # continuity.iloc[-1] = True

    # # check discontinuities in others properties.
    for variable in var_conditions:
        
        previous = df[variable].shift(1)
        same_val = (df[variable] == previous) | (df[variable].isna() & previous.isna())
        continuity &= same_val

        # continuity *= df[variable].shift(1, fill_value=df[variable].iloc[0]) == df[variable]

    continuity.iloc[0] = True #first file always starets a Dataset
    # Evaluate if the user wants the databse to be splitted in discontinuous parameters (see selected parameters).
    if not split:
        return continuity
    
    # every false value marks the beginning of a Dataset
    break_indices = [i for i in range(1,len(df)) if not continuity.iloc[i]]
    limits = [0] + break_indices + [len(df)]
    database_chunks = [df.iloc[start:end].copy() for start, end in zip(limits[:-1], limits[1:])]

    # else:
    #     database_chunks = []
    #     last_incident = 0

    #     for i in range(continuity.size):
    #         # Evaluate when there is a chunk.
    #         if (continuity.iloc[i] == False) or (i == continuity.size-1):

    #             database_chunks.append( df.iloc[last_incident:i] ) # we add the detected chunk.
    #             last_incident = i

    return database_chunks

def metadates_2_isoformat(dataframe, reverse=False):
    """Transforms the dates of the ``Dataset`` files metadata into isoformat.
        Convenient for saving since Timestamp objects can not be saved in JSONs.

    Parameters
    ----------
    dataframe : pandas.Dataframe
        DataFrame object to find times and manipulate.
    reverse : bool
        If ``True``, transforms from isoformat to Timestamp object.

    Returns
    -------
    dataframe : pandas.Dataframe
        Modified Dataframe object.

    """

    if reverse == False:
        dataframe["start_time"] = dataframe["start_time"].apply(lambda x: x.isoformat())
        dataframe["end_time"] = dataframe["end_time"].apply(lambda x: x.isoformat())

    else:
        # apply ISO8601 for accepting the fractional precision produced during calling .isoformat() function in obspy UTC
        dataframe["start_time"] = pd.to_datetime(dataframe["start_time"], format="ISO8601")
        dataframe["end_time"] = pd.to_datetime(dataframe["end_time"], format="ISO8601")

    return dataframe

def convert_types(value):
    """Convert NumPy values into JSON-compatible Python values.

    Parameters
    ----------
    value : object
        Value or nested metadata structure to convert.

    Returns
    -------
    object
        JSON-compatible value or metadata structure.
    """

    if isinstance(value, dict):
        return {key: convert_types(item) for key, item in value.items()}

    if isinstance(value, (list, tuple)):
        return [convert_types(item) for item in value]

    if isinstance(value, np.ndarray):
        return convert_types(value.tolist())

    if isinstance(value, np.integer):
        return int(value)

    if isinstance(value, np.floating):
        return float(value)

    if isinstance(value, np.bool_):
        return bool(value)

    return value


def df_time_filtering(df: pd.DataFrame, range: tuple, start_label: str = "start_time",
    end_label: str = "end_time", include_overlaps: bool = False) -> pd.DataFrame:
    """Filters a DataFrame by a requested time range.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to be filtered.
    range : tuple
        Start and end times in ISO format (or any datetime-parsable values).
        Example: ("2019-12-12T00:00:00Z", "2020-12-12T00:00:00Z")
    start_label : str, optional
        Column name containing row start times.
    end_label : str, optional
        Column name containing row end times.
    include_overlaps : bool, optional
        If False (default), keep rows fully contained in ``range``.
        If True, keep rows that overlap ``range``.

    Returns
    -------
    df : pd.DataFrame
        Filtered DataFrame.
    """

    if len(range) != 2:
        raise ValueError("⚠️ 'range' must contain exactly two values: (start_time, end_time).")

    missing_cols = [col for col in (start_label, end_label) if col not in df.columns]

    if missing_cols:
        raise KeyError(f"Missing required column(s): {missing_cols}")

    start_bound = pd.to_datetime(range[0], utc=True, errors="raise")
    end_bound = pd.to_datetime(range[1], utc=True, errors="raise")

    if start_bound > end_bound:
        raise ValueError("`range[0]` must be earlier than or equal to `range[1]`.")

    start_times = pd.to_datetime(df[start_label], utc=True, errors="coerce")
    end_times = pd.to_datetime(df[end_label], utc=True, errors="coerce")

    if include_overlaps:
        # Keep any row interval that intersects the requested range.
        condition = (start_times <= end_bound) & (end_times >= start_bound)
    else:
        # Keep only rows fully contained in the requested range.
        condition = (start_times >= start_bound) & (end_times <= end_bound)

    return df.loc[condition].copy()


def dump_metadatafile(meta : dict, filename : str, format : str = "fobench"):
    """Saves the metadata dictionary as json file and manages which fileds to be saved depending on the format chosen (FoBench or FDSN)

    Parameters
    ----------
    meta : dict
        dicitonary which usually is the metadata.
    filename : str
        file name (with path) for the name of the produced metadata file.
    format : str, optional
        format on how to save the metadata. Options are the "fdsn" and "fobench". 
        For working in future projects and exploring file, use fobench as fdsn do not share databases. Default value is "fobench"
    """
    
    format = format.lower()
    
    if format not in ("fobench","fdsn"):
        raise ValueError('format not recognised. Use either "fobench" or "fdsn".')
    
    dump_meta = convert_types(meta)
    
    if format == "fdsn":
        dump_meta = fdsn_proj_check(dump_meta)

    with open(filename, 'w') as file: # finally dumping the metadata in the JSON file.
        json.dump(dump_meta, file, indent=4)
        
        
"""Validation methods of metadata for FDSN""" # this is killing my brains


def fdsn_fibre_check(fiber):
    """Validate and clean Fibre metadata for FDSN export as fiber.

    Parameters
    ----------
    fiber : dict
        Fibre metadata to validate and clean.

    Returns
    -------
    fiber : dict
        Fibre metadata prepared for FDSN export.
    """
    
    fiber_rev_fields = ["fiber_type", "fiber_reflector_interval", "fiber_reflector_interval_unit"]
    required_fiber_fields = ["fiber_id", "fiber_geometry", "fiber_mode", "fiber_refraction_index"]
    
    missing_fiber_fields = [field for field in required_fiber_fields if fiber.get(field) in (None, "", [])]
    
    # check errors
    if missing_fiber_fields:
        raise ValueError(f"Missing required FDSN Fibre fields: "
                        f"{missing_fiber_fields}")

    # remove FoBench-only Fibre fields.
    for field in fiber_rev_fields:
        fiber.pop(field, None)

    # remove unavailable optional values together with their units.
    optional_fiber_pairs = (("fiber_winding_angle", "fiber_winding_angle_unit"),
                            ("fiber_start_location", "fiber_start_location_unit"),
                            ("fiber_end_location", "fiber_end_location_unit"),
                            ("fiber_optic_length", "fiber_optic_length_unit"),
                            ("fiber_one_way_attenuation", "fiber_one_way_attenuation_unit"))

    for value_field, unit_field in optional_fiber_pairs:
        if fiber.get(value_field) in (None, "", "NA"):
            fiber.pop(value_field, None)
            fiber.pop(unit_field, None)
            
    return fiber


def fdsn_cable_check(cable):
    """Validate and clean Cable metadata for FDSN export.

    Parameters
    ----------
    cable : dict
        Cable metadata to validate and clean.

    Returns
    -------
    cable : dict
        Cable metadata prepared for FDSN export.
    """
    
    cable_rev_fields = ["cable_length", "cable_layers", "cable_layers_length", "n_fibres"]
    required_cable_fields = ["cable_id", "cable_bounding_box", "cable_owner"]
    
    missing_cable_fields = [field for field in required_cable_fields if cable.get(field) in (None, "", [])]
    
    # check errors
    if missing_cable_fields:
        raise ValueError(f"Missing required FDSN Cable fields: "
                        f"{missing_cable_fields}")
    if len(cable["cable_bounding_box"]) != 4:
        raise ValueError("FDSN cable_bounding_box must contain exactly four values.")

    # remove FoBench-only Cable fields.
    for field in cable_rev_fields:
        cable.pop(field, None)

    # empty dates are invalid date values in FDSN.
    for field in ("cable_installation_date", "cable_removal_date"):
        if cable.get(field) in (None, ""):
            cable.pop(field, None)

    # remove diameter and its unit when diameter is unavailable.
    if cable.get("cable_outside_diameter") in (None, "", "NA"):
        cable.pop("cable_outside_diameter", None)
        cable.pop("cable_outside_diameter_unit", None)

    # empty fibers array is invalid because it requires at least one entry.
    if not cable.get("fibers"):
        cable.pop("fibers", None)

    for fiber in cable.get("fibers", []):
        fdsn_fibre_check(fiber)
        
    return cable


def fdsn_chgroup_check(ch_group):
    """Validate Channel Group metadata for FDSN export.

    Parameters
    ----------
    ch_group : dict
        Channel Group metadata to validate.

    Returns
    -------
    ch_group : dict
        Channel Group metadata prepared for FDSN export.
    """

    required_group_fields = ["channel_group_id", "cable_id", "fiber_id", "coordinate_generation_date", "coordinate_system", "reference_frame",
                                    "distance_along_fiber_unit", "x_coordinate_unit", "y_coordinate_unit"]
    required_ch_fields = ["channel_ids", "distances_along_fiber", "x_coordinates", "y_coordinates"]

    missing_group_fields = [field for field in required_group_fields if ch_group.get(field) in (None, "", [])]
    
    if missing_group_fields:
        raise ValueError(f"Missing required FDSN Channel Group fields: "
                        f"{missing_group_fields}")

    if ch_group["coordinate_system"] not in ("geographic", "UTM", "local"):
        raise ValueError("Invalid FDSN coordinate_system: "
                        f'{ch_group["coordinate_system"]}')

    # channels are optional, but their arrays are required when the channels block is included.
    if "channels" in ch_group:

        channels = ch_group["channels"]

        if not isinstance(channels, dict):
            raise ValueError("FDSN channels must be a dictionary.")

        missing_ch_fields = [field for field in required_ch_fields if field not in channels]

        if missing_ch_fields:
            raise ValueError(f"Missing required FDSN Channel fields: "
                            f"{missing_ch_fields}")
    
    # remove unavailable optional uncertainties and their units. If there are None values, they should not show up in the metadata. Not necessary.
    uncertainty_pairs = (
        ("uncertainty_in_x_coordinate", "uncertainty_in_x_coordinate_unit"),
        ("uncertainty_in_y_coordinate", "uncertainty_in_y_coordinate_unit"),
        ("uncertainty_in_elevation", "uncertainty_in_elevation_unit"),
        ("uncertainty_in_depth", "uncertainty_in_depth_unit"),
        ("uncertainty_in_strike", "uncertainty_in_strike_unit"),
        ("uncertainty_in_dip", "uncertainty_in_dip_unit"),
    )

    for value_field, unit_field in uncertainty_pairs:
        if ch_group.get(value_field) in (None, "", "NA"):
            ch_group.pop(value_field, None)
            ch_group.pop(unit_field, None)
    
    # convert channel identifiers to the FDSN string format.
    for field in ("first_usable_channel_id", "last_usable_channel_id"):
        if ch_group.get(field) not in (None, ""):
            ch_group[field] = str(ch_group[field])

    if "channels" in ch_group:
        channels = ch_group["channels"]
        channels["channel_ids"] = [str(channel_id) for channel_id in channels["channel_ids"]]
            
    return ch_group


def fdsn_acqui_check(acqui):
    """Validate and clean Acquisition metadata for FDSN export.

    Parameters
    ----------
    acqui : dict
        Acquisition metadata to validate and clean.

    Returns
    -------
    acqui : dict
        Acquisition metadata prepared for FDSN export.
    """
    
    acqui_rev_fields = ["company", "sensing", "time_stamp", "channel_offset", "database_path", "database"]
    required_acqui_fields = ["acquisition_id", "acquisition_start_time", "acquisition_end_time", "acquisition_sample_rate", "acquisition_sample_rate_unit",
                                "gauge_length", "gauge_length_unit", "unit_of_measure", "number_of_channels", "spatial_sampling_interval", "spatial_sampling_interval_unit"]
    
    valid_fdsn_units = ("count", "m/m", "m/m/s", "m/s", "rad/s", "rad/m/s") # to be honest, i don't agree with mucch that is here but ok
    
    # VALIDATION PART -> Dataset
    missing_acqui_fields = [field for field in required_acqui_fields if acqui.get(field) in (None, "")]

    # check errors and missing fields
    if missing_acqui_fields:
        raise ValueError(f"Missing required FDSN Acquisition fields: "
                        f"{missing_acqui_fields}")
    if acqui["unit_of_measure"] not in valid_fdsn_units:
        raise ValueError("Invalid FDSN unit_of_measure: "
                        f'{acqui["unit_of_measure"]}')
        
    for field2 in acqui_rev_fields:
        # remove Dataset level fobench fields
        acqui.pop(field2, None)
        
    if acqui.get("scale_factor") in (None, "", "NA"):
        acqui.pop("scale_factor", None)

    if acqui.get("pulse_rate") in (None, "", "NA"):
        acqui.pop("pulse_rate", None)
        acqui.pop("pulse_rate_unit", None)

    if acqui.get("pulse_width") in (None, "", "NA"):
        acqui.pop("pulse_width", None)
        acqui.pop("pulse_width_unit", None)
        
    # VALIDATION PART -> Channel groups
    for ch_group in acqui.get("channel_groups", []):
        fdsn_chgroup_check(ch_group)
    
    return acqui


def fdsn_inter_check(inter):
    """Validate and clean Interrogator metadata for FDSN export.

    Parameters
    ----------
    inter : dict
        Interrogator metadata to validate and clean.

    Returns
    -------
    inter : dict
        Interrogator metadata prepared for FDSN export.
    """
    
    inter_rev_fields = ["sensing", "interrogator_path", "earliest_usage", "latest_usage", "n_files"]
    required_inter_fields = ["interrogator_id", "manufacturer", "model"]
    
    # VALIDATION PART -> Interrogator
    missing_inter_fields = [field for field in required_inter_fields if inter.get(field) in (None, "")]

    if missing_inter_fields:
        raise ValueError(f"Missing required FDSN Interrogator fields:" 
                        f"{missing_inter_fields}")
    
    for field1 in inter_rev_fields:
        
        # remove Interrogator level fobench fields
        inter.pop(field1, None)
        
    for acqui in inter.get("acquisitions", []):
        fdsn_acqui_check(acqui)    
    
    return inter


def fdsn_proj_check(proj):
    """Validate and clean Project metadata for FDSN export.

    Parameters
    ----------
    proj : dict
        Project metadata to validate and clean.

    Returns
    -------
    proj : dict
        Project metadata prepared for FDSN export.
    """
    
    proj_rev_fields = ["country", "end_date", "digital_object_identifier"]
    required_proj_fields = ["schema_version", "network_code", "location", "principal_investigator", "point_of_contact",
                            "point_of_contact_email", "point_of_contact_address", "start_date"]
    required_pi_fields = ["name", "email", "address"]
    
    # VALIDATION PART -> Project
    missing_proj_fields = [field for field in required_proj_fields if proj.get(field) in (None, "", [])] # Mandatory for FDSN!
    
    # check errors and missing fields that are mandatory
    if missing_proj_fields:
        raise ValueError(f"Missing required FDSN Project fields:"
                        f"{missing_proj_fields}")
        
    for forscher in proj["principal_investigator"]:
        if not isinstance(forscher, dict):
            raise ValueError("Each FDSN principal_investigator must be a dictionary.")

        missing_pi_fields = [field for field in required_pi_fields if forscher.get(field) in (None, "")]

        if missing_pi_fields:
            raise ValueError(f"Missing required FDSN Principal Investigator fields: "
                            f"{missing_pi_fields}")
    
    # removing Project level fobench fields
    proj.pop("start_time", None)
    proj.pop("end_time", None)

    # to prevent empty optional values to reach the json format and making it invalid.
    for field in proj_rev_fields:
        if proj.get(field) in (None, "", "NA"):
            proj.pop(field, None)

    if not proj.get("cables"):
        proj.pop("cables", None)
        
    # VALIDATION PART -> Cable - Fibre
    for cable in proj.get("cables", []):
        fdsn_cable_check(cable)

    if not proj.get("interrogators"):
        proj.pop("interrogators", None)
    
    # VALIDATION PART -> Interrogator - Dataset - ...
    for inter in proj.get("interrogators", []):
        fdsn_inter_check(inter)
        
    return proj