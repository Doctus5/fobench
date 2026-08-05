"""
Method for managing, sorting and manipulating files and databases in folders.

Created on 2024-07-08 16:50:00
Last modification on 2024-07-08 16:50:00

:author:
	- Sergio Diaz-Meza (sergioad@gfz-potsdam.de)
:contributors:
	- Jonas Pätzel (jonas.patzel@ulb.be)
	- Christopher Wollin (wollin@gfz-potsdam.de)
:license:

"""

# Necessary packages to import
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
    from fobench.fiber import Fiber
except ModuleNotFoundError:
    from ..fobench.fiber import Fiber


def scan_folder(folder_path, format=None, storage_opts=None):
    '''
    Function that scans all within a folder to find all files inside.

    Parameters
    ----------
    folder_path : str
        Global folder path where files will be searched.
    format : str, optional
        If a format extension is specified (without the point), only the files with such an extension will be returned.
    storage_opts : dict, optional
        Options passed to the S3 filesystems (profile, endpoint, anonymous access settings, etc.).
    
    Returns
    -------
    list of str
        Sorted local paths or S3 UIRs.
    '''
    
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
    
    if not files:
        
        warnings.warn(f'No files with the specified format .{format} or files at all were found.', category=UserWarning)
                
    return files


def files2database(files, company, storage_opts=None):
    '''
    Co-authors: --
    Description: 
        Function reada all the FOS files with Fiber Class and reports a big Database from them with all its important metadata.
    :Params:
        - files(type:String): list of compelte paths from the files to be handled.
        - company(type:String): name of the manufacturer.
    :Return:
        - database(type:Dataframe): dataframe indicating paths and essential metadata from each file.  
	'''

    N = len(files) # number of data files.
    filtered_keys = ['start_time', 'end_time', 'dt', 'sampling_frequency', 'total_channels', 
                     'spatial_interval', 'gauge_length', 'channel_offset'] # attributes of interest for holding consistency.
    #filtered_keys = ['start_time', 'end_time', 'dt'] # attributes of interest for holding consistency.
    
    database = []
    
    # scan the files
    for i in range(N):
    
        file = files[i] # select file
        # by checking the start and end times in the files
        print('File: ' + file)
        d_file = Fiber(file, company=company, load_data=False, storage_opts=storage_opts)
        
        info = d_file.metadata(meta_dict=True) # getting public relevant attributes.
        info = {key: info[key] for key in filtered_keys if key in info}  # Filter
        info['file'] = d_file.__filepath__ # we also take the filepath.
        
        # correct format for start and end times into strings (UTC also can do the thing).
        info['start_time'] = pd.to_datetime(info['start_time'].isoformat())
        info['end_time'] = pd.to_datetime(info['end_time'].isoformat())
        
        database.append(info)
        del d_file # free memory.
        
    # handling it as a DataFrame
    database = pd.DataFrame.from_dict(database)
    
    return database


def init_dataframe(dictionary):
    '''
    Co-authors: --
    Description: 
        Initialize a DataFrame from dictionary and converts the dates in isoformat to datetime for management. 
    :Params:
        - dictionary(type:Dict): database from metadata as dictionary.
    :Return:
        - database(type:DataFrame): DataFrame version of the database.  
	'''

    database = pd.DataFrame(dictionary)
    database = metadates_2_isoformat(database, reverse=True)
    
    return database


def open_metadatafile(json_file):
    '''
    Co-authors: --
    Description: 
        Opens a metadata JSON file. 
    :Params:
        - json_file(type:String): dataframe indicating paths and essential metadata from each file.
    :Return:
        - metadata(type:Dict): metadata as dictionary.  
	'''

    with open(json_file, 'r') as file:
        
        metadata = json.load(file)
    
    return metadata


def chrono_order(database):
    '''
    Co-authors: --
    Description: 
        Function that aranges all fiber optic sensing files in a chronological way. Requires reading the file without loading the data.
    :Params:
        - database(type:Dataframe): dataframe indicating paths and essential metadata from each file.
    :Return:
        - database(type:Dataframe): dataframe indicating paths and essential metadata from each file.  
	'''

    database = database.sort_values(by=['start_time', 'end_time'])
    
    return database


def database_discontinuities(df, split=True):
    '''
    Co-authors: --
    Description: 
        Checks if there are time and parameters discontinuities along a reported Dataset. Necessary to separate chuncks of data.
        In case the user wants it, it returns the database splitted 
    :Params:
        - df(type:Dataframe): database indicating paths and essential metadata from each file.
    :Return:
        - continuity(type:DataFrame): if split == True, the function returns a boolean DataFrame where it shows where a
        dicontinuity happens. False objects mark the beginning of a dicontinuity along the data.
        - database_chunks(type:List): if split == False, it returns a List of dataframes indicating paths and essential 
        metadata from each file.
	'''
    
    # attributes of metadata to evaluate continuity.
    var_conditions = ['sampling_frequency', 'total_channels', 'spatial_interval', 'gauge_length', 'channel_offset']
    
    # Initialize continuity check by checking time.
    # End time of file i + dt must match the start time of file i+1 
    df['start_time'] = pd.to_datetime(df['start_time'])
    df['end_time'] = pd.to_datetime(df['end_time'])
    
    # check chrono continuity
    expected_start = df['end_time'] + pd.to_timedelta(df['dt'], unit='s')
    continuity = expected_start.shift(1) == df['start_time']
    continuity.iloc[0] = True
    continuity.iloc[-1] = True

    # # check discontinuities in others properties.
    for variable in var_conditions[1:]:
        
        continuity *= df[variable].shift(1, fill_value=df[variable].iloc[0]) == df[variable]
    
    # Evaluate if the user wants the databse to be splitted in discontinuous parameters (see selected parameters).
    if split == False:
        
        return continuity
    
    else:
        
        database_chunks = []
        last_incident = 0
        
        for i in range(continuity.size):
            
            # Evaluate when there is a chunk.
            if (continuity.iloc[i] == False) or (i == continuity.size-1):
                
                database_chunks.append( df.iloc[last_incident:i] ) # we add the detected chunk.
                last_incident = i
        
        return database_chunks
    

def metadates_2_isoformat(dataframe, reverse=False):
    '''
    Co-authors: --
    Description: 
        Transforms the dates of the Dataset files metadata into isoformat. 
        Convenient for saving since Timestamp objects can not be saved in JSONs.
    :Params:
        - dataframe(type:DataFrame): DataFrame object to find times and manipulate.
        - reverse(type:Boolean): if True, it transform from isoformat to Timestamp object.
    :Return:
        - dataframe(typeDataFrame): modified Dataframe object.
    '''
    
    # Find these start time and end time columns in the DataFrame.
    if reverse == False:
    
        dataframe['start_time'] = dataframe['start_time'].apply(lambda x: x.isoformat())
        dataframe['end_time'] = dataframe['end_time'].apply(lambda x: x.isoformat())
        
    else:
        
        dataframe['start_time'] = pd.to_datetime(dataframe['start_time'])
        dataframe['end_time'] = pd.to_datetime(dataframe['end_time'])
        
    return dataframe
    
    
# Recursive function to convert types in a dictionary
def convert_types(input_dict):
    
    output_dict = {}
    for key, value in input_dict.items():
        if isinstance(value, dict):
            # Recursively convert nested dictionaries
            output_dict[key] = convert_types(value)
        elif isinstance(value, np.int64):
            output_dict[key] = int(value)  # Convert to standard int
        elif value is None:
            output_dict[key] = None  # Keep as None (will be serialized as null in JSON)
        elif isinstance(value, list):
            # Convert any int64 or None in lists
            output_dict[key] = [convert_types(v) if isinstance(v, dict) else int(v) if isinstance(v, np.int64) else v for v in value]
        else:
            output_dict[key] = value  # Keep other types as is
    return output_dict


# Return filtered DataFrame in time
def df_time_filtering(df: pd.DataFrame, range: tuple, start_label: str = "start_time",
    end_label: str = "end_time", include_overlaps: bool = False) -> pd.DataFrame:
    """Filter a DataFrame by a requested time range.

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
    pd.DataFrame
        Filtered DataFrame.
    """

    if len(range) != 2:
        raise ValueError("`range` must contain exactly two values: (start_time, end_time).")

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
    
