"""
Method for managing, sorting and manipulating files in folders.

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
import pandas as pd

# inside packaged
from ..fobench.fiber import Fiber


def scan_folder(filepath, format=None):
    '''
    Co-authors: --
    Description: 
        Function that scans all within a folder to find all files inside.
    :Params:
        - filepath(type:String): compelte path fot he file to be read.
        - format(type:String, Optional): if a format extension is specified (without the point), only the files with such an extension will be returned.
    :Return:
        - files(type:List): list of paths of each of the files.  
	'''
    
    format = '*' if format == None else '*.' + format + '*'
    files = glob.glob(filepath + format)
    
    return files


def chrono_order(files, company):
    '''
    Co-authors: --
    Description: 
        Function that aranges all fiber optic sensing files in a chronological way. Requires reading the file without loading the data.
    :Params:
        - filepath(type:String): compelte path fot he file to be read.
        - company(type:String): name of the manufacturer.
    :Return:
        - files(type:List): list of paths of each of the files.  
	'''

    N = len(files) # number of data files.
    filtered_keys = ['start_time', 'end_time', 'dt', 'sampling_frequency', 'total_channels', 
                     'spatial_interval', 'gauge_lenght', 'channel_offset'] # attributes of interest for holding consistency.
    #filtered_keys = ['start_time', 'end_time', 'dt'] # attributes of interest for holding consistency.
    
    database = []
    
    # scan the files
    for i in range(N):
    
        file = files[i] # select file
        # by checking the start and end times in the files
        d_file = Fiber(file, company=company, load_data=False)
        
        info = d_file.metadata(meta_dict=True) # getting public relevant attributes.
        info = {key: info[key] for key in filtered_keys if key in info}  # Filter
        info['file'] = d_file.__filepath__ # we also take the filepath.
        
        # correct format for start and end times into strings (UTC also can do the thing).
        info['start_time'] = info['start_time'].isoformat()
        info['end_time'] = info['end_time'].isoformat()
        
        database.append(info)
        del d_file # free memory.
        
    # handling it as a DataFrame
    info = pd.DataFrame.from_dict(database)
    info = info.sort_values(by=['start_time', 'end_time'])
    
    return 0