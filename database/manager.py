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
    
    for i, file in enumerate(files):
    
        d_file = Fiber(file, company=company, load_data=True)
        file_start, file_end = d_file.start_time, d_file.end_time
        
        print(i, file_start, file_end)
    
    #return 0