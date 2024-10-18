"""
Class "Unit" for visualizing and handling data within a project.
a Project is understood as a field campaing in an specific location where several Datasets are collected
from deployments.
So far it recieves TDMS format (Silixa) and H5 format (Febus).

Created on 2022-08-19 12:07:17
Last modification on 2024-07-08 15:59:00

:author:
	- Sergio Diaz-Meza (sergioad@gfz-potsdam.de)
:contributors:
	- Jonas Pätzel (jonas.patzel@ulb.be)
	- Christopher Wollin (wollin@gfz-potsdam.de)
:license:

"""

# Necessary packages to import
import pandas as pd

# Fobench classes
from .dataset import Dataset
from .parallel import Parallel

# Inner functions
from . import manager as manager



class Unit(object):
	'''
	IMPORTANT INFO: Most of the methods perform changes within the class permanently. Therefore is usefull to make a copy of the class
	with the method copy() before performing any processing or changes.
	'''
	
	#Creates the basic variables of the DAS object with its characteristics
	def __init__(self, folder_path=None, metadata_file=None, sensing='das', company='silixa'):
		'''
		Co-authors: --
		Description: 
			Initializes a "Unit" Class which keeps track of Datasets associated to an interrogator or sensing system.
		:Params:
			- folder_path(type:String): compelte path where the single/multiple data files are located.
			- metadata_file(type:String): complete path where the metadata JSON file is located. This file is generated once
			the Unit class was run for the first time scanning files.
		:Return:
			- NA.  
		'''

		# In case class is initialized first time or empty (no metadata).
		self.datasets = [] # list of datasets. Each one is a Unit class that contains Datasets class.

		# Private attributes
		self.__folder_path__ = folder_path # central folder path of files
		self.__builded__ = False # is there a metadata file for it.

		# Public attributes
		self.sensing = sensing
		self.company = company

		# run initial build methods
		# if metadata_file == None:


		# else:
    


	'''
	######################################################
	Private Functions
	######################################################
	'''

	#Creates the basic variables of the DAS object with its characteristics
	def __build_from_metafile__(self, json_file=None, sensing='das'):
		'''
		Co-authors: --
		Description: 
			Builds from a given metadata file.
		:Params:
			- json_file(type:String): complete path where the single/multiple datasets are located.
		:Return:
			- NA.
		'''

		return 0
  

	'''
	######################################################
	Public Functions
	######################################################
	'''


	#Creates the basic variables of the DAS object with its characteristics
	def build(self, format=None, parallels=None):
		'''
		Co-authors: --
		Description: 
			Builds from a given metadata file.
		:Params:
			- parallel_params(type:Dict): dictionary containing the parameters for parallelization. If parameters are
			given, then the building method runs in parallel. If None, it runs in serial.
		:Return:
			- NA.
		'''

		files = manager.scan_folder(self.__folder_path__, format=format) # remove the limitations.

		# calculate in parallel mode
		if parallels != None:

			hpc = Parallel(params=parallels) # initialize parallel process.
			results = hpc.submit(manager.files2database, files, self.company) # run.
			
			# now lets joint the results
			database_files = pd.concat(results, ignore_index=True)

		# in serial mode
		else:

			database_files = manager.files2database(files, self.company) # organize it as a Dataframe. Use Fiber Class.

		database_files = manager.chrono_order(database_files) # aranges in a chronological order.
		chunks = manager.database_discontinuities(database_files, split=True) # splits based on discontinuities.

		# loop over the found chunks to initialize them as Datasets
		for chunk in chunks:
    
			self.datasets.append( Dataset(folder_path=self.__folder_path__, company=self.company, sensing=self.sensing, database=chunk) )

		self.__builded__ = True # its now builded.

		return self