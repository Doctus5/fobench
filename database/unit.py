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
import copy
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
			- folder_path(type:String): complete path where the single/multiple data files are located.
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
		self.earliest_usage = None # earliest start date of meassurements with the unit.
		self.latest_usage = None # latest end date of meassurements with the unit.

		self.metadata = self.__json_metadata__()
    


	'''
	######################################################
	Private Functions
	######################################################
	'''

	#Creates the basic variables of the DAS object with its characteristics
	def __build_from_metafile__(self, json_file=None):
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
    
    
    # Define metadata structure for JSON.
	def __json_metadata__(self):
		'''
		Co-authors: --
		Description: 
			Define the metadata structure.
		:Params:
			- NA.
		:Return:
			- metadata(type:Dict): dictionary which is the metadata parameters.
		'''

		metadata = {
                    "Attributes": {
						"interrogator_id": None,
						"manufacturer": 'NA',
						"sensing": 'NA',
						"model": 'NA',
						"serial_number": None,
						"firmware_version": None,
						"comment": None,
						"interrogator_path": 'NA'
						},
					"AttributeDefinitions": {
						"interrogator_id": "Unique identifier of the interrogator unit used in the experiment, assigned by data provider. Identifier should have a maximum of 8 alphanumeric characters with no special characters (e.g., underscores, period, dash).",
						"manufacturer": "Manufacturer name of the interrogator.",
						"sensing": 'Sensing technique of the unit. Determines the type of data.',
						"model": "Model number of the interrogator.",
						"serial_number": "Serial number of the interrogator.",
						"firmware_version": "Firmware version of the software used within the interrogator.",
						"comment": "Additional comments",
						"interrogator_path": "Folder path of the files adquired with this interrogator or unit."
						},
					"AttributeRequirements": {
						"interrogator_id": True,
						"manufacturer": True,
						"sensing": True,
						"model": True,
						"serial_number": False,
						"firmware_version": False,
						"comment": False,
						"interrogator_path": True
						},
					"Datasets": [] # list of metadata associated to datasets adquire at one interrogator unit.
                    }

		return metadata 
    

    # Define metadata structure for JSON.
	def __fill_metadata__(self):
		'''
		Co-authors: --
		Description: 
			Fill the arguments on the metadata.
		:Params:
			- NA. 
		:Return:
			- metadata(type:Dict): dictionary which is the metadata parameters.
		'''

		# Fill values in metadata file
		self.metadata['Attributes']["interrogator_id"] = None
		self.metadata['Attributes']["manufacturer"] = self.company
		self.metadata['Attributes']["sensing"] = self.sensing
		self.metadata['Attributes']["interrogator_path"] = self.__folder_path__
		# self.metadata['Attributes']["model"] = 'NA'
		# self.metadata['Attributes']["serial_number"] = 'NA'
		# self.metadata['Attributes']["firmware_version"] = 'NA'
		self.metadata['Datasets'] = [data_set.metadata for data_set in self.datasets] # populate with metadata


    # Define metadata structure for JSON.
	def __metadates_2_isoformat__(self, reverse=False):
		'''
		Co-authors: --
		Description: 
			Transforms the dates of the Dataset files metadata into isoformat. 
			Convenient for saving since Timestamp objects can not be saved in JSONs.
		:Params:
			- reverse(type:Boolean): if True, it transform from isoformat to Timestamp object.
		:Return:
			- NA.
		'''

		if self.datasets:

			for dataset_meta in self.metadata['Datasets']:

				dataset_meta['Database'] = manager.metadates_2_isoformat(dataset_meta['Database'], reverse=reverse)

		return self
    

	'''
	######################################################
	Public Functions
	######################################################
	'''
 

	#Return a deep copy of the object. Useful for instances where there is no wish to affect the original data while keeping notherone affected.	
	def copy(self):
		'''
		Co-authors: --
		Description: Returns a deep copy of the class in the moment of execution.
		:Params:
			- NA.
		:Return:
			- (type:Unit Class): Same Unit Class in the state when the method is called.  
		'''
	
		return copy.deepcopy(self)


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

		files = manager.scan_folder(self.__folder_path__, format=format)[:400] # remove the limitations.

		# calculate in parallel mode
		if parallels != None:

			hpc = Parallel(params=parallels) # initialize parallel process.
			results = hpc.submit(manager.files2database, files, self.company) # run.
			database_files = pd.concat(results, ignore_index=True) # joint the results

		# calculate in serial mode
		else:

			database_files = manager.files2database(files, self.company) # organize it as a Dataframe. Use Fiber Class.

		database_files = manager.chrono_order(database_files) # aranges in a chronological order.
		chunks = manager.database_discontinuities(database_files, split=True) # splits based on discontinuities.

		# loop over the found chunks to initialize them as Datasets
		for chunk in chunks:

			self.datasets.append( Dataset(folder_path=self.__folder_path__, company=self.company, sensing=self.sensing, database=chunk) )

		if self.datasets:

			earlier, later = [], []

			for dataset in self.datasets: # loop over existing datasets.

				earlier.append(dataset.start_time)
				later.append(dataset.end_time)
    
			self.earliest_usage = min(earlier)
			self.latest_usage = max(later)

		self.__fill_metadata__()
		self.__builded__ = True # its now builded.

		return self