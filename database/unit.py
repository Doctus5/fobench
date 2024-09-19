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
		self.__folder_path__ = folder_path

		# Public attributes
		self.sensing = sensing
		self.company = company

		# run initial build methods
		if metadata_file == None:

			self.__build_from_scratch__(parallel=True)

		else:
    
			self.__build_from_metafile__()


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


	#Creates the basic variables of the DAS object with its characteristics
	def __build_from_scratch__(self, parallel=False):
		'''
		Co-authors: --
		Description: 
			Builds from a given metadata file.
		:Params:
			- folder_path(type:String): compelte path where the raw files of a single unit are located.
		:Return:
			- NA.
		'''

		files = manager.scan_folder(self.__folder_path__)[100:400] # remove the limitations.
		# print('SERIAL!!!')
		# database_files = manager.files2database(files, self.company) # organize it as a Dataframe. Use Fiber Class.
		# print(database_files, type(database_files))

		# calculate in parallel mode
		print('PARALLEL!!!')
		if parallel == True:

			hpc = Parallel(params={'mode':'mpi', 'n_cores':50}) # initialize parallel process.
			results = hpc.submit(manager.files2database, files, self.company) # run.
			# print(results, type(results))
			# now lets joint the results

		return 0

		database_files = manager.chrono_order(database_files) # aranges in a chronological order.
		chunks = manager.database_discontinuities(database_files, split=True) # splits based on discontinuities.

		# loop over the found chunks to initialize them as Datasets
		for chunk in chunks:
    
			self.datasets.append( Dataset(folder_path=self.__folder_path__, company=self.company, sensing=self.sensing, database=chunk) )

		return self