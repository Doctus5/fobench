"""
Class "Project" for visualizing and handling data within a project.
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
import json
from obspy.core import UTCDateTime as UTC

# Fobench classes

# Inner functions
from .unit import Unit
from . import manager as manager



class Project(object):
	'''
	IMPORTANT INFO: Most of the methods perform changes within the class permanently. Therefore is usefull to make a copy of the class
	with the method copy() before performing any processing or changes.
	'''
	
	#Creates the basic variables of the Project object with its characteristics
	def __init__(self, folder_path=None, metadata_file=None):
		'''
		Co-authors: --
		Description: 
			Initializes a "Project" Class which keeps track of the organization of the single/multiple data collected.
		:Params:
			- metadata_file(type:json): JSON object which represents the metadata of the project. If this file is non existent
			then Project must be initialized without it and a metadata must be created from scratch (see Tutorial).
			- project_type(type:String): type of project to handle based on 
		:Return:
			- NA.  
		'''

		# Private attributes
		self.__folder_path__ = folder_path # central folder path of files
		self.__builded__ = False # is there a metadata file for it.

		# Public attributes
		
		# In case class is initialized first time or empty (no metadata).
		self.units = [] # list of units/interrogators used. Each one is a Unit class that contains Datasets class.
		self.metadata = self.__json_metadata__()

		if metadata_file is not None: # check the case in which a metadata is introduced. Initialization from here is needed.

			self.__build_from_metafile__(metadata_file)
			self.__builded__ = True

	'''
	######################################################
	Private Functions
	######################################################
	'''
    
    
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
						"network_code": "NA",
						"location": "NA",
						"country": "NA",
						"principal_investigator_name": "NA",
						"principal_investigator_email": "NA",
						"principal_investigator_address": "NA",
						"point_of_contact": "NA",
						"point_of_contact_email": "NA",
						"point_of_contact_address": "NA",
						"start_date": "NA",
						"end_date": "NA",
						"funding_agency": "NA",
						"project_number": "NA",
						"digital_object_identifier": "NA",
						"purpose_of_data_collection": "NA",
						"comment": None
						},
					"AttributeDefinitions": {
						"network_code": "Unique network name for the installation with a maximum of 8 alphanumeric characters with no special characters (e.g., underscores, period, dash).",
						"location": "Name of the geographic location of the installation.",
						"country": "Country where the installation is located. Use ISO 3166-1 alpha-3 three-letter country code.",
						"principle_investigator_name": "Name of principal investigator (last name, first name) for the installation.",
						"principle_investigator_email": "Email address of principal investigator.",
						"principle_investigator_address": "Physical address and institution of principal investigator.",
						"point_of_contact": "Point of contact (last name, first name) for the metadata.",
						"point_of_contact_email": "Email address of point of contact.",
						"point_of_contact_address": "Physical address and institution of point of contact.",
						"start_date": "Start date of data collection at the installation in UTC.",
						"end_date": "End date of data collection at the installation in UTC. If installation is still in operation, use a future date (e.g. 2999-01-01).",
						"funding_agency": "Name(s) of agency that funded the experiment.",
						"project_number": "Funding project number. Should be supplied if a number has been assigned by funding agency(s).",
						"digital_object_identifier": "Digital Object Identifier that uniquely identifies the metadata, this identifier may only become available following archiving.",
						"purpose_of_data_collection": "Brief explanation of the purpose of the experiment.",
						"comment": "Additional comments."
						},
					"AttributeRequirements": {
						"network_code": True,
						"location": True,
						"country": True,
						"principle_investigator_name": True,
						"principle_investigator_email": True,
						"principle_investigator_address": True,
						"point_of_contact": True,
						"point_of_contact_email": True,
						"point_of_contact_address": True,
						"start_date": True,
						"end_date": True,
						"funding_agency": True,
						"project_number": True,
						"digital_object_identifier": True,
						"purpose_of_data_collection": True,
						"comment": False
						},
					"Interrogator": []
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
		# self.metadata['Attributes']["model"] = 'NA'
		# self.metadata['Attributes']["serial_number"] = 'NA'
		# self.metadata['Attributes']["firmware_version"] = 'NA'
		self.metadata['Interrogator'] = [unit.metadata for unit in self.units] # populate with metadata


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

		# Check if just the path of the metadata is being indicated.
		if isinstance(json_file, str):

			meta_dict = manager.open_metadatafile(json_file)

		if isinstance(json_file, dict): # if the variable is already the dicitonary opened from Projects.

			meta_dict = json_file

		self.metadata = meta_dict
		self.__metadata_to_attributes__()

		if meta_dict['Interrogator']:

			for mes_unit in meta_dict['Interrogator']:

				ind_unit = Unit(self, metadata_file=mes_unit) # Initialize the Units.
				self.add_unit( ind_unit )

		return self


	def __metadata_to_attributes__(self):
		'''
		Co-authors: --
		Description: 
			Fills Project attributes (build) from metadata.
		:Params:
			- NA.
		:Return:
			- NA.
		'''

		# Fill values in attributes
		# self.__folder_path__ = self.metadata['Attributes']['interrogator_path']

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
			- (type:Project Class): Same Project Class in the state when the method is called.  
		'''
	
		return copy.deepcopy(self)

	# Adding Units to the Project.
	def add_unit(self, unit):
		'''
		Co-authors: --
		Description: 
			Add units class to the current Project object.
		:Params:
			- unit(type:Unit Class): Unit Class or Object to add to current Project class. 
		:Return:
			- NA.
		'''

		self.units.append(unit)

		return self

	# Builds the Project object and parameters.
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

		self.__fill_metadata__()
		self.__builded__ = True # its now builded.

		return self


	def save_metadata(self, filename='project_meta.json'):
		'''
		Co-authors: --
		Description: 
			Saves the metadata file for future usage and not building from scratch the project.
		:Params:
			- filename(type:String): file name with complete path and format of the metadata file. If not given, Default = 'project_meta.json',
			which means it is saved in the local folder of code execution.
		:Return: 
			- NA.
		'''

		# Save metadata.
		with open(filename, 'w') as file:

			dump_metadata = manager.convert_types(self.metadata)
			json.dump(dump_metadata, file, indent=4)


