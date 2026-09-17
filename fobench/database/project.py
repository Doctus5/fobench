"""Class ``Project`` for visualizing and handling data within a ``Project``.
A ``Project`` is understood as a field campaing in a specific location where
several ``Datasets`` are collected from deployments.

:Authors:
	- Sergio Diaz-Meza
	- Jonas Pätzel

:Contributors:
	- Christopher Wollin

"""

# Necessary packages to import
import copy
import json
from obspy.core import UTCDateTime as UTC

# Fobench classes

# Inner functions
from .interrogator import Interrogator
from . import manager as manager



class Project(object):
	"""Class keeps track of the organization of a single/multiple data collected.

	Note
	----
	Most of the methods perform changes within the class permanently.
	It can be useful to make a copy of the instance with the
	``.copy()`` method before performing any processing or changes.
	 """

	def __init__(self, folder_path = None, metadata_file:str = None):
		"""

		Parameters
		----------
		folder_path : str
			Path.
		metadata_file : json
			 JSON object which represents the metadata of the project. If this file is non existent
			 Project must be initialized without it and a metadata must be created from scratch.

		Returns
		-------
		None

		"""

		# Private attributes
		self.__folder_path__ = folder_path # central folder path of files
		self.__built__ = False # is there a metadata file for it.

		# Public attributes
		# In case class is initialized first time or empty (no metadata).
		self.network_code = ""
		self.location = ""
		self.country = ""
		self.inters : list[Interrogator] = [] # list of interrogators used. Each one is an Interrogator class that contains Datasets class.
		self.n_inters = len(self.inters)
		self.start_time = None
		self.end_time = None

		self.metadata = self.__json_metadata__()

		if metadata_file is not None: # check the case in which a metadata is introduced. Initialization from here is needed.

			self.__build_from_metafile__(metadata_file)
			self.__built__ = True

	"""Private Functions"""

	def __str__(self):
		attributes = ['location', 'n_inters', 'start_time', 'end_time']

		return ('Project class\n'
				'project parameters:\n'
				f'{"-" * 65}\n'+ "\n".join(f"{attr.ljust(25)} = {getattr(self, attr)}" for attr in attributes))

	def __json_metadata__(self):
		"""Define the metadata structure. Returns dict with the metadata parameters."""

		metadata = {
			# FDSN fields
			"schema":"https://www.fdsn.org/schemas/DAS-Metadata-FDSN/2.0",
			"schema_version": "2.0",
			"network_code": "",
			"location": "",
			"country": "", # string
			"principal_investigator": [{"name":"", "email":"", "address":""}],
			"point_of_contact": "",
			"point_of_contact_email": "",
			"point_of_contact_address": "",
			"start_date": "", # string
			"end_date": "", # string
			"funding_agency": "",
			"project_number": "",
			"digital_object_identifier": "",
			"purpose_of_data_collection": "",
			"comment": "",
			"interrogators": [],
			"cables": [],

			# Fobench fields
			"start_time": "",
			"end_time": ""
		}

		return metadata

	def __fill_metadata__(self):
		"""Define metadata structure for JSON. Fill the arguments on the metadata.
		Returns dictionary with metadata parameters.
		"""

		self.metadata["network_code"] = self.network_code
		self.metadata["location"] = self.location
		self.metadata["country"] = self.country
		self.metadata["start_date"] = self.start_time.date.isoformat()
		self.metadata["start_time"] = self.start_time.isoformat() + "Z"
		self.metadata["end_date"] = (self.end_time.date.isoformat() if self.end_time is not None else "")
		self.metadata["end_time"] = (self.end_time.isoformat() + "Z" if self.end_time is not None else "")
		self.metadata["interrogators"] = [inter.metadata for inter in self.inters]


	def __build_from_metafile__(self, json_file=None):
		"""Builds from a given metadata file. Needs  complete path where the
		single/multiple datasets are located.
		"""

		# Check if just the path of the metadata is being indicated.
		if isinstance(json_file, str):
			meta_dict = manager.open_metadatafile(json_file)
		elif isinstance(json_file, dict): # if the variable is already the dicitonary opened from Projects.
			meta_dict = json_file
		else:
			raise TypeError("metadata_file must be a path string or dictionary")

		self.metadata = meta_dict
		self.__metadata_to_attributes__()

		# Initialise the Interrogators
		self.inters = [Interrogator(metadata_file=item) for item in meta_dict.get("interrogators", [])]

		# if meta_dict['interrogators']:
		# 	for mes_inter in meta_dict['interrogators']:
		# 		ind_inter = Interrogator(self, metadata_file=mes_inter) # Initialize the Interrogators.
		# 		self.add_inter(ind_inter)

		self.n_inters = len(self.inters)

		return self


	def __metadata_to_attributes__(self):
		"""Fills Project attributes (build) from metadata."""

		# Fill values in attributes
		# self.__folder_path__ = self.metadata['Attributes']['interrogator_path']
		self.network_code = self.metadata["network_code"]
		self.location = self.metadata["location"]
		self.country = self.metadata.get("country") or ""
		self.start_time = UTC(self.metadata.get("start_time") or self.metadata["start_date"])
		end_time = self.metadata.get("end_time") or self.metadata.get("end_date")
		self.end_time = UTC(end_time) if end_time else None

		return self

	"""Public Functions"""

	def copy(self):
		"""Returns a deep copy of the class in the moment of execution."""

		return copy.deepcopy(self)

	def add_inter(self, inter: Interrogator):
		"""Add ``Interrogator`` class to the current ``Project`` object."""

		self.inters.append(inter)
		self.n_inters = len(self.inters)

		return self

	def build(self, parallels=None):
		"""Builds the ``Project`` object and parameters. Builds from a given metadata file.

		Parameters
		----------
		parallel_params : dict
			Dictionary containing the parameters for parallelization. If parameters are
			given, then the building method runs in parallel. If ``None``, it runs in serial.

		Returns
		-------
		None

		"""

		start_time_list, end_time_list = [], []

		for inter_index, inter in enumerate(self.inters):

			inter.metadata["interrogator_id"] = str(inter_index)
			inter.build(parallels=parallels)
			start_time_list.append(inter.earliest_usage)
			end_time_list.append(inter.latest_usage)

		self.start_time, self.end_time = min(start_time_list), max(end_time_list)

		self.__fill_metadata__()
		self.__built__ = True

		return self

	def save_metadata(self, filename : str = "project_meta.json", format : str = "fobench"):
		"""Saves the metadata file for future usage and toin order to having to
		build the project again.

		Parameters
		----------
		filename : str
			File name with complete path and format of the metadata file.
			If not given, Default = 'project_meta.json', which means it is saved
			in the local folder of code execution.
		format : str
			Format of the metadata file. If not given, Default = "fobench" is used, and useful for fobench utilities (recommended).
			Option "fdsn" gives the current standard fields required by FDSN. Useful for staying within the standard format, 
			but throws out filepaths necessary to keep track of the databases.

		Returns
		-------
		None
		"""

		with open(filename, 'w') as file:
			dump_metadata = manager.convert_types(self.metadata)
			json.dump(dump_metadata, file, indent=4)