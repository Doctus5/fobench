"""
Class "Unit" for visualizing and handling data within a project.
a Project is understood as a field campaing in an specific location where several Datasets are collected
from deployments.

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
from obspy.core import UTCDateTime as UTC

# Fobench classes
from .dataset import Dataset
from .parallel import Parallel

# Inner functions
from . import manager as manager
from .plotters import unit_plots as uni_plots
from .utils.windowing import unit_windowing



class Unit(object):
	'''
	IMPORTANT INFO: Most of the methods perform changes within the class permanently. Therefore is usefull to make a copy of the class
	with the method copy() before performing any processing or changes.
	'''
	
	#Creates the basic variables of the DAS object with its characteristics
	def __init__(self, folder_path=None, metadata_file=None, sensing='das', company="", format=""):
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
		self.datasets : list[Dataset] = [] # list of datasets. Each one is a Unit class that contains Datasets class.

		# Private attributes
		self.__folder_path__ = folder_path # central folder path of files
		self.__builded__ = False # is there a metadata file for it.

		# Public attributes
		self.sensing = sensing
		self.format = format
		self.company = company
		self.total_files = 0
		self.total_datasets = len(self.datasets) # total files that the unit produced.
		self.earliest_usage = None # earliest start date of meassurements with the unit.
		self.latest_usage = None # latest end date of meassurements with the unit.

		self.metadata = self.__json_metadata__()

		if metadata_file is not None: # check the case in which a metadata is introduced. Initialization from here is needed.

			self.__build_from_metafile__(metadata_file)
			self.__builded__ = True
    


	'''
	######################################################
	Private Functions
	######################################################
	'''
    
	def __str__(self):
		attributes = ['company', 'sensing', 'earliest_usage', 'latest_usage', 'total_datasets', 'total_files']

		return ('Unit class\n'
                'unit parameters:\n'
                f'{"-" * 65}\n'+ "\n".join(f"{attr.ljust(25)} = {getattr(self, attr)}" for attr in attributes))
    
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
						"earliest_usage": None,
						"latest_usage": None,
						"total_files": 0,
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
						"earliest_usage": "Earliest date of the datasets obtained with this unit for the project.",
						"latest_usage": "Latest date of the datasets obtained with this unit for the project.",
						"total_files": "Total number of files produced by this unit.",
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
						"earliest_usage": True,
						"latest_usage": True,
						"total_files": True,
						"model": True,
						"serial_number": False,
						"firmware_version": False,
						"comment": False,
						"interrogator_path": True
						},
					"Datasets": [] # list of metadata associated to datasets adquire at one interrogator unit.
                    }

		return metadata


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

		# Check if just the path of the metadata is being indicated.
		if isinstance(json_file, str):

			meta_dict = manager.open_metadatafile(json_file)

		if isinstance(json_file, dict): # if the variable is already the dicitonary opened from Projects.

			meta_dict = json_file

		# Initialize/Fill the Unit attributes.
		self.metadata = json_file
		self.__metadata_to_attributes__()

		# Initialize the Datasets.
		if meta_dict['Datasets']: # check Datasets.

			for meta_dataset in meta_dict['Datasets']:

				self.add_dataset( Dataset(self, metadata_file=meta_dataset) ) # Initialize the Units. )
    
			self.total_datasets = len(self.datasets)

		return self


	#Creates the basic variables of the DAS object with its characteristics
	def __metadata_to_attributes__(self):
		'''
		Co-authors: --
		Description: 
			Fills Unit attributes (build) from metadata.
		:Params:
			- NA.
		:Return:
			- NA.
		'''

		# Fill values in attributes
		self.__folder_path__ = self.metadata['Attributes']['interrogator_path'] # central folder path of files.

		self.sensing = self.metadata['Attributes']['sensing']
		self.company = self.metadata['Attributes']['manufacturer']
		self.total_files = self.metadata['Attributes']['total_files']
		self.earliest_usage = UTC(self.metadata['Attributes']['earliest_usage']) # earliest start date of meassurements with the unit.
		self.latest_usage = UTC(self.metadata['Attributes']['latest_usage']) # latest end date of meassurements with the unit.

		return self
    

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
		self.metadata['Attributes']["earliest_usage"] = self.earliest_usage.isoformat()
		self.metadata['Attributes']["latest_usage"] = self.latest_usage.isoformat()
		self.metadata['Attributes']["total_files"] = self.total_files
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


	# Adding Datasets to the Unit.
	def add_dataset(self, dataset):
		'''
		Co-authors: --
		Description: 
			Add Dataset objects to the current Unit Class.
		:Params:
			- dataset(type:Dataset Class): Dataset Class or Object to add to current Unit class. 
		:Return:
			- NA.
		'''

		self.datasets.append(dataset)

		return self


	#Creates the basic variables of the DAS object with its characteristics
	def build(self, parallels=None):
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

		files = manager.scan_folder(self.__folder_path__, format=self.format)#[:400] # remove the limitations.

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
			self.total_files = 0 # reset the variable to start summing.

			for dataset in self.datasets: # loop over existing datasets.

				earlier.append(dataset.start_time)
				later.append(dataset.end_time)
				self.total_files += dataset.total_files # adding to total number of files.
    
			self.earliest_usage = UTC(min(earlier))
			self.latest_usage = UTC(max(later))
			self.total_datasets = len(self.datasets)

		self.__fill_metadata__()
		self.__builded__ = True # its now builded.

		return self


	# Shows a table of the containing Datasets
	def dataset_table(self, return_table : bool = False, include_private: bool = False) -> None | pd.DataFrame:
		"""Prints or returns a summaized table of Datasets found for the unit and its properties.

		Parameters
		----------
		return_table : bool, optional
			if True, it returns the summarized table of datasets within Unit and its properties. If False, it prints it in terminal, by default False
		include_private : bool, optional
			includes private attributes of the datasets, by default False

		Returns
		-------
		None | DataFrame
			_description_
		"""

		rows = []

		for i, ds in enumerate(self.datasets):
    
			row = {}
 
			row["total_files"] = ds.total_files
			row["start_time"] = ds.start_time
			row["end_time"] = ds.end_time
			row["sampling_frequency"] = ds.sampling_frequency
			row["dt"] = ds.dt
			row["gauge_length"] = ds.gauge_length
			row["units"] = ds.units
			row["total_channels"] = ds.total_channels
			row["spatial_interval"] = ds.spatial_interval
			row["channel_offset"] = ds.channel_offset
    
			rows.append(row)

		df = pd.DataFrame(rows)

		if return_table:
			
			return df

		else:

			print(df.to_string(header=True))


	def trim(self, time_range: tuple, include_overlap : bool = True):
		"""Trim the unit and its datasets based on date ranges. See manager.fr_time_filtering for more details.

		Parameters
		----------
		time_range : tuple
			Start and end times in ISO format (or any datetime-parsable values).
			Example: ("2019-12-12T00:00:00Z", "2020-12-12T00:00:00Z")
		include_overlap : bool, optional
			If False (default), keep rows fully contained in ``range``.
			If True, keep rows that overlap ``range``. Deault is True.

		Returns
		-------
		_type_
			NA
		"""

		self.total_files = 0
		trimmed_datasets = []

		for ds in self.datasets:
    
			ds.trim(time_range, include_overlap)

			if ds.total_files > 0:
				trimmed_datasets.append(ds)
				self.total_files += ds.total_files

		self.datasets = trimmed_datasets
		self.total_datasets = len(self.datasets)

		# The following we can consider later. Early and latest usage mas night be related to the extend of the Dataset, but independent.
		if self.total_datasets > 0:
			self.earliest_usage = min(ds.start_time for ds in self.datasets)
			self.latest_usage = max(ds.end_time for ds in self.datasets)
		else:
			self.earliest_usage = None
			self.latest_usage = None

		return self


	def window_map(self, time_range: tuple = None, window_size: float = None, step: float = None,
		include_overlap: bool = True, min_overlap_s: float = 0.0,
		merge_datasets: bool = True, group_cols: list[str] = None,
		include_meta: bool = False,
		return_windows: bool = False):
		"""Build file-to-window mapping for this unit. Thsi is usually for the pipeline construction and parallel tasks.

		Parameters
		----------
		time_range : tuple, optional
			Requested time range as ``(start_time, end_time)``.
			If ``None``, the full unit range is used.
		window_size : float
			Window size in seconds.
		step : float, optional
			Step between consecutive windows in seconds. If ``None``,
			non-overlapping windows are used.
		include_overlap : bool, optional
			If ``True``, keep files that overlap the requested range.
			If ``False``, keep only files fully contained in the range.
		min_overlap_s : float, optional
			Minimum overlap (seconds) for file-window relations.
		merge_datasets : bool, optional
			If ``True``, all datasets in the unit are mapped into one single
			window timeline. If ``False``, each dataset is mapped independently
			and results are concatenated.
		group_cols : list of str, optional
			Columns to define compatibility groups.
		include_meta : bool, optional
			If ``True``, merge file and window metadata into the output map.
		return_windows : bool, optional
			If ``True``, return both mapping and windows table.

		Returns
		-------
		pandas.DataFrame or tuple[pandas.DataFrame, pandas.DataFrame]
			Window mapping for this unit. If ``return_windows=True``,
			returns ``(map_df, windows_df)``.
		"""

		return unit_windowing(
			unit=self,
			time_range=time_range,
			window_size=window_size,
			step=step,
			include_overlap=include_overlap,
			min_overlap_s=min_overlap_s,
			merge_datasets=merge_datasets,
			group_cols=group_cols,
			include_meta=include_meta,
			return_windows=return_windows
		)


	def apply_fiber(self, task, output_store: str, output_key: str = None,
		time_range: tuple = None, window_size: float = None, step: float = None,
		include_overlap: bool = True, min_overlap_s: float = 0.0,
		merge_datasets: bool = True, group_cols: list[str] = None,
		fiber_kwargs: dict = None, parallel_params: dict = None,
		submit_chunk_size: int = None, cpu_ratio: float = 0.85,
		show_progress: bool = True, reducer: str = "mean",
		output_adapter=None):
		"""Apply a Fiber task over unit windows and save to zarr.

		Parameters
		----------
		task : list, tuple, dict, str, or callable
			Fiber processing task definition.
		output_store : str
			Output zarr path.
		output_key : str, optional
			Output key to read from pipeline results.
		time_range : tuple, optional
			Time range as ``(start_time, end_time)``.
		window_size : float
			Window size in seconds.
		step : float, optional
			Window step in seconds. If ``None``, uses non-overlapping windows.
		include_overlap : bool, optional
			Keep files that overlap the requested range.
		min_overlap_s : float, optional
			Minimum overlap in seconds for file-window relations.
		merge_datasets : bool, optional
			If ``True``, map all datasets in one shared timeline.
		group_cols : list[str], optional
			Optional grouping columns for compatibility grouping.
		fiber_kwargs : dict, optional
			Keyword arguments passed to Fiber initialization.
		parallel_params : dict, optional
			Parallel backend parameters.
		submit_chunk_size : int, optional
			Number of files processed per submit wave.
		cpu_ratio : float, optional
			Automatic fraction of available CPUs used when ``n_cores`` is
			not explicitly set. Defaults to ``0.85``.
		show_progress : bool, optional
			Show progress bars for parallel execution.
		reducer : str, optional
			Reduction mode. Currently supports ``"mean"``.
		output_adapter : callable, optional
			Adapter ``adapter(value) -> 1D vector`` for operation output.

		Returns
		-------
		dict
			Execution summary and output metadata.
		"""

		# Lazy import to avoid requiring zarr unless apply_fiber is used.
		from .processing import run_window_pipeline_to_zarr

		# Thin wrapper: keep orchestration logic in shared processing module.
		return run_window_pipeline_to_zarr(
			source=self,
			task=task,
			output_store=output_store,
			output_key=output_key,
			time_range=time_range,
			window_size=window_size,
			step=step,
			include_overlap=include_overlap,
			min_overlap_s=min_overlap_s,
			group_cols=group_cols,
			merge_datasets=merge_datasets,
			fiber_kwargs=fiber_kwargs,
			parallel_params=parallel_params,
			submit_chunk_size=submit_chunk_size,
			cpu_ratio=cpu_ratio,
			show_progress=show_progress,
			reducer=reducer,
			output_adapter=output_adapter
		)


	'''
	####################################################
	Plotting functions below...
	####################################################
	'''


	def view_avail(self):
		'''
		Co-authors: --
		Description: 
			Function for plotting and viewing the available Datasets and their time coverage.
		:Params:
			- parallel_params(type:Dict): dictionary containing the parameters for parallelization. If parameters are
			given, then the building method runs in parallel. If None, it runs in serial.
		:Return:
			- NA.
		'''

		dataset_infos = [[dataset.start_time, dataset.end_time] for dataset in self.datasets]

		uni_plots.plot_data_coverage(dataset_infos, (self.earliest_usage, self.latest_usage))
