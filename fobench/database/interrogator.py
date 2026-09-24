"""Class ``Interrogator`` for visualizing and handling data within a ``Project``.
A ``Project`` is understood as a field campaing in an specific location where several
``Datasets`` are collected from deployments.

:Authors:
	- Sergio Diaz-Meza
	- Jonas Pätzel

:Contributors:
	- Christopher Wollin

"""

import copy
import pandas as pd
import numpy as np
from datetime import date
from obspy.core import UTCDateTime as UTC

from .dataset import Dataset
from .parallel import Parallel

from . import manager as manager
from .plotters import inter_plots as inter_plots
from .utils.windowing import inter_windowing
from .utils import utils



class Interrogator(object):
	"""Class keeps track of ``Datasets`` associated with an interrogator or sensing system.

	Note
	----
	Most of the methods perform changes within the class permanently.
	It can be useful to make a copy of the instance with the
	``.copy()`` method before performing any processing or changes.
	"""

	def __init__(self, folder_path: str = None, metadata_file: str = None,
			  sensing: str  = "das", company: str = "", format: str = "",
			  storage_opts=None):
		"""

		Parameters
		----------
		folder_path : str, optional
			Complete path where single/multiple data files are located.
		metadata_file : str, optional
			Complete path to where the metadata JSON file is located. This file is generated once
			the Interrogator class was run for the first time scanning files.
		sensing : str
			Fiber optic sensing technology. Currently only ``"das"``.
		company : str
			Interrogator manufacturer.
		format : str, optional
			Format specifier.
		storage_opts :
			DESCRIPTION.

		Returns
		-------
		None

		"""

		# In case class is initialized first time or empty (no metadata).
		self.datasets : list[Dataset] = [] # list of datasets. Each one is an Dataset class that contains files.

		# Private attributes
		self.__folder_path__ = folder_path # central folder path of files
		self.__built__ = False # is there a metadata file for it.
		self.__storage_opts__ = storage_opts # credentials of S3 to look at the bucket

		# Public attributes
		self.id = ""
		self.sensing = sensing
		self.format = format
		self.company = company
		self.n_files = 0
		self.n_datasets = len(self.datasets) # total datasets that the interrogator produced.
		self.earliest_usage = None # earliest start date of meassurements with the interrogator.
		self.latest_usage = None # latest end date of meassurements with the interrogator.

		self.metadata = self.__json_metadata__()

		if metadata_file is not None: # check the case in which a metadata is introduced. Initialization from here is needed.
			self.__build_from_metafile__(metadata_file)
			self.__built__ = True

	"""Private Methods"""

	def __str__(self):
		attributes = ['company', 'sensing', 'earliest_usage', 'latest_usage', 'n_datasets', 'n_files']

		return ('Interrogator class\n'
                'Interrogator parameters:\n'
                f'{"-" * 65}\n'+ "\n".join(f"{attr.ljust(25)} = {getattr(self, attr)}" for attr in attributes))
    # Define metadata structure for JSON.

	def __json_metadata__(self):
		"""Defines the metadata structure and returns dictionary with metadata parameters."""

		metadata = {
			# FDSN fields
			"interrogator_id": "",
			"manufacturer": "",
			"model": "",
			"serial_number": "",
			"firmware_version": "",
			"acquisitions": [],
			"comment": "",
   
			# Fobench fields
			"sensing": "",
			"interrogator_path": "",
			"earliest_usage": "",
			"latest_usage": "",
			"n_files": 0,
		}

		return metadata

	def __build_from_metafile__(self, json_file: str = None):
		"""Creates the basic variables of the DAS object with its characteristics.
		Builds from given metadata JSON file.
		"""

		# Check if just the path of the metadata is being indicated.
		if isinstance(json_file, str):
			meta_dict = manager.open_metadatafile(json_file)

		if isinstance(json_file, dict): # if the variable is already the dicitonary opened from Projects.
			meta_dict = json_file

		# Initialize/Fill the Interrogator attributes.
		self.metadata = meta_dict
		self.__metadata_to_attributes__()
  
		# Initialise the Datasets
		self.datasets = [Dataset(metadata_file=item) for item in meta_dict.get("acquisitions", [])]

		# Initialize the Datasets.
		# if meta_dict["acquisition"]: # check Datasets.

		# 	for meta_dataset in meta_dict["acquisition"]:
		# 		self.add_dataset( Dataset(self, metadata_file=meta_dataset) ) # Initialize the Interrogators. )
		# 	self.n_datasets = len(self.datasets)

		return self

	def __metadata_to_attributes__(self):
		"""Creates the basic variables of the FOS object with its characteristics.
		Fills Interrogator attributes (build) from metadata.
		"""

		# Fill values in attributes
		self.__folder_path__ = self.metadata.get("interrogator_path") # central folder path of files.
		self.id = self.metadata.get("interrogator_id")
		self.sensing = self.metadata.get("sensing") or "das"
		self.company = self.metadata["manufacturer"]
		self.n_files = self.metadata.get("n_files", 0)
		start, end  = self.metadata.get("earliest_usage"), self.metadata.get("latest_usage")
		self.earliest_usage = UTC(start) if start else None # earliest start date of meassurements with the interrogator.
		self.latest_usage = UTC(end) if end else None # latest end date of meassurements with the interrogator.

		return self

	def __fill_metadata__(self):
		"""Define metadata structure for JSON. Fills the metadata arguments and
		returns them as dictionary.
		"""

		# Fill values in metadata file
		self.metadata["interrogator_id"] = self.id
		self.metadata["manufacturer"] = self.company
		self.metadata["sensing"] = self.sensing
		self.metadata["earliest_usage"] = (self.earliest_usage.isoformat() + "Z" if self.earliest_usage is not None else "")
		self.metadata["latest_usage"] = (self.latest_usage.isoformat() + "Z" if self.latest_usage is not None else "")
		self.metadata["n_files"] = self.n_files
		self.metadata["interrogator_path"] = self.__folder_path__
		# self.metadata["model"] = 'NA'
		# self.metadata["serial_number"] = 'NA'
		# self.metadata["firmware_version"] = 'NA'
		self.metadata["acquisitions"] = [data_set.metadata for data_set in self.datasets] # populate with metadata

	def __metadates_2_isoformat__(self, reverse=False):
		"""Define metadata structure for JSON. Transforms the dates of the ``Dataset``
		files metadata into isoformat.Convenient for saving since Timestamp objects
		can not be saved in JSONs. If ``reverse=True``, transforma from isoformat
		to Timestamp object.
		"""

		if self.datasets:
			for dataset_meta in self.metadata["acquisition"]:
				dataset_meta["Database"] = manager.metadates_2_isoformat(dataset_meta["Database"], reverse=reverse)

		return self

	"""Public Functions"""

	def copy(self):
		"""Return a deep copy of the object. Useful for instances where there
		is no wish to affect the original data.
		"""
		return copy.deepcopy(self)


	def add_dataset(self, dataset):
		"""Adding ``Dataset`` to the ``Interrogator``."""
		self.datasets.append(dataset)

		return self

	def build(self, parallels=None):
		"""Creates the basic variables of the DAS object with its characteristics.
		Builds from a given metadata file.

		Parameters
		----------
		parallel_params : dict
			Dictionary containing the parameters for parallelization. If parameters are
			given, then the building method runs in parallel. If ``None``, it runs in serial.

		Returns
		-------
		None

		"""

		files = manager.scan_folder(self.__folder_path__, format=self.format, storage_opts=self.__storage_opts__)#[:400] # remove the limitations.

		# calculate in parallel mode
		if parallels != None:

			hpc = Parallel(params=parallels) # initialize parallel process.
			results = hpc.submit(manager.files2database, files, (self.company, self.__storage_opts__)) # run.
			database_files = pd.concat(results, ignore_index=True) # joint the results

		# calculate in serial mode
		else:

			database_files = manager.files2database(files, self.company, self.__storage_opts__) # organize it as a Dataframe. Use Fiber Class.

		database_files = manager.chrono_order(database_files) # aranges in a chronological order.
		chunks = manager.database_discontinuities(database_files, split=True) # splits based on discontinuities.

		# loop over the found chunks to initialize them as Datasets
		for chunk in chunks:

			self.datasets.append( Dataset(folder_path=self.__folder_path__, company=self.company, sensing=self.sensing, database=chunk, storage_opts=self.__storage_opts__) )

		if self.datasets:

			earlier, later = [], []
			self.n_files = 0 # reset the variable to start summing.

			for dataset_index, dataset in enumerate(self.datasets): # loop over existing datasets.

				dataset.id = str(dataset_index)
				dataset.update() # asignin id's
				earlier.append(dataset.start_time)
				later.append(dataset.end_time)
				self.n_files += dataset.n_files # adding to total number of files.

			self.earliest_usage = UTC(min(earlier))
			self.latest_usage = UTC(max(later))
			self.n_datasets = len(self.datasets)

		self.__fill_metadata__()
		self.__built__ = True # its now built.

		return self


	def update(self):
		"""Update Interrogator attributes and metadata without scanning files.

		Returns
		-------
		Interrogator
			Current updated Interrogator.
		"""

		self.n_files = 0
		self.n_datasets = len(self.datasets)

		for dataset_index, dataset in enumerate(self.datasets):
      
			dataset.id = str(dataset_index)

			for group_index, channel_group in enumerate(dataset.metadata["channel_groups"]):
       
				channel_group["channel_group_id"] = str(group_index)

			dataset.update()
			self.n_files += dataset.n_files

		if self.datasets:
      
			self.earliest_usage = min(dataset.start_time for dataset in self.datasets)
			self.latest_usage = max(dataset.end_time for dataset in self.datasets)
		else:
      
			self.earliest_usage = None
			self.latest_usage = None

		self.__fill_metadata__()
		self.__built__ = True

		return self


	def merge_datasets(self, max_gap:float):
		"""Merge consecutive compatible Datasets separated by short gaps.

		Parameters
		----------
		max_gap : float
			Maximum accepted interruption in seconds.

		Returns
		-------
		Interrogator
			Current Interrogator with its compatible Datasets merged.
		"""
		
		# secutiry checks
		if max_gap < 0:
			raise ValueError("max_gap cannot be negative.")
		if len(self.datasets) < 2: # no datasets to merge
			return self
    
		# these are the constants that must always remain so Dataset can be merged.
		constants = (
			"sampling_rate",
			"n_channels",
			"spatial_interval",
			"gauge_length",
			"channel_offset",
			"units",
			"scale_factor",
		)

		groups = [[self.datasets[0]]]

		# grouping adjacent compatible Datasets
		for dataset in self.datasets[1:]:
			
			previous = groups[-1][-1]
			gap = dataset.start_time - (previous.end_time + previous.dt)
			condition = all([getattr(dataset, constant) == getattr(previous, constant) for constant in constants]) # all must be True.

			if 0 <= gap <= max_gap and condition:
				groups[-1].append(dataset)
			else:
				groups.append([dataset])
    
		# concatenate
		merged_datasets = []
		for group in groups:

			dataset = group[0]
   
			if len(group) > 1:
				
				dataset.database = pd.concat([item.database for item in group], ignore_index=True)

			merged_datasets.append(dataset)
   
		self.datasets = merged_datasets
		self.update()

		return self


	def dataset_table(self, return_table : bool = False, include_private: bool = False) -> None | pd.DataFrame:
		"""Prints or returns a summarized table of ``Datasets`` found for the
		``Interrogator`` and its properties.

		Parameters
		----------
		return_table : bool, optional
			If ``True``, returns the summarized table of ``Datasets`` within ``Interrogator``
			and its properties. If ``False``, prints them in terminal.
		include_private : bool, optional
			Include private attributes of the ``Datasets``.

		Returns
		-------
		None | pandas.DataFrame

		"""

		rows = []
		for i, ds in enumerate(self.datasets):

			row = {}

			row["n_files"] = ds.n_files
			row["start_time"] = ds.start_time
			row["end_time"] = ds.end_time
			row["sampling_rate"] = ds.sampling_rate
			row["dt"] = ds.dt
			row["gauge_length"] = ds.gauge_length
			row["units"] = ds.units
			row["n_channels"] = ds.n_channels
			row["spatial_interval"] = ds.spatial_interval
			row["channel_offset"] = ds.channel_offset

			rows.append(row)

		df = pd.DataFrame(rows)

		if return_table:
			return df

		else:
			print(df.to_string(header=True))


	def trim(self, time_range: tuple, include_overlap : bool = True):
		"""Trim the interrogator and its datasets based on date ranges.
		See :func:`~fobench.database.manager.df_time_filtering`.

		Parameters
		----------
		time_range : tuple
			Start and end times in ISO format (or any datetime-parsable values).
			Example: ``("2019-12-12T00:00:00Z", "2020-12-12T00:00:00Z")``
		include_overlap : bool, optional
			If ``False`` keep rows fully contained in ``range``.
			If ``True``, keep rows that overlap ``range``.

		Returns
		-------
		None

		"""

		self.n_files = 0
		trimmed_datasets = []

		for ds in self.datasets:
			ds.trim(time_range, include_overlap)
			if ds.n_files > 0:
				trimmed_datasets.append(ds)
				self.n_files += ds.n_files

		self.datasets = trimmed_datasets
		self.n_datasets = len(self.datasets)

		# The following we can consider later. Early and latest usage mas night be related to the extend of the Dataset, but independent.
		if self.n_datasets > 0:
			self.earliest_usage = min(ds.start_time for ds in self.datasets)
			self.latest_usage = max(ds.end_time for ds in self.datasets)
		else:
			self.earliest_usage = None
			self.latest_usage = None

		return self


	'''Tools'''

	def append_coord(self, n_ch, x_ch, y_ch, z_ch, system, ref, which="all", coord_date=None):
		"""Attaches channel coordinates for later plotting. Takes 1D arrays of
		channel number (n_ch), longitude and latitude (x_ch and y_ch) and elevation in m (z_ch).
		"""

		datasets = utils.select_datasets(self, which)

		n_ch = np.asarray(n_ch, dtype=int)
		x_ch = np.zeros_like(n_ch, dtype=float) if x_ch is None else np.asarray(x_ch)
		y_ch = np.zeros_like(n_ch, dtype=float) if y_ch is None else np.asarray(y_ch)
		z_ch = np.zeros_like(n_ch, dtype=float) if z_ch is None else np.asarray(z_ch)
  
		if not (n_ch.size == x_ch.size == y_ch.size == z_ch.size):
			raise ValueError("Channel and coordinate arrays must have the same length.")
		if n_ch.size == 0:
			raise ValueError("At least one channel must be provided.")
		if np.unique(n_ch).size != n_ch.size:
			raise ValueError("Channel indices cannot contain duplicates.")

		# sorting channels and coordinates together. Just in case.
		order = np.argsort(n_ch)
		n_ch = n_ch[order]
		x_ch = x_ch[order]
		y_ch = y_ch[order]
		z_ch = z_ch[order]

		# now we split the coordinates wherever channel indices are discontinuous.
		split_indices = np.where(np.diff(n_ch) > 1)[0] + 1
		ch_sections = np.split(n_ch, split_indices)
		x_sections = np.split(x_ch, split_indices)
		y_sections = np.split(y_ch, split_indices)
		z_sections = np.split(z_ch, split_indices)

		coord_opts = {
			"decimal": ("geographic", "degree"),
			"geographic": ("geographic", "degree"),
			"utm": ("UTM", "m"),
			"local": ("local", "m")
		}

		system_key = system.lower()
		if system_key not in coord_opts:
			raise ValueError("system must be 'decimal', 'geographic', 'utm', or 'local'.")
		coord_system, coord_unit = coord_opts[system_key]

		if coord_date is None:
			coord_date = date.today().isoformat()

		for dataset in datasets:

			invalid_channels = n_ch[(n_ch < 0) | (n_ch >= dataset.n_channels)]

			if invalid_channels.size > 0:
				raise ValueError(f"Some provided channels do not match channels in Dataset {dataset.id}.")

			# preserving existing infrastructure references when replacing coordinates.
			cable_id, fiber_id = "", ""

			if dataset.metadata["channel_groups"]:
				cable_id, fiber_id = dataset.metadata["channel_groups"][0].get("cable_id", ""), dataset.metadata["channel_groups"][0].get("fiber_id", "")

			dataset.metadata["channel_groups"] = []

			# creating one Channel Group for every continuous channel section.
			for ch_section, x_section, y_section, z_section in zip(ch_sections, x_sections, y_sections, z_sections):

				dataset.add_ch_group(n_ch=ch_section, cable_id=cable_id, fiber_id=fiber_id)

				channel_group = dataset.metadata["channel_groups"][-1]
				channels = channel_group["channels"]

				channels["x_coordinates"], channels["y_coordinates"] = x_section, y_section
				channels["elevations_above_sea_level"] = z_section

				channel_group["coordinate_generation_date"] = coord_date
				channel_group["coordinate_system"] = coord_system
				channel_group["x_coordinate_unit"], channel_group["y_coordinate_unit"] = coord_unit, coord_unit
				channel_group["reference_frame"] = ref

		self.metadata["acquisitions"] = [dataset.metadata for dataset in self.datasets]

		return self


	def georeference(self, n_ch, x_ch, y_ch, z_ch, system="decimal", ref="WGS84", which="all", err=None, coord_date=None):
		"""Takes known channel locations, e.g. from tap tests and interpolates channel locations
		inbetween, attaches new coordinates.
		takes 1D arrays of channel number (n_ch), longitude and latitude (x_ch and y_ch)
		and elevation in m (z_ch), coordinate system can be for lon and lat can be "decimal" or "utm"
		"err" is maximum accepted interpolation error between original metadata location and new interpolated
		location
		"""
  
		datasets = utils.select_datasets(self, which)

		n_ch = np.asarray(n_ch)
		x_ch = np.zeros(n_ch.size) if x_ch is None else x_ch
		y_ch = np.zeros(n_ch.size) if y_ch is None else y_ch
		z_ch = np.zeros(n_ch.size) if z_ch is None else z_ch

		coords = None

		for dataset in datasets:

			coords = utils.interpolate_channels(n_ch, x_ch, y_ch, z_ch, system, err, dataset.spatial_interval)

		if coords is not None:

			self.append_coord(*coords, system=system, which=which, ref=ref, coord_date=coord_date)

		return self


	def window_map(self, time_range: tuple = None, window_size: float = None, step: float = None,
		include_overlap: bool = True, min_overlap_s: float = 0.0,
		merge_datasets: bool = True, group_cols: list[str] = None,
		include_meta: bool = False,
		return_windows: bool = False):
		"""Build file-to-window mapping for this ``Interrogator``. This is usually for the
		pipeline construction and parallel tasks.

		Parameters
		----------
		time_range : tuple, optional
			Requested time range as ``(start_time, end_time)``.
			If ``None``, the full interrogator range is used.
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
			If ``True``, all datasets in the interrogator are mapped into one single
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
		pandas.DataFrame | tuple(pandas.DataFrame, pandas.DataFrame)
			Window mapping for this interrogator. If ``return_windows=True``,
			returns ``(map_df, windows_df)``.

		"""

		return inter_windowing(
			inter=self,
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
		"""Apply a ``Fiber`` task over ``Interrogator`` windows and save to zarr.

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

	"""Plotting Functions"""

	def view_avail(self):
		"""Function for plotting and viewing the available Datasets and their time coverage."""

		dataset_infos = [[dataset.start_time, dataset.end_time] for dataset in self.datasets]
		inter_plots.plot_data_coverage(dataset_infos, (self.earliest_usage, self.latest_usage))