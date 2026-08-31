"""``Dataset`` class for visualizing, processing and handling a ``Dataset``.
A Dataset is understood as group of fiber optic sensing data ("ddss/das","dss","dts")
files part of a campaign. Can be continuous data in many files, and may have discontinuities.

:Authors:
    - Sergio Diaz-Meza
    - Jonas Pätzel

:Contributors:
    - Christopher Wollin

"""

import copy
import pandas as pd
from obspy.core import UTCDateTime as UTC

from .parallel import Parallel
from . import manager as manager
from .utils.windowing import dataset_windowing

class Dataset(object):
    """
    ``Dataset`` class which keeps track of the organization of the sfile(s)
    collected.

    Note
    ----
    Most of the methods perform changes within the class permanently.
    It can be useful to make a copy of the class instance with the
    ``.copy()`` method before performing any processing or changes.
    """

    def __init__(self, folder_path, company="silixa", sensing="das", database=None, metadata_file=None, storage_opts=None):
        """

        Parameters
        ----------
        filepath : str
            Folder/complete path where the data file(s) are located.
        company : str
                Manufacturer or the instrument that generates the data.
        sensing : str
            Specifies the type of fiber optic sensing technique of the data.

        """

        # internal attributed.
        self.__filepath__ = folder_path # filepath where the data is located (a folder).
        self.__built__ = False # is there a metadata file for it.
        self.__storage_opts__ = storage_opts

        # This variable might be redundant with the variables in metadata. Check this to reduce memory!
        self.database : pd.DataFrame = database # DataFrame format of available files corresponding to the curent Dataset.
        self.n_files = 0
        self.company = company # company of the manufacturer where the data comes from. Important to know how to read.
        self.sensing = sensing # sensing target of the dataset.
        self.start_time = None
        self.end_time = None
        self.sampling_rate = None
        self.dt = None
        self.gauge_length = None
        self.units = None #units of meassure.
        self.n_channels = None
        self.spatial_interval = None
        self.pulse_rate = None
        self.pulse_width = None
        self.channel_offset = 0

        self.metadata = self.__json_metadata__()

        if database is not None: # build from database (initial scanning).
            self.__fill_metadata__()
            self.__built__ = True # its now built.

        if metadata_file is not None: # build from metadata (metadata open).
            self.__build_from_metafile__(metadata_file)
            self.__built__ = True

    """Private Functions    """

    def __json_metadata__(self):
        """Define the metadata structure. Returns dict with the metadata parameters."""

        metadata = {
                        "Attributes": {
                            "acquisition_id": None,
                            "interrogator_id": None,
                            "acquisition_start_time": None,
                            "acquisition_end_time": None,
                            "acquisition_sample_rate": None,
                            "acquisition_sample_rate_unit": "Hz",
                            "time_stamp": None,
                            "gauge_length": None,
                            "gauge_length_unit": "meter",
                            "unit_of_measure": None,
                            "number_of_channels": None,
                            "spatial_sampling_interval": None,
                            "spatial_sampling_interval_unit": "meter",
                            "pulse_rate": "NA",
                            "pulse_rate_unit": "Hz",
                            "pulse_width": None,
                            "pulse_width_unit": "meter",
                            "comment": "NA",
                            "database_path": "NA"
                        },
                        "AttributeDefinitions": {
                            "acquisition_id": "Unique identifier of the data acquisition, assigned by data provider. Identifier should have a maximum of 8 alphanumeric characters with no special characters (e.g., underscores, period, dash). One identifier per acquisition settings.",
                            "interrogator_id": "Unique identifier of the interrogator unit used in this data acquisition. The acquisition_id must nest within this interrogator_id.",
                            "acquisition_start_time": "Start time of this data acquisition in UTC.",
                            "acquisition_end_time": "End time of this data acquisition in UTC. if data acquisition is still in operation, use a date in the future (e.g. 2999-01-01T00:00:00.000Z).",
                            "acquisition_sample_rate": "The rate at which the interrogator provides output data.",
                            "acquisition_sample_rate_unit": "Unit of acquisition sample rate.",
                            "time_stamp": "time stamp length.",
                            "gauge_length": "The averaging length along the fiber for a measurement, determined at experiment setup and used during acquisition.",
                            "gauge_length_unit": "Unit of gauge length.",
                            "unit_of_measure": "Unit of measure of archived data set. This may be the same unit as the Interrogator Unit of Measure if the data are raw.",
                            "number_of_channels": "The total number of sampling points along the fiber as output from the interrogator, referred to as NumberOfLoci in PRODML.",
                            "spatial_sampling_interval": "The channel spacing, or offset, between channels.",
                            "spatial_sampling_interval_unit": "Unit of spatial sampling interval.",
                            "pulse_rate": "Rate at which the interrogator unit interrogates the fiber sensor.",
                            "pulse_rate_unit": "Unit of pulse rate.",
                            "pulse_width": "Width of the pulse sent down the fiber in unit of time.",
                            "pulse_width_unit": "Unit of pulse width.",
                            "comment": "Additional comments.",
                            "database_path": "Central path fot he files of the Databse"
                        },
                        "AttributeRequirements": {
                            "acquisition_id": True,
                            "interrogator_id": True,
                            "acquisition_start_time": True,
                            "acquisition_end_time": True,
                            "acquisition_sample_rate": True,
                            "acquisition_sample_rate_unit": True,
                            "time_stamp": True,
                            "gauge_length": True,
                            "gauge_length_unit": True,
                            "unit_of_measure": True,
                            "number_of_channels": True,
                            "spatial_sampling_interval": True,
                            "spatial_sampling_interval_unit": True,
                            "pulse_rate": False,
                            "pulse_rate_unit": False,
                            "pulse_width": False,
                            "pulse_width_unit": False,
                            "comment": False,
                            "database_path": True,
                        },
                        "Database": None # DataFrame of parameters and location of files. Only for processing.
                    }

        return metadata

    def __metadates_2_isoformat__(self, reverse=False):
        """Transforms the dates of the ``Dataset`` files metadata into isoformat.
            Convenient for saving since Timestamp objects can not be saved in JSONs.

        Parameters
        ----------
        dataframe : pandas.Dataframe
            DataFrame object to find times and manipulate.
        reverse : bool
            If ``True``, transforms from isoformat to Timestamp object.

        """

        self.metadata["Database"] = manager.metadates_2_isoformat(self.metadata["Database"], reverse=reverse)

        return self

    def __build_from_metafile__(self, json_file):
        """Creates the basic variables of the DAS object with its characteristics.
        Builds from given metadata JSON file.

        Parameters
        ----------
        metadata : dict
            Dictionary of metadata.
        """

        if isinstance(json_file, str):
            meta_dict = manager.open_metadatafile(json_file)

        if isinstance(json_file, dict): # if the variable is already the dicitonary opened from Projects.
            meta_dict = json_file

        # Initialize Database
        self.metadata = meta_dict
        self.database = manager.init_dataframe(meta_dict["Database"]) # initialize the Dataframe of the database.
        self.__database_to_attributes__()

        # Fill attributes
        self.__filepath__ = meta_dict["Attributes"]["database_path"]

        return 0

    def __database_to_attributes__(self):
        """Fills Dataset attributes (build) from database DataFrame."""

        # Fill values in attributes
        self.n_files = self.database["file"].size
        self.start_time = UTC(self.database["start_time"].iloc[0])
        self.end_time = UTC(self.database["end_time"].iloc[-1])
        self.sampling_rate = self.database["sampling_rate"].iloc[0]
        self.dt = self.database["dt"].iloc[0]
        self.gauge_length = self.database["gauge_length"].iloc[0]
        self.units = None # units of meassure.
        self.n_channels = self.database["n_channels"].iloc[0]
        self.spatial_interval = self.database["spatial_interval"].iloc[0]
        self.pulse_rate = None
        self.pulse_width = None
        self.channel_offset = self.database["channel_offset"].iloc[0]

        return self

    def __fill_metadata__(self):
        """Define metadata structure for JSON. Fills the metadata arguments and
        returns them as dictionary.
        """

        self.__database_to_attributes__()

        self.metadata["Attributes"]["interrogator_id"] = None
        self.metadata["Attributes"]["acquisition_start_time"] = self.start_time.isoformat()
        self.metadata["Attributes"]["acquisition_end_time"] = self.end_time.isoformat()
        self.metadata["Attributes"]["acquisition_sample_rate"] = self.sampling_rate
        self.metadata["Attributes"]["acquisition_sample_rate_unit"] = "Hz"
        self.metadata["Attributes"]["time_stamp"] = self.dt
        self.metadata["Attributes"]["gauge_length"] = self.gauge_length
        self.metadata["Attributes"]["gauge_length_unit"] = "meter"
        self.metadata["Attributes"]["unit_of_measure"] = None
        self.metadata["Attributes"]["number_of_channels"] = self.n_channels
        self.metadata["Attributes"]["spatial_sampling_interval"] = self.spatial_interval
        self.metadata["Attributes"]["spatial_sampling_interval_unit"] = "meter"
        self.metadata["Attributes"]["pulse_rate"] = "NA"
        self.metadata["Attributes"]["pulse_rate_unit"] = "Hz"
        self.metadata["Attributes"]["pulse_width"] = self.pulse_width
        self.metadata["Attributes"]["channel_offset"] = self.channel_offset
        self.metadata["Attributes"]["pulse_width_unit"] = "meter"
        self.metadata["Attributes"]["comment"] = "NA"
        self.metadata["Attributes"]["database_path"] = self.__filepath__
        self.metadata["Database"] = manager.metadates_2_isoformat(self.database, reverse=False).to_dict(orient="list")

    """Public Functions"""

    def copy(self):
        """Returns a deep copy of the class in the moment of execution."""

        return copy.deepcopy(self)

    def build(self, format=None, parallels=None):
        """Builds from a given metadata file.

        Parameters
        ----------
        parallel_params : dict
            Dictionary containing the parameters for parallelization. If parameters are given
            then the building method runs in parallel. If ``None``, it runs in serial.

        Returns
        -------
        None

        """

        files = manager.scan_folder(self.__folder_path__, format=format, storage_opts=self.__storage_opts__) # remove the limitations.

        # calculate in parallel mode
        if parallels != None:
            hpc = Parallel(params=parallels) # initialize parallel process.
            results = hpc.submit(manager.files2database, files, (self.company, self.__storage_opts__)) # run.
            # now lets joint the results
            database_files = pd.concat(results, ignore_index=True)

        # in serial mode
        else:
            database_files = manager.files2database(files, self.company, self.__storage_opts__) # organize it as a Dataframe. Use Fiber Class.
        database_files = manager.chrono_order(database_files) # aranges in a chronological order.

        self.database = database_files
        self.__fill_metadata__()
        self.__built__ = True # its now built.

        return self

    def trim(self, time_range: tuple, include_overlap : bool = True):
        """Trim the dataset based on date ranges.
        See :func:`~fobench.database.manager.df_time_filtering`.

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
        None

        """

        self.database = manager.df_time_filtering(df=self.database, range=time_range, include_overlaps=include_overlap)
        self.n_files = len(self.database)

        if self.n_files == 0:

            self.start_time = None
            self.end_time = None
            return self

        self.start_time = self.database["start_time"].iloc[0]
        self.end_time = self.database["end_time"].iloc[-1]

        return self


    def window_map(self, time_range: tuple = None, window_size: float = None, step: float = None,
        include_overlap: bool = True, min_overlap_s: float = 0.0,
        group_cols: list[str] = None, include_meta: bool = False,
        return_windows: bool = False):
        """Build file-to-window mapping for this dataset.

        Parameters
        ----------
        time_range : tuple, optional
            Requested time range as ``(start_time, end_time)``.
            If ``None``, the full dataset range is used.
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
        group_cols : list of str, optional
            Optional grouping columns for compatibility grouping.
        include_meta : bool, optional
            If ``True``, merge file and window metadata into the output map.
        return_windows : bool, optional
            If ``True``, return both mapping and windows table.

        Returns
        -------
        out : pandas.DataFrame or tuple(pandas.DataFrame, pandas.DataFrame)
            Window mapping for this dataset. If ``return_windows=True``,
            returns ``(map_df, windows_df)``.

        """

        return dataset_windowing(
            dataset=self,
            time_range=time_range,
            window_size=window_size,
            step=step,
            include_overlap=include_overlap,
            min_overlap_s=min_overlap_s,
            group_cols=group_cols,
            include_meta=include_meta,
            return_windows=return_windows
        )


    def apply_fiber(self, task, output_store: str, output_key: str = None,
        time_range: tuple = None, window_size: float = None, step: float = None,
        include_overlap: bool = True, min_overlap_s: float = 0.0,
        group_cols: list[str] = None, fiber_kwargs: dict = None,
        parallel_params: dict = None, submit_chunk_size: int = None,
        cpu_ratio: float = 0.85, show_progress: bool = True,
        reducer: str = "mean", output_adapter=None):
        """Apply a Fiber task over the dataset windows and save to zarr.

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
        dictionary : dict
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
            merge_datasets=True,
            fiber_kwargs=fiber_kwargs,
            parallel_params=parallel_params,
            submit_chunk_size=submit_chunk_size,
            cpu_ratio=cpu_ratio,
            show_progress=show_progress,
            reducer=reducer,
            output_adapter=output_adapter
        )