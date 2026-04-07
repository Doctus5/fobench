"""
Parallel wrapper methods for parallel computing and HPC for Fiber Optic projects.

Created on 2024-09-07 11:42:00
Last modification on 2024-09-07 11:42:00

:author:
	- Sergio Diaz-Meza (sergioad@gfz-potsdam.de)
:contributors:
	- Jonas Pätzel (jonas.patzel@ulb.be)
	- Christopher Wollin (wollin@gfz-potsdam.de)
:license:

"""

# Necessary libraries
import subprocess as terminal
import numpy as np
import sys


class Parallel(object):
    """Parallel wrapper class for FoBench database and processing tasks.

    Notes
    -----
    This class centralizes multiple backends under one interface:
    ``mpi``, ``multiprocessing``, ``dask`` and ``processpool``.
    """
    
    def __init__(self, params = None):
        """Initialize the parallel wrapper and selected backend.

        Parameters
        ----------
        params : dict, optional
            Configuration dictionary for parallel execution. Supported keys
            include:

            - ``mode`` : {"mpi", "multiprocessing", "dask", "processpool"}
            - ``n_cores`` : int
            - ``batch_size`` : int or None
            - Dask-specific options such as ``scheduler_address``,
              ``n_workers``, ``threads_per_worker``, ``processes``,
              ``memory_limit``, ``dashboard_address``, ``retries``.

        Raises
        ------
        ValueError
            If the selected ``mode`` is not supported.
        """
        
        params = {} if params is None else params

        # public attributes
        self.mode = params.get('mode', 'multiprocessing') # communication method / parallelisation method.
        self.n_cores = int(params.get('n_cores', 1)) # number of cores to use in parallelisation.
        self.batch_size = params.get('batch_size', None) # optional explicit task batch size. Amount of works that one worker/core must do.
        self.params = params

        # private attributes
        self._MODE_ = None
        self._dask_client = None
        self._dask_cluster = None
        self._dask_retries = 0

        # initialize any module and cores for parallel process.
        if self.mode == 'mpi':
            
            self.__import_mpi__()
            self.__init_mpi__()
        
        elif self.mode == 'multiprocessing':
            
            self.__import_multiprocessing__()
            # self.__init_multiprocessing__()

        elif self.mode == 'dask':

            self.__import_dask__()
            self.__init_dask__()

        elif self.mode in ('processpool', 'concurrent'):

            self.__import_concurrent__()

        else:

            raise ValueError(f'Unsupported parallel mode "{self.mode}". Use one of: mpi, multiprocessing, dask, processpool.')
    
    '''
	######################################################
	Private Functions
	######################################################
	'''

    def __import_mpi__(self):
        """Import and validate MPI dependencies.

        Raises
        ------
        ImportError
            If ``mpi4py`` is not available.
        EnvironmentError
            If ``mpirun`` (OpenMPI/MPICH runtime) is not available.
        """
        
        # See if mpi4py is installed in the Python environment.
        try:
            
            from mpi4py import MPI
            self._MODE_ = MPI
            
        except ImportError:
            
            raise ImportError("mpi4py is not installed. Please install it by running 'pip install mpi4py' or 'conda install mpi4py'.")
        
        # See if OpenMPI is installed/configured correctly in the system.
        try:
            
            terminal.run(["mpirun", "--version"], check=True, stdout=terminal.PIPE, stderr=terminal.PIPE)
            
        except terminal.CalledProcessError:
            
            raise EnvironmentError("OpenMPI is not installed or not configured correctly on your system.")
                
    
    def __init_mpi__(self):
        """Initialize MPI communicator and rank/size metadata.

        Notes
        -----
        After initialization, ``self.n_cores`` is set to the detected MPI
        world size.
        """
        
        self.mpi_comm = self._MODE_.COMM_WORLD # Initialize parallelization communications.
        self.mpi_rank = self.mpi_comm.Get_rank() # recognize cores.
        self.mpi_size = self.mpi_comm.Get_size() # recognize available cores.
        self.n_cores = self.mpi_size
        
        print(f"MPI initialized. Core: {self.mpi_rank}, Size: {self.mpi_size}")
        
    
    def __import_multiprocessing__(self):
        """Import Python multiprocessing backend.

        Raises
        ------
        ImportError
            If ``multiprocessing`` is unavailable.
        """
        
        # See if mpi4py is installed in the Python environment.
        try:
            
            # import concurrent.futures as multip
            import multiprocessing as multip
            self._MODE_ = multip
            
        except ImportError:
            
            raise ImportError("multiprocessing is not installed. Please install it by running 'pip install multiprocessing' or 'conda install multiprocessing'.")


    def __import_concurrent__(self):
        """Import ``concurrent.futures`` ProcessPool backend.

        Raises
        ------
        ImportError
            If ``concurrent.futures`` is unavailable.
        """

        try:

            from concurrent import futures as concurrent_futures
            self._MODE_ = concurrent_futures

        except ImportError:

            raise ImportError("concurrent.futures is not available in your Python environment.")


    def __import_dask__(self):
        """Import Dask distributed classes.

        Raises
        ------
        ImportError
            If Dask distributed is unavailable.
        """

        try:

            from dask.distributed import Client, LocalCluster
            self._MODE_ = {'Client': Client, 'LocalCluster': LocalCluster}

        except ImportError:

            raise ImportError("Dask distributed is not installed. Please install it by running 'pip install \"dask[distributed]\"' or 'conda install dask distributed'.")


    def __init_dask__(self):
        """Initialize Dask client.

        Notes
        -----
        If ``scheduler_address`` is provided, the class connects to an
        external scheduler. Otherwise, a local ``LocalCluster`` is created.
        """

        client_cls = self._MODE_['Client']
        cluster_cls = self._MODE_['LocalCluster']

        scheduler_address = self.params.get('scheduler_address', None)
        self._dask_retries = int(self.params.get('retries', 0))

        if scheduler_address is not None:

            self._dask_client = client_cls(scheduler_address)

        else:

            n_workers = int(self.params.get('n_workers', self.n_cores))
            threads_per_worker = int(self.params.get('threads_per_worker', 1))
            processes = bool(self.params.get('processes', True))
            memory_limit = self.params.get('memory_limit', 'auto')
            dashboard_address = self.params.get('dashboard_address', None)

            self._dask_cluster = cluster_cls(
                n_workers=n_workers,
                threads_per_worker=threads_per_worker,
                processes=processes,
                memory_limit=memory_limit,
                dashboard_address=dashboard_address
            )
            self._dask_client = client_cls(self._dask_cluster)
            self.n_cores = n_workers

        print(f'Dask initialized. Scheduler: {self._dask_client.scheduler_info()["address"]}')


    def __close_dask__(self):
        """Close Dask client and local cluster resources safely.

        Returns
        -------
        None
            Operates in-place and releases resources when available.
        """

        if self._dask_client is not None:

            self._dask_client.close()
            self._dask_client = None

        if self._dask_cluster is not None:

            self._dask_cluster.close()
            self._dask_cluster = None


    def __del__(self):
        """Finalize Dask resources during object destruction.

        Returns
        -------
        None
            Attempts best-effort cleanup.
        """

        if getattr(self, 'mode', None) == 'dask':

            try:

                self.__close_dask__()

            except Exception:

                pass


    def __normalize_fixed_args__(self, fixed_args):
        """Normalize fixed arguments into a tuple.

        Parameters
        ----------
        fixed_args : object, tuple, list, or None
            Fixed arguments that should be passed to each task invocation.

        Returns
        -------
        tuple
            Tuple-form fixed arguments.
        """

        if fixed_args is None:

            return tuple()

        if isinstance(fixed_args, tuple):

            return fixed_args

        if isinstance(fixed_args, list):

            return tuple(fixed_args)

        return (fixed_args,)


    def __chunk_arguments__(self, var_args, n_groups=None):
        """Split variable arguments into task batches.

        Parameters
        ----------
        var_args : array-like
            Variable arguments to split for parallel processing.
        n_groups : int, optional
            Target number of groups. If ``None``, ``self.n_cores`` is used.
            Ignored when ``batch_size`` is explicitly set.

        Returns
        -------
        list of list
            Batched arguments, each inner list corresponding to one task.
        """

        n_groups = self.n_cores if n_groups is None else int(n_groups)
        n_groups = max(1, n_groups)

        if not isinstance(var_args, np.ndarray):

            var_args = np.array(list(var_args), dtype=object)

        else:

            var_args = var_args.astype(object)

        if var_args.size == 0:

            return []

        if self.batch_size is not None:

            batch_size = max(1, int(self.batch_size))
            return [var_args[i:i+batch_size].tolist() for i in range(0, var_args.size, batch_size)]

        if var_args.size < n_groups:

            n_groups = var_args.size

        return [chunk.tolist() for chunk in np.array_split(var_args, n_groups) if chunk.size > 0]
        
    
    '''
	######################################################
	Public Functions
	######################################################
	'''

    # method for processing tasks.
    def submit(self, task, var_args, fixed_args, to_return=True):
        """Run a task in parallel using the configured backend.

        Parameters
        ----------
        task : callable
            Function to execute in parallel. It must accept a batch as first
            argument and fixed arguments afterwards.
        var_args : array-like
            Variable arguments to be partitioned into batches.
        fixed_args : object, tuple, list, or None
            Fixed arguments passed unchanged to each task invocation.
        to_return : bool, optional
            If ``True``, return aggregated backend results.

        Returns
        -------
        list or None
            List of backend task outputs when ``to_return`` is ``True``,
            otherwise ``None``.
        """
        
        # split the data to process into batches and rearrange variables.
        if not isinstance(var_args, np.ndarray):

            var_args = np.array(list(var_args), dtype=object)

        else:

            var_args = var_args.astype(object)

        fixed_args = self.__normalize_fixed_args__(fixed_args)

        if var_args.size == 0:

            return [] if to_return else None

        arguments_chunk = self.__chunk_arguments__(var_args)

        # for MPI Protocol
        if self.mode == 'mpi':

            arguments_chunk = self.__chunk_arguments__(var_args, n_groups=self.n_cores)

            # distribute tasks to cores.
            var_batch = arguments_chunk[self.mpi_rank] if self.mpi_rank < len(arguments_chunk) else []
            # task_id_batch = task_id_chunk[self.mpi_rank]

            batch_result = None
            if len(var_batch) != 0: # check that core is not processing empty tasks. Could bring problems.

                print(f'Processing in core {self._MODE_.Get_processor_name()}..{self.mpi_rank} a batch of {len(var_batch)} tasks.')

                batch_result = task(var_batch, *fixed_args) # parallel task running.

            self.mpi_comm.Barrier() # wait for all to finish.
            results = self.mpi_comm.gather(batch_result, root=0)
            self.mpi_comm.Barrier() # wait for all to finish.
            sys.stdout.flush() # Flush all output of other cores that are not the central.

            if self.mpi_rank != 0:

                self._MODE_.Finalize()
                exit()

            results = [result for result in results if result is not None]
            sys.stdout.flush() # flush out the central core messages
            
        # for Multiprocessing.
        elif self.mode == 'multiprocessing':
        
            results = []

            for i, batch in enumerate(arguments_chunk):

                print(f'Processing in core {i} a batch of {len(batch)} tasks.')

            with self._MODE_.Pool(processes=min(self.n_cores, len(arguments_chunk))) as pool:

                results = pool.starmap(task, [(batch, *fixed_args) for batch in arguments_chunk])

        # for Dask.
        elif self.mode == 'dask':

            for i, batch in enumerate(arguments_chunk):

                print(f'Processing in dask-worker-slot {i} a batch of {len(batch)} tasks.')

            futures = [
                self._dask_client.submit(
                    task,
                    batch,
                    *fixed_args,
                    retries=self._dask_retries,
                    pure=False
                )
                for batch in arguments_chunk
            ]
            results = self._dask_client.gather(futures)

        # for concurrent.futures ProcessPool backend.
        elif self.mode in ('processpool', 'concurrent'):

            results = []

            for i, batch in enumerate(arguments_chunk):

                print(f'Processing in processpool-worker-slot {i} a batch of {len(batch)} tasks.')

            with self._MODE_.ProcessPoolExecutor(max_workers=min(self.n_cores, len(arguments_chunk))) as executor:

                future_list = [executor.submit(task, batch, *fixed_args) for batch in arguments_chunk]
                results = [future.result() for future in future_list]

        # Is there any output from the parallel computing for further process?
        if to_return == True:

            return results
        
    
    def __run_task__(self, task_id, task):
        """Return task metadata as-is.

        Parameters
        ----------
        task_id : object
            Identifier associated with a task.
        task : object
            Task payload or result.

        Returns
        -------
        tuple
            Pair ``(task_id, task)``.
        """
        
        return task_id, task
