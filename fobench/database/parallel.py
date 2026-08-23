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
import time
from tqdm import tqdm


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

            from dask.distributed import Client, LocalCluster, as_completed
            self._MODE_ = {'Client': Client, 'LocalCluster': LocalCluster, 'as_completed': as_completed}

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
    
    def __submit_once__(self, task, var_args, fixed_args, to_return=True,
        show_progress=True, progress_desc='Parallel submit'):
        """Run one submission wave with current backend.

        Parameters
        ----------
        task : callable
            Function to execute in parallel.
        var_args : array-like
            Variable arguments for this wave.
        fixed_args : object, tuple, list, or None
            Fixed arguments passed to each task.
        to_return : bool, optional
            If ``True``, return wave results.
        show_progress : bool, optional
            If ``True``, display one progress bar for the wave.
        progress_desc : str, optional
            Description label for progress bar.

        Returns
        -------
        list or None
            Results for this submission wave when ``to_return=True``.
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

            batch_result = None
            if len(var_batch) != 0: # check that core is not processing empty tasks. Could bring problems.

                print(f'Processing in core {self._MODE_.Get_processor_name()}..{self.mpi_rank} a batch of {len(var_batch)} tasks.')

                batch_result = task(var_batch, *fixed_args) # parallel task running.

            self.mpi_comm.Barrier() # wait for all to finish.
            results = self.mpi_comm.gather(batch_result, root=0)
            self.mpi_comm.Barrier() # wait for all to finish.
            sys.stdout.flush() # flush output from non-root cores.

            if self.mpi_rank != 0:

                self._MODE_.Finalize()
                exit()

            results = [result for result in results if result is not None]
            sys.stdout.flush() # flush root core messages

        # for Multiprocessing.
        elif self.mode == 'multiprocessing':

            results = []

            if not show_progress:
                for i, batch in enumerate(arguments_chunk):
                    print(f'Processing in core {i} a batch of {len(batch)} tasks.')

            with self._MODE_.Pool(processes=min(self.n_cores, len(arguments_chunk))) as pool:

                # Submit each backend batch asynchronously to track global progress.
                async_results = [pool.apply_async(task, (batch, *fixed_args)) for batch in arguments_chunk]

                if show_progress:
                    pbar = tqdm(total=len(async_results), desc=progress_desc, leave=True, disable=not show_progress)
                    done_count = 0
                    while done_count < len(async_results):
                        new_done = sum(1 for result in async_results if result.ready())
                        if new_done > done_count:
                            pbar.update(new_done - done_count)
                            done_count = new_done
                        if done_count < len(async_results):
                            time.sleep(0.1)
                    pbar.close()

                results = [result.get() for result in async_results]

        # for Dask.
        elif self.mode == 'dask':

            if not show_progress:
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
            if show_progress:
                # Keep result order consistent with submission order.
                results = [None] * len(futures)
                future_idx = {future: i for i, future in enumerate(futures)}
                pbar = tqdm(total=len(futures), desc=progress_desc, leave=True, disable=not show_progress)
                for future, wave_result in self._MODE_['as_completed'](futures, with_results=True):
                    results[future_idx[future]] = wave_result
                    pbar.update(1)
                pbar.close()
            else:
                results = self._dask_client.gather(futures)

        # for concurrent.futures ProcessPool backend.
        elif self.mode in ('processpool', 'concurrent'):

            results = []

            if not show_progress:
                for i, batch in enumerate(arguments_chunk):
                    print(f'Processing in processpool-worker-slot {i} a batch of {len(batch)} tasks.')

            with self._MODE_.ProcessPoolExecutor(max_workers=min(self.n_cores, len(arguments_chunk))) as executor:

                future_map = {
                    executor.submit(task, batch, *fixed_args): i
                    for i, batch in enumerate(arguments_chunk)
                }
                results = [None] * len(future_map)

                if show_progress:
                    pbar = tqdm(total=len(future_map), desc=progress_desc, leave=True, disable=not show_progress)
                    for future in self._MODE_.as_completed(future_map):
                        idx = future_map[future]
                        results[idx] = future.result()
                        pbar.update(1)
                    pbar.close()
                else:
                    for future, idx in future_map.items():
                        results[idx] = future.result()

        # Is there any output from the parallel computing for further process?
        if to_return is True:

            return results
        
    
    '''
	######################################################
	Public Functions
	######################################################
	'''

    # method for processing tasks.
    def submit(self, task, var_args, fixed_args, to_return=True,
        submission_mode='single_submit', submit_chunk_size=None,
        show_progress=True):
        """Run a task in parallel using the configured backend.

        Parameters
        ----------
        task : callable or dict
            Task definition for parallel execution.

            - callable:
              Function to execute in parallel. It must accept a batch as
              first argument and fixed arguments afterwards.
            - dict:
              Pipeline configuration. Expected keys are:
              ``pipeline`` (or ``steps``), optional ``fiber_kwargs``, and
              optional ``validate_pipeline``.
        var_args : array-like
            Variable arguments to be partitioned into batches.
        fixed_args : object, tuple, list, or None
            Fixed arguments passed unchanged to each task invocation.
        to_return : bool, optional
            If ``True``, return aggregated backend results.
        submission_mode : {"single_submit", "chunked_submit"}, optional
            Submission strategy:
            - ``"single_submit"`` runs all ``var_args`` in one parallel call.
            - ``"chunked_submit"`` splits ``var_args`` into submission waves. Efficient for hughe memory limitations.
        submit_chunk_size : int, optional
            Number of ``var_args`` items per submission wave when
            ``submission_mode="chunked_submit"``. If ``None``, value is
            taken from ``params["submit_chunk_size"]``. If still ``None``,
            defaults to ``n_cores``.
        show_progress : bool, optional
            If ``True``, display one global progress bar for submit.

        Returns
        -------
        list or None
            List of backend task outputs when ``to_return`` is ``True``,
            otherwise ``None``.
        """

        # Resolve task mode:
        # - callable task: keep legacy behavior.
        # - dict task: treat as pipeline configuration.
        if callable(task):
            worker_task = task
            worker_fixed_args = fixed_args

        elif isinstance(task, dict):

            # Import only when needed to avoid unnecessary dependencies at import time.
            from .utils.pipeline import run_batch_pipeline, pipeline_check

            # Keep dictionary API simple and explicit.
            if ("pipeline" in task) and ("steps" in task):
                raise ValueError('Use only one of "pipeline" or "steps" in task dictionary.')

            if "pipeline" in task:
                pipeline = task["pipeline"]
            else:
                pipeline = task.get("steps", None)

            if pipeline is None:
                raise ValueError('Dictionary task requires "pipeline" (or "steps").')

            fiber_kwargs = task.get("fiber_kwargs", None)
            validate_pipeline = bool(task.get("validate_pipeline", False))

            # Validate once here, then disable repeated validation in worker batches.
            if validate_pipeline:
                pipeline = pipeline_check(pipeline, allow_callables=True)
                validate_pipeline = False

            # In dictionary task mode, fixed arguments are controlled by the dictionary.
            if fixed_args not in (None, (), []):
                raise ValueError(
                    '"fixed_args" is not used when "task" is a dictionary. '
                    'Put pipeline options inside task dict.'
                )

            worker_task = run_batch_pipeline
            worker_fixed_args = (pipeline, fiber_kwargs, validate_pipeline)

        else:
            raise TypeError('"task" must be callable or dictionary.')

        # Allow default strategy from params for easier user control.
        if submission_mode is None:
            submission_mode = self.params.get("submission_mode", "single_submit")

        if submission_mode not in ("single_submit", "chunked_submit"):
            raise ValueError(
                f'Unsupported submission_mode "{submission_mode}". '
                'Use one of: single_submit, chunked_submit.'
            )

        # Standard behavior: one submit call for all args.
        if submission_mode == "single_submit":
            return self.__submit_once__(
                task=worker_task,
                var_args=var_args,
                fixed_args=worker_fixed_args,
                to_return=to_return,
                show_progress=show_progress,
                progress_desc='Parallel submit'
            )

        # Chunked submission mode: run in waves to reduce memory pressure.
        if self.mode == "mpi":
            raise ValueError('chunked submittions are not currently supported for MPI mode.')

        if not isinstance(var_args, np.ndarray):
            var_args = np.array(list(var_args), dtype=object)
        else:
            var_args = var_args.astype(object)

        if var_args.size == 0:
            return [] if to_return else None

        # Resolve chunk size with a single public name.
        if submit_chunk_size is None:
            submit_chunk_size = self.params.get("submit_chunk_size", None)
        if submit_chunk_size is None:
            submit_chunk_size = self.n_cores

        submit_chunk_size = max(1, int(submit_chunk_size))

        # Build submission waves from the full argument list.
        wave_chunks = [
            var_args[i:i+submit_chunk_size]
            for i in range(0, var_args.size, submit_chunk_size)
        ]

        all_results = []
        pbar = tqdm(total=var_args.size, desc='Parallel submit', leave=True, disable=not show_progress)

        for i, wave in enumerate(wave_chunks):

            if not show_progress:
                print(f'Running submit wave {i+1}/{len(wave_chunks)} with {wave.size} tasks.')

            wave_result = self.__submit_once__(
                task=worker_task,
                var_args=wave,
                fixed_args=worker_fixed_args,
                to_return=to_return,
                show_progress=False
            )

            pbar.update(wave.size)

            if to_return and wave_result is not None:
                all_results.extend(wave_result)

        pbar.close()

        if to_return:
            return all_results
        
    
