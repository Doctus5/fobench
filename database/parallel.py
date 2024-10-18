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
from itertools import chain
from functools import partial


class Parallel(object):
    '''
	IMPORTANT INFO: Nothing bro. We take it from here :).
	'''
    
    def __init__(self, params=None):
        '''
        Co-authors: --
        Description: 
            Initializes a "Parallel" Class which works as parallel wrapper for processes
        :Params:
            - params(type:Dictionary): dicitonary stating the variables for parallel programming.
        :Return:
            - NA.  
        '''
        
        # public attributes
        self.mode = params['mode'] # communication method / parallelisation method.
        self.n_cores = params['n_cores'] # number of cores to use in parallelisation.
        
        # private attributes
        
        # initialize any module and cores for parallel process.
        if self.mode == 'mpi':
            
            self.__import_mpi__()
            self.__init_mpi__()
        
        if self.mode == 'multiprocessing':
            
            self.__import_multiprocessing__()
            # self.__init_multiprocessing__()
    
    '''
	######################################################
	Private Functions
	######################################################
	'''

    def __import_mpi__(self):
        '''
        Co-authors: --
        Description: 
            Test and import necessary libraries for MPI.
        :Params:
            - NA.
        :Return:
            - NA.  
        '''
        
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
        '''
        Co-authors: --
        Description: 
            Initialize MPI modules
        :Params:
            - NA.
        :Return:
            - NA.  
        '''
        
        self.mpi_comm = self._MODE_.COMM_WORLD # Initialize parallelization communications.
        self.mpi_rank = self.mpi_comm.Get_rank() # recognize cores.
        self.mpi_size = self.mpi_comm.Get_size() # recognize available cores.
        self.n_cores = self.mpi_size
        
        print(f"MPI initialized. Core: {self.mpi_rank}, Size: {self.mpi_size}")
        
    
    def __import_multiprocessing__(self):
        '''
        Co-authors: --
        Description: 
            Test and import necessary libraries for multiprocessing.
        :Params:
            - NA.
        :Return:
            - NA.  
        '''
        
        # See if mpi4py is installed in the Python environment.
        try:
            
            # import concurrent.futures as multip
            import multiprocessing as multip
            self._MODE_ = multip
            
        except ImportError:
            
            raise ImportError("multiprocessing is not installed. Please install it by running 'pip install multiprocessing' or 'conda install multiprocessing'.")
        
    
    '''
	######################################################
	Public Functions
	######################################################
	'''

    # method for processing tasks.
    def submit(self, task, var_args, fixed_args, to_return=True):
        '''
        Co-authors: --
        Description: 
            Run a given task in parallel. The way it is done depends whereas the module used is MPI or Multiprocessing.
        :Params:
            - task(type:Function): function to use for calculations.
            - arguments(type:List or Tuple): arguments to use for the task to be parallelized.
        :Return:
            - NA.  
        '''
        
        # split the data to process into batches and rearange variables.
        var_args = np.array(var_args).astype(object)
        task_id = np.arange(var_args.size) # track id and order of processes.
        arguments_chunk = np.array_split(var_args, self.n_cores)
        task_id_chunk = np.array_split(var_args, self.n_cores)
        fixed_args = [fixed_args]
        
        # for MPI Protocol
        if self.mode == 'mpi':
            
            # distribute tasks to cores.
            var_batch = arguments_chunk[self.mpi_rank]
            # task_id_batch = task_id_chunk[self.mpi_rank]
            
            if var_batch.size != 0: # check that core is not processing empty tasks. Could bring problems.
            
                print(f'Processing in core {self._MODE_.Get_processor_name()}..{self.mpi_rank} a batch of {var_batch.size} tasks.')
                
                batch_result = task(var_batch, *fixed_args) # parallel task running.
                
                self.mpi_comm.Barrier() # wait for all to finish.
                results = self.mpi_comm.gather(batch_result, root=0)
                self.mpi_comm.Barrier() # wait for all to finish.
                
                if self.mpi_rank != 0:
                    
                    self._MODE_.Finalize()
                    exit()
            
        # for Multiprocessing.
        if self.mode == 'multiprocessing':
        
            results = []
            ids = np.arange(self.n_cores)
            
            with self._MODE_.Pool() as pool:
                
                results = pool.starmap(task, [(arguments_chunk[i], *fixed_args) for i in range(self.n_cores)])

        # Is there any output from the parallel computing for further process?
        if to_return == True:

            return results
        
    
    def __run_task__(self, task_id, task):
        
        return task_id, task