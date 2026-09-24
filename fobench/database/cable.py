"""
Class ``Cable`` for listing the infrastructure to use in ``Project``.
A ``Project`` is understood as a field campaing in a specific location where
several ``Datasets`` are collected from deployments using the ``Cable``.

:Authors:
	- Sergio Diaz-Meza
	- Jonas Pätzel

:Contributors:
	- Christopher Wollin

"""

# Necessary packages to import
import copy
from obspy.core import UTCDateTime as UTC

# Fobench classes

# Inner functions
from .fibre import Fibre
from . import manager as manager



class Cable(object):
    """Class keeps track of the organization of a single/multiple data collected.

    Note
    ----
    Most of the methods perform changes within the class permanently.
    It can be useful to make a copy of the instance with the
    ``.copy()`` method before performing any processing or changes.
    """

    def __init__(self, metadata_file:str = None):
        """

        Parameters
        ----------
        metadata_file : json
                JSON object which represents the metadata of the project. If this file is non existent
                Cable must be initialized without it and a metadata must be created from scratch.

        Returns
        -------
        None

        """
        
        # private attributes
        self.__built__ = False
        
        
        self.id = ""
        self.owner = ""
        self.bounding_box = []
        self.start_time = None
        self.end_time = None
        self.model = ""
        self.env = ""
        self.diameter = None
        self.fibres : list[Fibre] = [] # list of cables used in the project
        self.n_fibres = len(self.fibres)
        
        self.length = None
        self.layers = []
        self.layers_len = []
        
        self.metadata = self.__json_metadata__()
        
        if metadata_file is not None: # check the case in which a metadata is introduced. Initialization from here is needed.

            self.__build_from_metafile__(metadata_file)
            self.__built__ = True
            
    
    """Private Functions"""
    
    def __str__(self):
        attributes = ["id", "start_time", "end_time", "model", "env", "length", "n_fibres"]

        return ('Cable class\n'
                'project parameters:\n'
                f'{"-" * 65}\n'+ "\n".join(f"{attr.ljust(25)} = {getattr(self, attr)}" for attr in attributes))
    
        
    def __json_metadata__(self):
        """Defines the metadata structure and returns dictionary with metadata parameters."""

        metadata = {
            # FDSN fields
            "cable_id": "",
            "cable_bounding_box": [],
            "cable_owner": "",
            "cable_installation_date": "",
            "cable_removal_date": "",
            "cable_characteristics": "",
            "cable_environment": "",
            "cable_installation_environment": "",
            "cable_model": "",
            "cable_outside_diameter": None,
            "cable_outside_diameter_unit": "",
            "fibers": [],
            "comment": "",

            # FoBench fields
            "cable_length": None,
            "cable_layers": [],
            "cable_layers_length": [],
            "n_fibres": 0
        }

        return metadata
    
    def __fill_metadata__(self):
        """Define metadata structure for JSON. Fill the arguments on the metadata.
        Returns dictionary with metadata parameters.
        """

        self.metadata["cable_id"] = str(self.id)
        self.metadata["cable_bounding_box"] = list(self.bounding_box)
        self.metadata["cable_owner"] = self.owner
        self.metadata["cable_installation_date"] = (self.start_time.date.isoformat() if self.start_time is not None else "")
        self.metadata["cable_removal_date"] = (self.end_time.date.isoformat() if self.end_time is not None else "")
        self.metadata["cable_environment"] = self.env
        # self.metadata["cable_installation_environment"] = ""
        self.metadata["cable_model"] = self.model
        self.metadata["cable_outside_diameter"] = self.diameter
        self.metadata["fibers"] = [glass.metadata for glass in self.fibres] # populate with fibres inside the cable
        
        # self.metadata["cable_outside_diameter_unit"] = "m"
        # self.metadata["comment"] = ""
        
        self.metadata["cable_length"] = self.length
        self.metadata["cable_layers"] = self.layers
        self.metadata["cable_layers_length"] = self.layers_len
        self.metadata["n_fibres"] = self.n_fibres
    
    
    def __build_from_metafile__(self, json_file=None):
        """Builds from a given metadata file. Needs  complete path where the single/multiple datasets are located."""

        # Check if just the path of the metadata is being indicated.
        if isinstance(json_file, str):
            meta_dict = manager.open_metadatafile(json_file)
        elif isinstance(json_file, dict): # if the variable is already the dicitonary opened from Projects.
            meta_dict = json_file
        else:
            raise TypeError("metadata_file must be a path string or dictionary")

        self.metadata = meta_dict
        self.__metadata_to_attributes__()

        # Initialise the fibers
        self.fibres = [Fibre(metadata_file=item) for item in meta_dict.get("fibers", [])]

        self.n_fibres = len(self.fibres)

        return self


    def __metadata_to_attributes__(self):
        """Fills Project attributes (build) from metadata."""

        # Fill values in attributes
        self.id = self.metadata["cable_id"]
        self.bounding_box = list(self.metadata.get("cable_bounding_box", []))
        self.owner = self.metadata["cable_owner"]

        start = self.metadata.get("cable_installation_date")
        end = self.metadata.get("cable_removal_date")
        self.start_time = UTC(start) if start else None
        self.end_time = UTC(end) if end else None

        self.model = self.metadata.get("cable_model", "")
        self.env = self.metadata.get("cable_environment", "")
        self.diameter = self.metadata.get("cable_outside_diameter")

        self.length = self.metadata.get("cable_length")
        self.layers = self.metadata.get("cable_layers", [])
        self.layers_len = self.metadata.get("cable_layers_length", [])

        return self
    
    
    """Public Functions"""
    
    
    def copy(self):
        """Return a deep copy of the object. Useful for instances where there
        is no wish to affect the original data.
        """
        return copy.deepcopy(self)


    def add_fibre(self, fibre):
        """Adding ``Fibre`` to the ``Cable``."""
        
        self.fibres.append(fibre)
        self.n_fibres = len(self.fibres)

        return self
    
    
    def build(self):
        """Builds the ``Cable`` object and parameters. Builds from a given metadata file.

        Parameters
        ----------
        parallel_params : dict
            Dictionary containing the parameters for parallelization. If parameters are
            given, then the building method runs in parallel. If ``None``, it runs in serial.

        Returns
        -------
        None

        """
        
        for fibre_index, fibre in enumerate(self.fibres):
            
            fibre.id = str(fibre_index)
            fibre.build()

        self.n_fibres = len(self.fibres)
        self.__fill_metadata__()
        self.__built__ = True

        return self
    
    