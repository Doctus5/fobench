"""
Class ``Fibre`` for listing the infrastructure to use in ``Cable``

:Authors:
    - Sergio Diaz-Meza
    - Jonas Pätzel

:Contributors:
    - Christopher Wollin

"""

import copy

from . import manager as manager
from .utils import utils



class Fibre(object):
    """
    ``Fibre`` class which keeps track of the fibres present in a cable.

    Note
    ----
    Most of the methods perform changes within the class permanently.
    It can be useful to make a copy of the class instance with the
    ``.copy()`` method before performing any processing or changes.
    """

    def __init__(self, metadata_file:str = None):
        """

        Parameters
        ----------
        metadata_file : json
                JSON object which represents the metadata of the project. If this file is non existent
                Fibre must be initialized without it and a metadata must be created from scratch.

        Returns
        -------
        None

        """

        # private attributes
        self.__built__ = False

        self.id = ""
        self.geometry = "linear"
        self.mode = "" # single-mode, multi-mode, other
        self.refr_index = 1.468 # fibre refractive index
        self.length = None
        
        self.type = "standard"
        self.refr_interval = None # just in case the type is engineered
        
        self.metadata = self.__json_metadata__()
        
        if metadata_file is not None: # check the case in which a metadata is introduced. Initialization from here is needed.

            self.__build_from_metafile__(metadata_file)
            self.__built__ = True

    """Private Functions    """
    
    def __str__(self):

        attributes = ['id', 'geometry', 'mode', 'type', 'refr_index', 'length', 'refr_interval']
        return ('Fibre class\n'
                'Fibre parameters:\n'
                f'{"-" * 65}\n'+ "\n".join(f"{attr.ljust(25)} = {getattr(self, attr)}" for attr in attributes))
    # Define metadata structure for JSON.

    def __json_metadata__(self):
        """Define the metadata structure. Returns dict with the metadata parameters."""

        metadata = {
            # FDSN fields
            "fiber_id": "",
            "fiber_geometry": "",
            "fiber_mode": "",
            "fiber_refraction_index": None,
            "fiber_winding_angle": None,
            "fiber_winding_angle_unit": "deg",
            "fiber_start_location": None,
            "fiber_start_location_unit": "m",
            "fiber_end_location": None,
            "fiber_end_location_unit": "m",
            "fiber_optic_length": None,
            "fiber_optic_length_unit": "m",
            "fiber_one_way_attenuation": None,
            "fiber_one_way_attenuation_unit": "dB/m",
            "comment": "",
            
            # Fobench fields
            "fiber_type": "", # standard, engineered
            "fiber_reflector_interval": None, # in case of engineered fibres
            "fiber_reflector_interval_unit": "m"
        }

        return metadata

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
        elif isinstance(json_file, dict): # if the variable is already the dicitonary opened from Projects.
            meta_dict = json_file
        else:
            raise TypeError("metadata_file must be a path string or a dictionary")

        # Fill attributes
        self.metadata = meta_dict

        self.id = self.metadata.get("fiber_id")
        self.geometry = self.metadata.get("fiber_geometry")
        self.mode = self.metadata.get("fiber_mode")
        self.refr_index = self.metadata.get("fiber_refraction_index")
        self.length = self.metadata.get("fiber_optic_length")

        # fobench fields
        self.type = self.metadata.get("fiber_type", "standard")
        self.refr_interval = self.metadata.get("fiber_reflector_interval")
        
        return self


    def __fill_metadata__(self):
        """Define metadata structure for JSON. Fills the metadata arguments and
        returns them as dictionary.
        """

        self.metadata["fiber_id"] = str(self.id)
        self.metadata["fiber_geometry"] = self.geometry
        self.metadata["fiber_mode"] = self.mode
        self.metadata["fiber_refraction_index"] = self.refr_index
        self.metadata["fiber_optic_length"] = self.length
        
        self.metadata["fiber_type"] = self.type
        self.metadata["fiber_reflector_interval"] = self.refr_interval
        
        self.metadata["fiber_one_way_attenuation"] = utils.power_loss(self.mode, self.type, self.refr_interval) # calculations assuming some values for the moment.
        
        return self

    """Public Functions"""

    def copy(self):
        """Returns a deep copy of the class in the moment of execution."""

        return copy.deepcopy(self)

    def build(self):
        """Builds from a given metadata file.

        Parameters
        ----------

        Returns
        -------
        None

        """

        self.__fill_metadata__()
        self.__built__ = True # its now built.

        return self