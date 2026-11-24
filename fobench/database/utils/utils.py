"""Some utilities to handle things in Project and subsequent levels."""

import numpy as np

from fobench.core.tools.utils import interpolate_channels



def select_datasets(inter, which="all"):
    """Select Datasets using their acquisition identifiers.

    Parameters
    ----------
    inter : Interrogator
    	Interrogator containing the Datasets to select.
    which : str, int or array-like, optional
        Acquisition IDs to select.

    Returns
    -------
    list
        Selected Dataset instances.
    """

    if isinstance(which, str) and which == "all":
        return inter.datasets

    if np.isscalar(which):
        selected_ids = {str(which)}
    else:
        selected_ids = {str(dataset_id) for dataset_id in which}

    datasets = [dataset for dataset in inter.datasets if dataset.metadata["acquisition_id"] in selected_ids]

    found_ids = {dataset.metadata["acquisition_id"] for dataset in datasets}
    missing_ids = selected_ids - found_ids

    if missing_ids:
        raise ValueError(f"Dataset IDs not found: {sorted(missing_ids)}")

    return datasets


def power_loss(mode:str = "single-mode", fibre_type:str = "standard", refr_interval:float = 0):
	"""power_loss _summary_

	Parameters
	----------
	fibre_type : str, optional
		Type of fibre. Can be ``"standard"`` or ``"engineered"``.
	refr_interval : float, optional
		In case of engineered fibre, it specifies the spacing between engineered reflectors.

	Returns
	-------
	loss : float
		Total power loss per metre

	"""

	# security checks
	if mode not in ("single-mode","multi-mode"):
		raise ValueError("mode must be 'single-mode' or 'multi-mode'")
	if fibre_type not in ("standard","engineered"):
		raise ValueError("fibre_type must be 'standard' or 'engineered'")
	if fibre_type == "engineered" and (refr_interval <= 0 or refr_interval == None):
		raise ValueError("refr_interval must be greater than 0")

	loss = 2e-4 if mode == "single-mode" else 2.5e-3 # dB/m, loss of standard telecom fibre (single-mode or multi-mode)
	eng_loss = 4.34e-4 # power loss of transmitted light due to engineered reflector

	if fibre_type == "engineered":

		loss += eng_loss/refr_interval

	return loss