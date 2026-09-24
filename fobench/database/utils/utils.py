"""Some utilities to handle things in Project and subsequent levels."""

import numpy as np

from pyrocko.orthodrome import distance_accurate50m_numpy as deg2m



def select_datasets(inter, which="all"):
    """Select Datasets using their acquisition identifiers.

    Parameters
    ----------
    inter : Interrogator
    	Interrogator containing the Datasets to select.
    which : str, int or array-like, optional
        Acquisition IDs to select. Default is ``"all"``.

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


def interpolate_channels(n_ch: np.ndarray, x_ch: np.ndarray, y_ch: np.ndarray,
                         z_ch: np.ndarray, system: str = "decimal", err: float = None,
                         spacing: int | float = None)-> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:

	"""Do a linear interpolation between sections of georeferenced channels to
	georeference the non-located channels. Inputs are the georeferenced channels
	in ascending order.

	Parameters
	----------
	n_ch : np.ndarray
		1D array of channel number.
	x_ch : np.ndarray
		1D array of X (longitude) coordinates of the channels specified in ``"n_ch"``.
	y_ch : np.ndarray
		1D array of Y (latitude) coordinates of the channels specified in ``"n_ch"``.
	z_ch : np.ndarray
		1D array of Z (depth - meters) coordinates of the channels specified in ``"n_ch"``.
	system : str
		Defines the output coordinate system for X and Y. It can be ``'decimal'``
		for decimal degrees or ``'utm'`` for Universal Transverse Mercator.
	err: float
		The maximum accepted error of gauge length from interpolation in decimals 0.0 = 0% and 1.0 = 100%.
	spacing : int, float
		Real channel spacing value in meters.

	Returns
	-------
	new_ch : np.ndarray
		-
	new_x : np.ndarray
		X coordinates.
	new_y : np.ndarray
		Y coordinates.
	new_z : np.ndarray
		Z coordinates.

	"""

	new_ch, new_x, new_y, new_z =  [], [], [], []
	for i in range(n_ch.size-1):
		ch_start, ch_end = n_ch[i], n_ch[i+1]

		# deltas or differences between extremes.
		dch =  n_ch[i+1] - n_ch[i]
		dx = x_ch[i+1] - x_ch[i]
		dy = y_ch[i+1] - y_ch[i]
		dz = z_ch[i+1] - z_ch[i]

		# channels = np.linspace(ch_start, ch_end, ch_end - ch_start + 1)
		N = np.linspace(0, 1, ch_end - ch_start + 1) # number of channels in between, counting extremes.
		CH = n_ch[i] + dch * N
		X = x_ch[i] + dx * N
		Y = y_ch[i] + dy * N
		Z = z_ch[i] + dz * N

		ii = 0 if i == 0 else 1 # index for start aappending new georefereced channels.

		new_ch += list(CH[ii:].astype(int))
		new_x += list(X[ii:])
		new_y += list(Y[ii:])
		new_z += list(Z[ii:])

		if err != None:
			if system == "decimal":
				r1 = deg2m(y_ch[i+1], x_ch[i+1], y_ch[i], x_ch[i], implementation="python")[0] # lat/lon distance in meters.
				r2 = np.sqrt(r1**2 + dz**2) # total calculated distance between known locations in 3D.
			elif system == "utm":
				r2 = np.sqrt(dx**2 + dy**2 + dz**2) # total calculated distance between known locations in 3D with UTM.

			calc_spacing = (1/(ch_end - ch_start)) * r2 # calculated channel spacing from georeferencing.
			dev = (calc_spacing - spacing) / spacing # calculated error between new channel spacing and real.
			print(f"Fiber section {i+1} -> channel spacing of {calc_spacing} m. Original spacing is {spacing} m.")

			if np.abs(dev) > err:
				print(f"\n⚠️ Error of calculated channels spacing in fiber section "
					f"{i+1} is of {dev*100}%.\nCheck control points.")
				return None # Calculation is interrupted for not fulfilling standards.

	return np.array(new_ch).astype(int), np.array(new_x), np.array(new_y), np.array(new_z)