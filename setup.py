"""
INSTALL & SETUP FOR THE FOBENCH PYTHON PACKAGE.

Created on 28-01-2024 12:07:17
Last modification on 2023-09-14 14:51:00

:author:
	- Sergio Diaz (sergioad@gfz-potsdam.de)
:contributors:
	- Jonas Pätzel (jonas.patzel@ulb.be)
	- Christopher Wollin (wollin@gfz-potsdam.de)
:license:

"""

# Necessary packages
from setuptools import setup, find_packages


setup(
    name ='fobench',
    author='Sergio Diaz-Meza, Jonas Pätzel',
    author_email='sergioad@gfz.de, jonas.patzel@ulb.be',
    description='A toolbox for basic signal processing of fibre optic sensing data, and data/file management.',
    version = '0.0.31',
    packages = find_packages(),
    install_requires = [
        # Here goes the dependencies !!
        'numpy<2.0.0', # we can check later why is not working with upper versions of numpy. Conflict with other packages. Which one?
        'h5py',
        'nptdms',
        'pandas',
        'matplotlib',
        'scipy',
        'obspy',
        'pyrocko',
        'pyqtgraph',
        'PyQt5',
        'tqdm',
        'zarr',
        'rich' # dependency of zarr
    ],
    extras_require={
        # Here goes the optional dependencies, usually for paralle runinng.
        "mpi":["mpi4py"],
        "dask":["dask[distributed]"],
        "hpc":["mpi4py","dask[distributed]"]
    }
)