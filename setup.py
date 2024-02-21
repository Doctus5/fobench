"""
INSTALL & SETUP FOR THE FOBENCH PYTHON PACKAGE.

Created on 28-01-2024 12:07:17
Last modification on 2023-09-14 14:51:00

:author:
	- Sergio Diaz (sergioad@gfz-potsdam.de)
:contributors:
	- Christopher Wollin (wollin@gfz-potsdam.de)
:license:

"""

# Necessary packages
from setuptools import setup


setup(
    name ='fobench',
    version = '0.0.1',
    packages = ['fobench'],
    install_requires = [
        # Here goes the dependencies !!
        'h5py',
        'nptdms',
        'numpy',
        'matplotlib',
        'scipy',
        'obspy',
        'pyrocko'
    ]
)