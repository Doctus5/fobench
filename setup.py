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
from setuptools import setup


setup(
    name ='fobench',
    author='Sergio Diaz-Meza',
    author_email='sergioad@gfz-potsdam.de',
    description='A fiber optic sensing toolbox',
    version = '0.0.21',
    packages = ['fobench'],
    install_requires = [
        # Here goes the dependencies !!
        'numpy<2.0.0', # we can check later why is not working with upper versions of numpy. Conflict with other packages. Which one?
        'h5py',
        'nptdms',
        'pandas',
        'matplotlib',
        'scipy',
        'obspy',
        'pyrocko'
        'pyqtgraph'
    ]
)