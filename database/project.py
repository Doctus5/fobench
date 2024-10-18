"""
Class "Project" for visualizing and handling data within a project.
a Project is understood as a field campaing in an specific location where several Datasets are collected
from deployments.
So far it recieves TDMS format (Silixa) and H5 format (Febus).

Created on 2022-08-19 12:07:17
Last modification on 2024-07-08 15:59:00

:author:
	- Sergio Diaz-Meza (sergioad@gfz-potsdam.de)
:contributors:
	- Jonas Pätzel (jonas.patzel@ulb.be)
	- Christopher Wollin (wollin@gfz-potsdam.de)
:license:

"""

# Necessary packages to import

# Fobench classes

# Inner functions



class Project(object):
	'''
	IMPORTANT INFO: Most of the methods perform changes within the class permanently. Therefore is usefull to make a copy of the class
	with the method copy() before performing any processing or changes.
	'''
	
	#Creates the basic variables of the DAS object with its characteristics
	def __init__(self, metadata_file=None):
		'''
		Co-authors: --
		Description: 
			Initializes a "Project" Class which keeps track of the organization of the single/multiple data collected.
		:Params:
			- metadata_file(type:json): JSON object which represents the metadata of the project. If this file is non existent
			then Project must be initialized without it and a metadata must be created from scratch (see Tutorial).
			- project_type(type:String): type of project to handle based on 
		:Return:
			- NA.  
		'''

		# In case class is initialized first time or empty (no metadata).
		self.units = [] # list of units/interrogators used. Each one is a Unit class that contains Datasets class.

		# Private attributes
		