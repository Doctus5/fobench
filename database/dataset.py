"""
Class "Dataset" for visualizing, processing and handling a dataset.
a Dataset is understood as group of fiber optic sensing data ('ddss/das','dss','dts') files part of a campaign.
Can be continuous data in many files, and may have discontinuities.

Created on 2024-07-08 16:26:00
Last modification on 2024-07-08 16:26:00

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



class Dataset(object):
    '''
    IMPORTANT INFO: Most of the methods perform changes within the class permanently. Therefore is usefull to make a copy of the class
    with the method copy() before performing any processing or changes.
    '''
	
	#Creates the basic variables of the DAS object with its characteristics
    def __init__(self, filepath, company='silixa', sensing='das', scan=True):
        '''
        Co-authors: --
        Description: 
            Initializes a "Dataset" Class which keeps track of the organization of the single/multiple files collected.
        :Params:
            - filepath(type:String): folder/compelte path where the single/multiple files of data are located.
			- company(type:String): manufacturer or the instrument that generates the data. Currently supporting "silixa" (Default), "febus", and "bam".
			- sensing(type:String): specifies the type of fiber optic sensing technique of the data. Default is 'das'.
        :Return:
            - NA.  
        '''

        # internal attributed.
        self.__filepath__ = filepath # filepath where the data is located (a folder).
        
        self.company = company # company of the manufacturer where the data comes from. Important to know how to read.
        self.sensing = sensing # sensing target of the dataset.
        self.files = 0
        