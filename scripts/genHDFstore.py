#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2019-08-01

Author: Jonathan Robinson
"""

import sys
import numpy as np
import pandas as pd
import Omics.OmicsData as OD

def prepCancerTypeDict(inFile, outFile):
    """
    This function loads the entire PSN cancer dataset from  'allcancerdata.csv' 
    and returns a disctionaly of dataframes, where each dataframe corresponds 
    to a cancer types. 
    """ 
    
    # Import data from csv to a data frame
    if '.csv' in inFile:
        df = pd.read_csv('data/' + inFile)
    else:
        df = pd.read_csv('data/' + inFile + ".csv")
    # There is a column with name Unnamed: 0. 
    # Dropping it here.
    if 'Unnamed: 0' in df.columns:
        df.drop('Unnamed: 0', axis=1, inplace=True)

    df = df.dropna(subset = ['Project'])
    projects = df['Project'].unique()
    arr = []
    for project in projects:
        arr.append(project)
    arr = np.array(arr)
    
    # Create a dictionary of data frames separated by cancer type
    cancerTypesDic = dict()
    for project in arr:
        ClassVar = 'Project'
        toKeep = [project]
        cancerTypesDic[project]= OD.FilterLevels(df, ClassVar, toKeep, printStats='no')
    
    # write the dictionary to an hdfStore.
    if '.h5' in outFile:
        CancerDataStore = pd.HDFStore('data/' + outFile)
    else:    
        CancerDataStore = pd.HDFStore('data/' + outFile + '.h5')

    for (key, value) in cancerTypesDic.items():
        # keys are names of cancers, e.g., TCGA-BRCA. Using split to ignore the TCGA- part and use
        # the rest as the name. With prefix TCGA-, it is not a valid Python identifier.
        CancerDataStore.put(key.split('-')[1], value)
        print("{0} successfully saved in store!".format(key))
        
    #print(CancerDataStore)
    CancerDataStore.close()
    
    return cancerTypesDic


if len(sys.argv) < 2:
    inFile='allcancerdata'
else:
    inFile=sys.argv[1]

if len(sys.argv) < 3:
    outFile='CancerDataStore'
else:
    outFile=sys.argv[2]

prepCancerTypeDict(inFile, outFile)


