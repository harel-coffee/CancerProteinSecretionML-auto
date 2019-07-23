#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 10:43:47 2017

@author: azams

We want to analyze how well different cancer types are separated from each other.
We filter out samples where 'CancerStatus' = 'Primary solid Tumor' and ClassVar = 'Project'.
Then we chose which CancerTypes to compare against each other for the analysis purposes. 

"""

import Omics.OmicsData as OD
import numpy as np
from sklearn.manifold import TSNE
from sklearn.ensemble import IsolationForest

#import matplotlib.pyplot as plt
import pandas as pd
#import itertools as itr


#%% Read the PSN cancer data
dfAllOD = OD.ReadOMICSdataCSV("allcancerdata")

dfdna_pst, CancerTypesSorted = OD.dfCancerTypesOrdered(dfAllOD)
ClassVar = 'Project'
#%%
# Once the data is loaded after running the above lines, then only run the
# following code to run it for different cancer types.

## Belows are the three options to select the cancer types to compare against
## each other.
#CancerTypes = CancerTypesSorted[1:5]
CancerTypes = [CancerTypesSorted[i] for i in [3,5]]
#CancerTypes = ['TCGA-HNSC' ,'TCGA-LUSC']

dfdna_pst_fl = OD.FilterLevels(dfdna_pst, ClassVar, CancerTypes)
dfdna_pst_fl, ClassVarEncOrder = OD.mapClassVar(dfdna_pst_fl,ClassVar)
dfdna_pst_fl_cd = OD.CleanData(dfdna_pst_fl,2)
    
  
X, y = OD.fitLogTransform(dfdna_pst_fl_cd)
#X, y = OD.fitScalarTransform(dfdna_pst_fl_cd)
# Now let's run the t-SNE algorithm on the dataset.
# Random state.
RS = 20170628
omics_proj = TSNE(random_state=RS).fit_transform(X)

OD.tSNEscatter(omics_proj, y, ClassVarEncOrder, len(CancerTypes))
OD.plotPCA(X, y, 2, CancerTypes, save=False)
if len(CancerTypes) >2:
    OD.plotLDA(X, y, 2, CancerTypes, save=False)
    
#%%
## Here in this cell we take the same selection as from the previous cell
## and remove outliers first and then plot.
    
isoForest = IsolationForest(n_estimators=1000, contamination=0.40)

for i in np.unique(y):
    Xi = X[y == i]
    Yi = y[y == i]
    isoForest.fit(Xi)
    outlierIsoFor = isoForest.predict(Xi)
    pd.DataFrame(outlierIsoFor)[0].value_counts()
#outlierTF = outlier == 1
#dftesting_clean =dftesting[outlierTF]
    Xicleaned = Xi[outlierIsoFor == 1]
    Yicleaned = Yi[outlierIsoFor == 1]
    if i==0:
        Xcleaned = Xicleaned
        Ycleaned = Yicleaned
    else:
        
        Xcleaned = np.append(Xcleaned, Xicleaned, axis=0)
        Ycleaned = np.append(Ycleaned, Yicleaned, axis=0)


omics_projcleaned = TSNE(random_state=RS).fit_transform(Xcleaned)

OD.tSNEscatter(omics_projcleaned, Ycleaned, ClassVarEncOrder, len(CancerTypes))
OD.plotPCA(Xcleaned, Ycleaned, 2, CancerTypes, save=False)
if len(CancerTypes) >2:
    OD.plotLDA(Xcleaned, Ycleaned, 2, CancerTypes, save=False)
#OD.CancerTypesDiscAnalysis(dfAllOD, CancerTypes, nComp = 2, save=False)
#%%
"""
#Following is a count of the number of cases corresponding to each cancer type.
Specify the index (indices) of the cancer types (Project-IDs) which should be included
in the comparison. 0:10 means Project-IDs corresponding to indices 0-9 will be included.
0:1 means Project-ID corresponding to index 0 will be included.
6:7 means Project-ID corresponding to index 6 will be included.
 
Project-ID  Frequency  Index
TCGA-BRCA    1101       0
TCGA-UCEC     551       1
TCGA-KIRC     538       2
TCGA-LUAD     533       3
TCGA-LGG      510       4
TCGA-LUSC     502       5
TCGA-THCA     502       6
TCGA-HNSC     500       7
TCGA-PRAD     498       8
TCGA-COAD     476       9
TCGA-BLCA     414       10
TCGA-STAD     375       11 
TCGA-OV       374       12
TCGA-LIHC     371       13
TCGA-CESC     304       14
TCGA-KIRP     288       15
TCGA-SARC     259       16
TCGA-PCPG     178       17
TCGA-PAAD     177       18
TCGA-READ     165       19
TCGA-ESCA     161       20
TCGA-GBM      155       21
TCGA-TGCT     134       22
TCGA-THYM     119       23
TCGA-SKCM     103       24
TCGA-MESO      86       25
TCGA-UVM       80       26
TCGA-ACC       79       27
TCGA-KICH      65       28
TCGA-UCS       56       29
TCGA-DLBC      48       30
TCGA-CHOL      36       31

Following is the list of all cancer types. Choose the ones that you want to compare and pass to the above function.
['TCGA-ACC', 'TCGA-BLCA', 'TCGA-BRCA', 'TCGA-CESC', 'TCGA-CHOL', 'TCGA-COAD', 'TCGA-DLBC', 
'TCGA-ESCA', 'TCGA-GBM', 'TCGA-HNSC', 'TCGA-KICH', 'TCGA-KIRC', 'TCGA-KIRP', 'TCGA-LGG', 
'TCGA-LIHC', 'TCGA-LUAD', 'TCGA-LUSC', 'TCGA-MESO', 'TCGA-OV', 'TCGA-PAAD', 'TCGA-PCPG', 
'TCGA-PRAD', 'TCGA-READ', 'TCGA-SARC', 'TCGA-SKCM', 'TCGA-STAD', 'TCGA-TGCT', 'TCGA-THCA', 
'TCGA-THYM', 'TCGA-UCEC', 'TCGA-UCS', 'TCGA-UVM']

"""