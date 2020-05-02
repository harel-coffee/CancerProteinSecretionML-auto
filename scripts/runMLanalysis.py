#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: azams and Jonathan Robinson
"""

#%%
import omicsAnalysisFunctions as OF
import os
import pandas as pd

RS = 20170628
proj_dir = os.path.dirname(os.getcwd())

#%%
#==============================================================================
# Analysis Parameters
#==============================================================================

# all available cancer types
allCancerTypes = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC',
                  'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LGG', 
                  'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG',
                  'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA',
                  'THYM', 'UCEC', 'UCS', 'UVM']

# ClassVar options: 'CancerStatus','TumorStage','TumorStageMerged','TumorStageBinary',
#                   'OverallSurvival','Race','Gender','Barcode','Mutations',
#                   'HyperMut','HyperMutBinary'
ClassVar = 'CancerStatus'

# Select which levels of the class variable to keep.
# VarLevelsToKeep = ['Low','Hypermutant']
VarLevelsToKeep = ['Solid Tissue Normal', 'Primary solid Tumor']
# VarLevelsToKeep = ['FALSE', 'TRUE']
#VarLevelsToKeep = ['stage i','stage iii']
# VarLevelsToKeep = ['stage i-iii','stage iv']

# specify offset to add to TPM values before log-transforming (to handle zeros)
logTransOffset = 1  # transformed TPM = log(TPM + offset)

# dimensionality reduction options
dimReduction = False  # True or False

# if dimReduction is False, the following two variables are NOT used:
dimRedMethod = 'numSigCancers' # 'signifDEgenes' or 'numSigCancers'
numSigCancers = 10  # Number of cancers in which gene must be significant 

# if dimReduction is False, there is option to remove low-TPM genes by
# specifying the med_tpm_threshold parameter:
#  'none' - don't remove any genes.
#  'zero' - remove genes with all zeros.
#  'X%' - where X is a number from 0 to 100, removes genes with median TPM
#         in the bottom X-percentile.
#   X - where X is a number, removes genes with median TPM below X
med_tpm_threshold = 0.1


#%%
        
# Loop through each cancer type, performing the analysis on each type
for CancerType in allCancerTypes:
    
#    CancerType = 'COAD'
    
    CancerDataStore = pd.HDFStore(proj_dir + '/data/CancerDataStore_psp.h5')
    dfCancerType = CancerDataStore.get(CancerType)
    CancerDataStore.close()

    print('Cancer Type: ' + '\033[1m{:10s}\033[0m'.format(CancerType))
    
    colnames = list(dfCancerType)  # get list of all class variables available

    if ClassVar == 'Mutations':
        all_mutClassVars = [s for s in colnames if 'mut' == s[0:3]]  # extract mutation variables
        for mutClassVar in all_mutClassVars:                        
            if (CancerType) in os.listdir('results'):
                if any([True for x in os.listdir(os.getcwd() + '/results/' + CancerType) if mutClassVar + '_GenesRanking' in x]):
                    print('Already analyzed; skipping.')
                    continue
            # filter samples from data
            dfAnalysis_fl, ClassVarLevelsFreqTab = OF.filterSamplesFromData(dfCancerType, mutClassVar, VarLevelsToKeep)
            
            # check if there are at least 10 samples in each class, and at least 2 classes
            if ((ClassVarLevelsFreqTab['Frequency'].min() < 10) or (ClassVarLevelsFreqTab.shape[0] < 2)):
                print('Insufficient samples to perform analysis; skipping.')
                continue
            
            # filter genes from data
            dfAnalysis_fl_cd = OF.filterGenesFromData(dfAnalysis_fl, CancerType, mutClassVar, dimReduction, med_tpm_threshold)
            
            # fit models, rank genes, and perform cross-validation
            dfRanks, dfCVscores_accuracy, dfCVscores_ROC = OF.performGeneRanking(dfAnalysis_fl_cd, mutClassVar, VarLevelsToKeep)
            
            # write results to file
            OF.writeResultsToFile(dfRanks, dfCVscores_accuracy, dfCVscores_ROC, CancerType, mutClassVar, VarLevelsToKeep)
            
    else: 
        
        if (CancerType) in os.listdir(proj_dir + '/results'):
                if any([True for x in os.listdir(proj_dir + '/results/' + CancerType) if ClassVar + '_GenesRanking' in x]):
                    print('Already analyzed; skipping.')
                    continue
        
        # filter samples from data
        dfAnalysis_fl, ClassVarLevelsFreqTab = OF.filterSamplesFromData(dfCancerType, ClassVar, VarLevelsToKeep)
        
        
        # check if there are at least 10 samples in each class, and at least 2 classes
        if ((ClassVarLevelsFreqTab['Frequency'].min() < 10) or (ClassVarLevelsFreqTab.shape[0] < 2)):
            print('Insufficient samples to perform analysis; skipping.')
            continue
        
        # filter genes from data
        dfAnalysis_fl_cd = OF.filterGenesFromData(dfAnalysis_fl, CancerType, ClassVar, dimReduction, med_tpm_threshold)
        
        # fit models, rank genes, and perform cross-validation
        dfRanks, dfCVscores_accuracy, dfCVscores_ROC = OF.performGeneRanking(dfAnalysis_fl_cd, ClassVar, VarLevelsToKeep, logTransOffset, RS)
        
        # write results to file
        resultsPath = proj_dir + '/results/'
        OF.writeResultsToFile(dfRanks, dfCVscores_accuracy, dfCVscores_ROC, CancerType, ClassVar, VarLevelsToKeep, resultsPath)
        
        
    
    
