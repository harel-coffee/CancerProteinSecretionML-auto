#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: azams and Jonathan Robinson
"""

#%%
from os import path, listdir, getcwd, mkdir
import pandas as pd

# must be in the project "scripts" directory to import OF
import omicsAnalysisFunctions as OF

RS = 20170628
proj_dir = path.dirname(getcwd())


#%%
#==============================================================================
# Analysis Parameters
#==============================================================================

# ClassVar options: 'mutTP53', 'CancerStatus', 'TumorStageMerged',
#                   'Race', 'Gender', 'Mutations', 'AllStageCombos'
ClassVar = 'CancerStatus'

# Select which levels of the class variable to keep. Not needed for classVar='AllStageCombos'
VarLevelsToKeep = ['Solid Tissue Normal', 'Primary solid Tumor']  # for 'CancerStatus'
# VarLevelsToKeep = ['FALSE', 'TRUE']  # for 'mutTP53' or 'Mutations'
# VarLevelsToKeep = ['stage iv','stage x']  # for 'TumorStageMerged'
# VarLevelsToKeep = ['stage i', 'stage ii', 'stage iii', 'stage iv']  # for 'TumorStageMerged' regression

# specify offset to add to TPM values before log-transforming (to handle zeros)
logTransOffset = 1  # transformed TPM = log(TPM + offset)

# Optional removal of low-TPM genes by specifying the med_tpm_threshold parameter:
#  'none' - don't remove any genes.
#  'zero' - remove genes with all zeros.
#  'X%' - where X is a number from 0 to 100, removes genes with median TPM
#         in the bottom X-percentile.
#   X - where X is a number, removes genes with median TPM below X
med_tpm_threshold = 0.1

# Specify the scoring metric that the model will try to maximize
#   Examples for binary classification: 'accuracy', 'average_precision', 'f1', 'roc_auc'
#   Examples for regression: 'explained_variance', 'neg_mean_squared_error', 'r2'
score_metric = 'roc_auc'


#==============================================================================
# Define some paths and other run settings
#==============================================================================

dataStoreFile = 'CancerDataStore_psp.h5'
output_dir = 'results'
overwrite_results = False
allCancerTypes = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC',
                  'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LGG', 
                  'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG',
                  'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA',
                  'THYM', 'UCEC', 'UCS', 'UVM']

#%%
#==============================================================================
# Main analysis section
#==============================================================================

# create output directory if it does not yet exist
if not path.isdir(path.join(proj_dir, output_dir)):
    mkdir(path.join(proj_dir, output_dir))

# Loop through each cancer type, performing the analysis on each type
for CancerType in allCancerTypes:
    
#    CancerType = 'COAD'
    
    CancerDataStore = pd.HDFStore(proj_dir + '/data/' + dataStoreFile)
    dfCancerType = CancerDataStore.get(CancerType)
    CancerDataStore.close()

    print('Cancer Type: ' + '\033[1m{:10s}\033[0m'.format(CancerType))
    
    colnames = list(dfCancerType)  # get list of all class variables available

    if ClassVar == 'Mutations':
        all_mutClassVars = [s for s in colnames if 'mut' == s[0:3]]  # extract mutation variables
        for mutClassVar in all_mutClassVars:                        
            if (not overwrite_results) and (CancerType) in listdir(proj_dir + '/' + output_dir):
                if any([True for x in listdir(proj_dir + '/' + output_dir + '/' + CancerType) if mutClassVar + '_GenesRanking' in x]):
                    print('Already analyzed; skipping.')
                    continue
            # filter samples from data
            dfAnalysis_fl, ClassVarLevelsFreqTab = OF.filterSamplesFromData(dfCancerType, mutClassVar, VarLevelsToKeep)
            
            # check if there are at least 10 samples in each class, and at least 2 classes
            if ((ClassVarLevelsFreqTab['Frequency'].min() < 10) or (ClassVarLevelsFreqTab.shape[0] < 2)):
                print('Insufficient samples to perform analysis; skipping.')
                continue
            
            # filter genes from data
            dfAnalysis_fl_cd = OF.filterGenesFromData(dfAnalysis_fl, CancerType, mutClassVar, med_tpm_threshold)
            
            # fit models, rank genes, and perform cross-validation
            dfRanks, dfCVscores = OF.performGeneRanking(dfAnalysis_fl_cd, mutClassVar, VarLevelsToKeep, logTransOffset, RS, score_metric)
            
            # write results to file
            resultsPath = proj_dir + '/' + output_dir + '/'
            OF.writeResultsToFile(dfRanks, dfCVscores, CancerType, mutClassVar, VarLevelsToKeep, resultsPath)

    elif ClassVar == 'AllStageCombos':
        all_tumor_combinations = [['stage i', 'stage ii'], ['stage i', 'stage iii'], ['stage i', 'stage iv'], \
                                  ['stage i', 'stage x'], ['stage ii', 'stage iii'], ['stage ii', 'stage iv'], \
                                  ['stage ii', 'stage x'], ['stage iii', 'stage iv'], ['stage iii', 'stage x'], \
                                  ['stage iv', 'stage x']]
        for stage_combo in all_tumor_combinations:
            ClassVar = 'TumorStageMerged'
            VarLevelsToKeep = stage_combo
            if (not overwrite_results) and (CancerType) in listdir(proj_dir + '/' + output_dir):
                file_name_piece = '_'.join(['TumorStage'] + VarLevelsToKeep)
                file_name_piece = file_name_piece.replace(' ','')
                if any([True for x in listdir(proj_dir + '/' + output_dir + '/' + CancerType) if file_name_piece + '_GenesRanking' in x]):
                    print('Already analyzed; skipping.')
                    continue
            
            # filter samples from data
            dfAnalysis_fl, ClassVarLevelsFreqTab = OF.filterSamplesFromData(dfCancerType, ClassVar, VarLevelsToKeep)
            
            # check if there are at least 10 samples in each class, and at least 2 classes
            if ((ClassVarLevelsFreqTab['Frequency'].min() < 10) or (ClassVarLevelsFreqTab.shape[0] < 2)):
                print('Insufficient samples to perform analysis; skipping.')
                continue
            
            # filter genes from data
            dfAnalysis_fl_cd = OF.filterGenesFromData(dfAnalysis_fl, CancerType, ClassVar, med_tpm_threshold)
            
            # fit models, rank genes, and perform cross-validation
            dfRanks, dfCVscores = OF.performGeneRanking(dfAnalysis_fl_cd, ClassVar, VarLevelsToKeep, logTransOffset, RS, score_metric)
            
            # write results to file
            resultsPath = proj_dir + '/' + output_dir + '/'
            OF.writeResultsToFile(dfRanks, dfCVscores, CancerType, ClassVar, VarLevelsToKeep, resultsPath)
        
        # re-assign class variable after looping
        ClassVar = 'AllStageCombos'
        
    else: 
        if (not overwrite_results) and (CancerType) in listdir(proj_dir + '/' + output_dir):
            if any([True for x in listdir(proj_dir + '/' + output_dir + '/' + CancerType) if ClassVar + '_GenesRanking' in x]) or \
            (len(VarLevelsToKeep) > 2 and any([True for x in listdir(proj_dir + '/' + output_dir + '/' + CancerType) if 'regression' in x])):
                print('Already analyzed; skipping.')
                continue
        
        # filter samples from data
        dfAnalysis_fl, ClassVarLevelsFreqTab = OF.filterSamplesFromData(dfCancerType, ClassVar, VarLevelsToKeep)
        
        
        # check if there are at least 10 samples in each class, and that all classes are present
        if ((ClassVarLevelsFreqTab['Frequency'].min() < 10) or (ClassVarLevelsFreqTab.shape[0] < len(VarLevelsToKeep))):
            print('Insufficient samples to perform analysis; skipping.')
            continue
        
        # filter genes from data
        dfAnalysis_fl_cd = OF.filterGenesFromData(dfAnalysis_fl, CancerType, ClassVar, med_tpm_threshold)
        
        # fit models, rank genes, and perform cross-validation
        dfRanks, dfCVscores = OF.performGeneRanking(dfAnalysis_fl_cd, ClassVar, VarLevelsToKeep, logTransOffset, RS, score_metric)
        
        # write results to file
        resultsPath = proj_dir + '/' + output_dir + '/'
        OF.writeResultsToFile(dfRanks, dfCVscores, CancerType, ClassVar, VarLevelsToKeep, resultsPath)
        
        
    
    
