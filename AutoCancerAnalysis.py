#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: azams and Jonathan Robinson
"""

#%%
import Omics.OmicsData as OD
import os
import warnings

import numpy as np
import pandas as pd

#from sklearn.cross_validation import cross_val_score
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
#from sklearn.linear_model import RidgeClassifier, ElasticNet, ElasticNetCV #, Ridge, Lasso, LinearRegression
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier, AdaBoostClassifier

from xgboost import XGBClassifier
from sklearn import svm
#from sklearn.feature_selection import f_regression #, RFE, VarianceThreshold
#from sklearn.model_selection import KFold, StratifiedKFold, cross_val_score

RS = 20170628

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
ClassVar = 'mutTP53'

# Select which levels of the class variable to keep.
# VarLevelsToKeep = ['Low','Hypermutant']
#VarLevelsToKeep = ['Solid Tissue Normal', 'Primary solid Tumor']
VarLevelsToKeep = ['FALSE', 'TRUE']
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
med_tpm_threshold = 1



#%%
def filterSamplesFromData(dfCancerType, ClassVar, VarLevelsToKeep):
    """
    Remove NaNs and "not reported" values from dataset.
    In addition, if ClassVar is not "CancerStatus", then only keep "Primary 
    solid Tumor" samples in the dataset.
    """
    
    totalsamples = dfCancerType.shape[0]

    dfCancerType = OD.dropNaNs(dfCancerType, ClassVar)    

    if totalsamples > dfCancerType.shape[0]:
        print('Number of samples in the dataset after removing missing values: {0}' \
              .format(dfCancerType.shape[0]))
    
    dfAnalysis = dfCancerType.copy()
    
    ClassVarLevelsFreqTab, ClassVarLevelsSorted = OD.returnVarLevelsSorted(dfAnalysis, ClassVar)
    totalsamples = dfAnalysis.shape[0]
    print('Variable for analysis: ' + '\033[1m{:10s}\033[0m'.format(ClassVar))
    print('Total samples: ' + '\033[1m{:d}\033[0m\n'.format(totalsamples))
    print(ClassVarLevelsFreqTab)
    
    # Keep samples related to Tumor cells only if CancerStatus is not the ClassVar
    if ClassVar != 'CancerStatus':
        toKeep = ['Primary solid Tumor']
#            print('\nKeeping samples concerning "Primary solid Tumor" only.')
        dfAnalysis = OD.FilterLevels(dfAnalysis, 'CancerStatus', toKeep, printStats='no')
    
    # print updated stats if ClassVar was not CancerStatus
    if totalsamples > dfAnalysis.shape[0]:
#            print('Updated, number of samples in the dataset:' + '\033[1m{:d}\033[0m'.format(dfAnalysis.shape[0])) 
        print('\nRemoved {0} samples where CancerStatus was not "Primary solid Tumor".'.format(totalsamples - dfAnalysis.shape[0]))
        ClassVarLevelsFreqTab, ClassVarLevelsSorted = OD.returnVarLevelsSorted(dfAnalysis,ClassVar)
        ClassVarLevelsFreqTab
        
    # sometimes ClassVar is 'not reported' for some samples. We need to remove those as well.
    # and print the updated stats and also update the dataset.
    if 'not reported' in ClassVarLevelsSorted:
        notReported = sum(ClassVarLevelsFreqTab[ClassVarLevelsFreqTab[ClassVar] == 'not reported']['Frequency'])
        print('\nRemoved {0} samples where "{1}" is "not reported".'.format(notReported, ClassVar))
        dfAnalysis.drop(dfAnalysis.index[dfAnalysis[ClassVar] == 'not reported'], inplace= True)
        print('Now, there are '
              + '\033[1m'
              + str(dfAnalysis.shape[0])
              + '\033[0m'
              + ' samples in the dataset.')#.format(dfAnalysis.shape[0]))
        ClassVarLevelsFreqTab, ClassVarLevelsSorted = OD.returnVarLevelsSorted(dfAnalysis, ClassVar)
        ClassVarLevelsFreqTab
    
    # Keep samples only for the values in VarLevelsToKeep while samples corresponding to the rest are filtered out.
    dfAnalysis_fl = OD.FilterLevels(dfAnalysis, ClassVar, VarLevelsToKeep, printStats='no')
    
    ClassVarLevelsFreqTab, ClassVarLevelsSorted = OD.returnVarLevelsSorted(dfAnalysis_fl, ClassVar)
    print(ClassVarLevelsFreqTab)
    
    dfAnalysis_fl = OD.prepareDF(dfAnalysis_fl, ClassVar)

    return dfAnalysis_fl, ClassVarLevelsFreqTab



#%%
def filterGenesFromData(dfAnalysis_fl, CancerType, ClassVar, dimReduction, med_tpm_threshold):
    """
    Remove genes from dataset according to specified parameters.
    """
    
    if dimReduction: # step into dim reduction if dimReduction is set to True
        signifDEgenes = pd.read_excel('PSN_genes_signifDE.xlsx')
        signifDECancerTypes = signifDEgenes.columns[3:].tolist()
        signifDECancerTypes = [s.split('_')[1] for s in signifDECancerTypes]
        if CancerType.split('-')[1] in signifDECancerTypes: # Make sure that the selected cancer type exists
            if dimRedMethod == 'numSigCancers':
                signifDEgenes = signifDEgenes.loc[signifDEgenes['num sig cancers']>=numSigCancers, 'gene name'].tolist()
            elif dimRedMethod == 'signifDEgenes':
                signifDEgenes = signifDEgenes.loc[signifDEgenes['Padj_' + CancerType.split('-')[1]]<=0.01, 'gene name'].tolist()
            signifDEgenes.insert(0,ClassVar)
            dfAnalysis_fl = dfAnalysis_fl[signifDEgenes]
            print('Size of the dataframe after filtering signifDEgenes: {0}'.format(dfAnalysis_fl.shape))
            dfAnalysis_fl_cd = dfAnalysis_fl
        else:
            print('Dim reduction cannot be performed because the cancer type' \
                  '"{0}" does not have paired samples.' \
                  .format(CancerType))
            
    elif med_tpm_threshold != 'none': # remove low-TPM genes if specified, and dim reduction is not requested
        # Look at the list low_tpm_genes, these are the genes which will be removed.
        data_stats, low_tpm_genes = OD.GeneExpression(dfAnalysis_fl,med_tpm_threshold)
        print('\n********************************************************************')
        if type(med_tpm_threshold) == 'str':
            if med_tpm_threshold == 'zero':
                print('The following {0} genes are removed because all their' \
                      'TPM values in the set are zero:' \
                      .format(len(low_tpm_genes)))
            else:
                print('The following {0} genes are removed because their' \
                      'median TPM values lie in the lower {1} percentile of' \
                      'the entire set:' \
                      .format(len(low_tpm_genes),med_tpm_threshold[0:-1]))
        else:
            print('The following {0} genes are removed because their median' \
                  'TPM values are less than {1}:' \
                  .format(len(low_tpm_genes), med_tpm_threshold))
        print(low_tpm_genes)
        
        # Remove low-TPM genes
        dfAnalysis_fl_cd = OD.CleanData(dfAnalysis_fl, med_tpm_threshold)
        print('\nSize of the dataframe after filtering low-TPM genes: {0}' \
              .format(dfAnalysis_fl_cd.shape))
        
    else:
        # Don't remove any genes
        print('No genes were removed from the dataset.')
        dfAnalysis_fl_cd = dfAnalysis_fl
    
    return dfAnalysis_fl_cd



#%%
def performGeneRanking(dfAnalysis_fl_cd, ClassVar, varLevelsToKeep):
    """
    Fit classification models, rank genes (features) based on feature 
    importance scores, and perform a cross-fold validation analysis to assess
    the accuracy and ROC AUC of each model.
    """
    
    # Perform label encoding for the ClassVar and log-transform data
    dfAnalysis_fl_cd = OD.mapClassVar(dfAnalysis_fl_cd, ClassVar, VarLevelsToKeep)
    X, y = OD.fitLogTransform(dfAnalysis_fl_cd, logTransOffset)    
    
    print('Performing ranking of the genes...\n')
    
    geneNames = dfAnalysis_fl_cd.columns[1:].tolist()
    ranks = {}
    
    
    # for random forest methods, use floor(sqrt(numfeats)) as the number of estimators
    num_est = int(X.shape[1]**0.5)
        
    
    extc = ExtraTreesClassifier(n_estimators=num_est, random_state=RS)
    extc.fit(X,y)
    ranks['ExtraTreesClassifier'] = OD.Ranks2Dict(extc.feature_importances_, geneNames)
    print('- ExtraTreesClassifier complete.')
    
    rfc = RandomForestClassifier(n_estimators=num_est, random_state=RS)
    rfc.fit(X,y)
    ranks['RandomForestClassifier'] = OD.Ranks2Dict(rfc.feature_importances_, geneNames)
    print('- RandomForestClassifier complete.')
        
    AdabCLF = AdaBoostClassifier(n_estimators=num_est)
    AdabCLF.fit(X,y)
    ranks['AdaBoostClassifier'] = OD.Ranks2Dict(AdabCLF.feature_importances_, geneNames)
    print('- AdaBoostClassifier complete.')
        
    xgb = XGBClassifier()
    xgb.fit(X, y)
    ranks['XGBClassifier'] = OD.Ranks2Dict(xgb.feature_importances_, geneNames)
    print('- XGBClassifier complete.')
        
    lda =  LinearDiscriminantAnalysis(solver='eigen', shrinkage='auto')
    lda.fit(X, y)
    ranks['LinearDiscriminantAnalysis'] = OD.Ranks2Dict(np.abs(lda.coef_[0]), geneNames)
    print('- LinearDiscriminantAnalysis complete.')
        
    svmSVC = svm.SVC(kernel='linear')
    svmSVC.fit(X,y)
    ranks['SVC'] = OD.Ranks2Dict(np.abs(svmSVC.coef_[0]), geneNames)
    print('- SVC complete.')
    
    # Run a logistic regression using Elastic Net regularization. Run a CV
    # analysis to determine optimal values of the l1_ratio and C parameters.
    logReg = LogisticRegressionCV(cv=10, penalty='elasticnet', scoring='roc_auc',
                                  solver='saga', random_state=RS, n_jobs=-1,
                                  l1_ratios=[0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # ignore the convergence warnings
        logReg.fit(X, y)  # Note: this is quite slow, ~15 min run time
    ranks['LogisticRegression'] = OD.Ranks2Dict(np.abs(logReg.coef_[0]), geneNames)
    print('- LogisticRegression complete.')
    
    # Calculation of individual gene performance (very slow!)
    
#        # Computing genes ranking based on their individual performance using the SVM model.
#        indGeneScores = []
#        for i in range(X.shape[1]):
#             score = cross_val_score(svmSVC, X[:, i:i+1], y, scoring='accuracy',
#                                     cv=10, n_jobs=-1)#StratifiedKFold(n_splits=10, shuffle=True))#ShuffleSplit(len(X), 3, .3))
#             indGeneScores.append((round(np.mean(score), 10)))#, geneNames[i]))
#        ranks['SVMlinearindvGenes'] = dict(zip(geneNames, indGeneScores ))
#        print('- SVMlinearindvGenes complete.')
#    
#        # Computing genes ranking based on their individual performance using the LDA model.
#        indGeneScores = []
#        for i in range(X.shape[1]):
#             score = cross_val_score(lda, X[:, i:i+1], y, scoring='accuracy',
#                                     cv=10, n_jobs=-1)#StratifiedKFold(n_splits=10, shuffle=True))#ShuffleSplit(len(X), 3, .3))
#             indGeneScores.append((round(np.mean(score), 10)))#, geneNames[i]))
#        ranks['LDAindvGenes'] = dict(zip(geneNames, indGeneScores ))
#        print('LDAindvGenes complete.')


    # calculate average rank for each gene    
    r = {}
    for name in geneNames:
        r[name] = round(np.mean([ranks[method][name] for method in ranks.keys()]), 10)
    ranks['Average'] = r
    
    
    # Recursive feature elimination (RFE) methods

#    rfeSVM = RFE(svm.SVC(kernel='linear'), n_features_to_select=1)
#    rfeSVM.fit(X,y)
#    ranks['rfeSVM'] = dict(zip(geneNames, rfeSVM.ranking_ ))
#    print('- rfeSVM complete.')
#    methods.append('SVMrfe')
#    
#    rfeET = RFE(extc, n_features_to_select=1)
#    rfeET.fit(X,y)
#    ranks['rfeExtraTree'] = dict(zip(geneNames, rfeET.ranking_ ))
#    methods.append('rfeExtraTree')
#    
#    rfeRFC = RFE(RandomForestClassifier(n_estimators=200, random_state=RS), n_features_to_select=1)
#    rfeRFC.fit(X,y)
#    ranks['rfeRFC'] = dict(zip(geneNames, rfeRFC.ranking_ ))#OD.Ranks2Dict(np.array(rfeRFC.ranking_, dtype=float), geneNames)#, order=1)
#    print('- rfeRFC complete.')    

    # organize and sort ranks
    dfRanks = pd.DataFrame.from_dict(ranks)
    dfRanks.reset_index(inplace=True)
    dfRanks.rename(columns={'index':'GeneNames'}, inplace=True)
    dfRanks.sort_values(by='Average', inplace=True, ascending=False)
    
    print('\nDone!\n')
    print('\n********************************************************************')
    
    
    # Run model cross-validation and determine accuracy and ROC AUC
    models = [ExtraTreesClassifier(n_estimators=num_est, random_state=RS), # 0
              RandomForestClassifier(n_estimators=num_est, random_state=RS), # 1
              AdaBoostClassifier(n_estimators=num_est), # 2
              XGBClassifier(), # 3
              LinearDiscriminantAnalysis(), # 4
              svm.SVC(kernel='linear'), # 5
              LogisticRegression(penalty='elasticnet', C=logReg.C_[0],
                                 l1_ratio=logReg.l1_ratio_[0], solver='saga',
                                 random_state=RS)  # 6
             ]

    CV = 'Validation: SKF'
    shuffle = True
    scoring = 'accuracy'
    folds = 10

    print('Performing models CV analysis using accuracy...\n')
    dfCVscores_accuracy = OD.CVScorer(models, CV, X, y, scoring, shuffle, folds)
    print('\nDone!\n')
    
    if len(VarLevelsToKeep) == 2:
        scoring = 'roc_auc'
        print('Performing models CV analysis using area under the ROC curve...\n')
        dfCVscores_ROC = OD.CVScorer(models, CV, X, y, scoring, shuffle, folds)
        print('\nDone!\n')
    else:
        print('Skipping CV analysis using area under the ROC curve. ' \
              'This is possible for binary problems only.')
    
    return dfRanks, dfCVscores_accuracy, dfCVscores_ROC



#%%
def writeResultsToFile(dfRanks, dfCVscores_accuracy, dfCVscores_ROC, CancerType, VarLevelsToKeep):
    """
    Export gene ranks and model scores (accuracy and ROC AUC) to .csv files
    """

    parent_dir_name = 'results/'
    
    print('Writing dataset, genees ranking and CV analysis results to a ' \
          'directory named "{0}"'.format(CancerType))
    
    os.makedirs(parent_dir_name + CancerType , exist_ok=True)

    if ClassVar in ['TumorStage', 'TumorStageMerged', 'TumorStageBinary']:
        file_name_piece = '_'.join(['TumorStage'] + VarLevelsToKeep)
        file_name_piece = file_name_piece.replace(' ','')
    else:
        file_name_piece = ClassVar

    dfRanks.to_csv(parent_dir_name + CancerType + '/' + CancerType + '_' \
                   + file_name_piece + '_GenesRanking.csv', index=False)    
    dfCVscores_accuracy.to_csv(parent_dir_name + CancerType + '/' + CancerType \
                               + '_' + file_name_piece + '_CVscores_Accuracy.csv', index=False)
    if len(VarLevelsToKeep) == 2:
        dfCVscores_ROC.to_csv(parent_dir_name + CancerType + '/' + CancerType \
                              + '_' + file_name_piece + '_CVscores_ROCAUC.csv', index=False)
    
    print('\nDone!\n')



#%%
    
    
# Loop through each cancer type, performing the analysis on each type
for CancerType in allCancerTypes:
    
#    CancerType = 'COAD'
    
    CancerDataStore = pd.HDFStore('data/CancerDataStore_psn.h5')
    dfCancerType = CancerDataStore.get(CancerType)
    CancerDataStore.close()

    print('Cancer Type: ' + '\033[1m{:10s}\033[0m'.format(CancerType))
    
    colnames = list(dfCancerType)  # get list of all class variables available

    if ClassVar == 'Mutations':
        all_mutClassVars = [s for s in colnames if 'mut' == s[0:3]]  # extract mutation variables
        for ClassVar in all_mutClassVars:                        
            if (CancerType) in os.listdir('results'):
                if any([True for x in os.listdir(os.getcwd() + '/results/' + CancerType) if ClassVar + '_GenesRanking' in x]):
                    print('Already analyzed; skipping.')
                    continue
            # filter samples from data
            dfAnalysis_fl, ClassVarLevelsFreqTab = filterSamplesFromData(dfCancerType, ClassVar, VarLevelsToKeep)
            
            # check if there are at least 10 samples in each class, and at least 2 classes
            if ((ClassVarLevelsFreqTab['Frequency'].min() < 10) or (ClassVarLevelsFreqTab.shape[0] < 2)):
                print('Insufficient samples to perform analysis; skipping.')
                continue
            
            # filter genes from data
            dfAnalysis_fl_cd = filterGenesFromData(dfAnalysis_fl, CancerType, ClassVar, dimReduction, med_tpm_threshold)
            
            # fit models, rank genes, and perform cross-validation
            dfRanks, dfCVscores_accuracy, dfCVscores_ROC = performGeneRanking(dfAnalysis_fl_cd, ClassVar, VarLevelsToKeep)
            
            # write results to file
            writeResultsToFile(dfRanks, dfCVscores_accuracy, dfCVscores_ROC, CancerType, VarLevelsToKeep)
            
    else: 
        
        if (CancerType) in os.listdir('results'):
                if any([True for x in os.listdir(os.getcwd() + '/results/' + CancerType) if ClassVar + '_GenesRanking' in x]):
                    print('Already analyzed; skipping.')
                    continue
        
        # filter samples from data
        dfAnalysis_fl, ClassVarLevelsFreqTab = filterSamplesFromData(dfCancerType, ClassVar, VarLevelsToKeep)
        
        
        # check if there are at least 10 samples in each class, and at least 2 classes
        if ((ClassVarLevelsFreqTab['Frequency'].min() < 10) or (ClassVarLevelsFreqTab.shape[0] < 2)):
            print('Insufficient samples to perform analysis; skipping.')
            continue
        
        # filter genes from data
        dfAnalysis_fl_cd = filterGenesFromData(dfAnalysis_fl, CancerType, ClassVar, dimReduction, med_tpm_threshold)
        
        # fit models, rank genes, and perform cross-validation
        dfRanks, dfCVscores_accuracy, dfCVscores_ROC = performGeneRanking(dfAnalysis_fl_cd, ClassVar, VarLevelsToKeep)
        
        # write results to file
        writeResultsToFile(dfRanks, dfCVscores_accuracy, dfCVscores_ROC, CancerType, VarLevelsToKeep)
        
        
    
    
