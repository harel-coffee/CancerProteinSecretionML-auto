#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: azams and Jonathan Robinson
"""

#%%
import Omics.OmicsData as OD
import os

import numpy as np
import pandas as pd

#from sklearn.cross_validation import cross_val_score
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.linear_model import RidgeClassifier, RandomizedLasso #, RandomizedLogisticRegression, Ridge, Lasso, LinearRegression
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier, AdaBoostClassifier

from xgboost import XGBClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn import svm
from sklearn.feature_selection import RFE, f_regression #, VarianceThreshold
#from sklearn.model_selection import KFold, StratifiedKFold, cross_val_score

from minepy import MINE
RS = 20170628
#%%
# Run these two lines if there are import errors related to our library: OmicsData
# import importlib
# importlib.reload(OD)

#%%
#==============================================================================
# Analysis Parameters
#==============================================================================

# all available cancer types
allCancerTypes = ['TCGA-ACC', 'TCGA-BLCA', 'TCGA-BRCA', 'TCGA-CESC', 'TCGA-CHOL',
                  'TCGA-COAD', 'TCGA-DLBC', 'TCGA-ESCA', 'TCGA-GBM', 'TCGA-HNSC',
                  'TCGA-KICH', 'TCGA-KIRC', 'TCGA-KIRP', 'TCGA-LGG', 'TCGA-LIHC',
                  'TCGA-LUAD', 'TCGA-LUSC', 'TCGA-MESO', 'TCGA-OV', 'TCGA-PAAD',
                  'TCGA-PCPG', 'TCGA-PRAD', 'TCGA-READ', 'TCGA-SARC', 'TCGA-SKCM',
                  'TCGA-STAD', 'TCGA-TGCT', 'TCGA-THCA', 'TCGA-THYM', 'TCGA-UCEC',
                  'TCGA-UCS', 'TCGA-UVM']

# ClassVar options: 'CancerStatus','TumorStage','TumorStageMerged','TumorStageBinary',
#                   'OverallSurvival','Race','Gender','Barcode','Mutations',
#                   'HyperMut','HyperMutBinary'
ClassVar = 'TumorStageMerged'

# Select which levels of the class variable to keep.
#VarLevelsToKeep = ['Low','Hypermutant']
#VarLevelsToKeep = ['Primary solid Tumor','Solid Tissue Normal']
#VarLevelsToKeep = [True, False]
VarLevelsToKeep = ['stage i','stage iv']

# specify offset to add to TPM values before log-transforming (to handle zeros)
logTransOffset = 1  # transformed TPM = log(TPM + offset)

# dimensionality reduction options
dimReduction = False  # True or False
dimRedMethod = 'numSigCancers' # 'signifDEgenes' or 'numSigCancers'
numSigCancers = 10  # Number of cancers in which gene must be significant 

# if dimReduction is false, there is option to remove low-TPM genes by
# specifying the med_tpm_threshold parameter:
#  'none' - don't remove any genes.
#  'zero' - remove genes with all zeros.
#  'X%' - where X is a number from 0 to 100, removes genes with median TPM
#         in the bottom X-percentile.
#   X - where X is a number, removes genes with median TPM below X
med_tpm_threshold = 1


#%%
# Loop through each cancer type, performing the analysis on each

for CancerType in allCancerTypes:
    
    #CancerType = 'TCGA-COAD'
    CancerDataStore = pd.HDFStore('data/CancerDataStore_psn.h5')
    dfCancerType = CancerDataStore.get(CancerType.split('-')[1])
    CancerDataStore.close()
#    print('Number of samples in the dataset before removing missing values: {0}' \
#          .format(dfCancerType.shape[0])) 
    print('Cancer Type: ' + '\033[1m{:10s}\033[0m'.format(CancerType))
    totalsamples = dfCancerType.shape[0]
    #print('Number of samples in the dataset: {0}'.format(totalsamples)) 
    
    # get list of all class variables available, and narrow to mutation variables
    colnames = list(dfCancerType)
#    all_mutClassVars = [s for s in colnames if 'mut' == s[0:3]]
    
#    for ClassVar in all_mutClassVars:
    #check if the combination of CancerType and ClassVar have already been analyzed
#    if (CancerType) in os.listdir():
#        if any([True for x in os.listdir(os.getcwd() + '/' + CancerType) if ClassVar.split(sep='_')[1] + 'CVscores' in x]):
#            continue
        
    dfCancerType = OD.dropNaNs(dfCancerType,ClassVar)
    
    if totalsamples > dfCancerType.shape[0]:
        print('Number of samples in the dataset after removing missing values: {0}' \
              .format(dfCancerType.shape[0]))
    
    dfAnalysis = dfCancerType.copy()
    #ClassVar = 'CancerStatus'
    ClassVarLevelsFreqTab, ClassVarLevelsSorted = OD.returnVarLevelsSorted(dfAnalysis,ClassVar)
    totalsamples = dfAnalysis.shape[0]
    print('Cancer Type: ' + '\033[1m{:10s}\033[0m'.format(CancerType))
    print('Variable for analysis: ' + '\033[1m{:10s}\033[0m'.format(ClassVar))
    print('Total samples: ' + '\033[1m{:d}\033[0m'.format(totalsamples))
    ClassVarLevelsFreqTab
    #print('{:_<17} : {:5d}'.format('Total samples',ClassVarLevelsFreqTab.Frequency.sum()))
    
    # Keep samples related to Tumor cells only if CancerStatus is not the class var. 
    if ClassVar != 'CancerStatus':
        toKeep = ['Primary solid Tumor']
#            print('\nKeeping samples concerning "Primary solid Tumor" only.')
        dfAnalysis = OD.FilterLevels(dfAnalysis, 'CancerStatus', toKeep, printStats='no')
    
    # if previous if has been True, we need to print updated stats    
    if totalsamples > dfAnalysis.shape[0]:
#            print('Updated, number of samples in the dataset:' + '\033[1m{:d}\033[0m'.format(dfAnalysis.shape[0])) 
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
        ClassVarLevelsFreqTab, ClassVarLevelsSorted = OD.returnVarLevelsSorted(dfAnalysis,ClassVar)
        ClassVarLevelsFreqTab
        
#        print('\nVarLevelsToKeep = {0}'.format(ClassVarLevelsSorted))
    
    # Keep samples only for the values in VarLevelsToKeep while samples corresponding to the rest are filtered out.
    dfAnalysis_fl = OD.FilterLevels(dfAnalysis, ClassVar, VarLevelsToKeep, printStats='no')
    
    ClassVarLevelsFreqTab, ClassVarLevelsSorted = OD.returnVarLevelsSorted(dfAnalysis_fl,ClassVar)
    
    #print(ClassVarLevelsFreqTab)
#        ClassVarLevelsFreqTab
    
    # check if there are at least 10 (or 10%, whichever is higher) of samples in each class
#    if ((ClassVarLevelsFreqTab['Frequency'].min() < max(10,dfAnalysis.shape[0]/10)) or (ClassVarLevelsFreqTab.shape[0] < 2)):
#    if ((ClassVarLevelsFreqTab['Frequency'].min() < 10) or (ClassVarLevelsFreqTab.shape[0] < 2)):
#        continue
    
    #print('Finally, number of samples in the dataset:' + '\033[1m{:d}\033[0m'.format(dfAnalysis_fl.shape[0])) 
    dfAnalysis_fl = OD.prepareDF(dfAnalysis_fl, ClassVar)
    
    #print('Size of the dataframe: {0}'.format(dfAnalysis_fl.shape))
    
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
        print('********************************************************************')
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
                  .format(len(low_tpm_genes),med_tpm_threshold))
        print(low_tpm_genes)
        # Remove low-TPM genes
        dfAnalysis_fl_cd = OD.CleanData(dfAnalysis_fl,med_tpm_threshold)
        print('\nSize of the dataframe after filtering low-TPM genes: {0}' \
              .format(dfAnalysis_fl_cd.shape))
    else:
        # Don't remove any genes
        print('No genes were removed from the dataset.')
        dfAnalysis_fl_cd = dfAnalysis_fl
    
    # Perform label encoding for the ClassVar and log-transform data
    dfAnalysis_fl_cd, ClassVarEncOrder = OD.mapClassVar(dfAnalysis_fl_cd,ClassVar)
    X, y = OD.fitLogTransform(dfAnalysis_fl_cd,logTransOffset)
    
    
    
    print('Performing ranking of the genes.')
    
    geneNames = dfAnalysis_fl_cd.columns[1:].tolist()
    ranks = {}
    
    #lr = LinearRegression(normalize=True)
    #lr.fit(X, y)
    #ranks['LinearReg'] = OD.Ranks2Dict(np.abs(lr.coef_), geneNames)
    
    f, pval  = f_regression(X, y, center=True)
    ranks['LinearCorr'] = OD.Ranks2Dict(f, geneNames)
    print('LinearCorr complete.')
    
    mine = MINE()
    mic_scores = []
    for i in range(X.shape[1]):
        mine.compute_score(X[:,i], y)
        m = mine.mic()
        mic_scores.append(m)
    ranks['MaxInfoCont'] =  OD.Ranks2Dict(mic_scores, geneNames) 
    print('MaxInfoCont complete.')
    
    # for random forest methods, use floor(sqrt(numfeats)) as the number of estimators
    num_est = int(X.shape[1]**0.5)
    
    rfc = RandomForestClassifier(n_estimators=num_est, random_state=RS) 
    rfc.fit(X,y)
    ranks['RandomForest'] = OD.Ranks2Dict(rfc.feature_importances_, geneNames)
    print('RandomForest complete.')
    
    svmSVC = svm.SVC(kernel='linear')
    svmSVC.fit(X,y)
    ranks['SVMlinear'] = OD.Ranks2Dict(np.abs((svmSVC.coef_.T).T[0]), geneNames)
    print('SVMlinear complete.')
     
#        # Computing genes ranking based on their individual performance using the SVM model.
#        indGeneScores = []
#        for i in range(X.shape[1]):
#             score = cross_val_score(svmSVC, X[:, i:i+1], y, scoring='accuracy',
#                                     cv=10, n_jobs=-1)#StratifiedKFold(n_splits=10, shuffle=True))#ShuffleSplit(len(X), 3, .3))
#             indGeneScores.append((round(np.mean(score), 10)))#, geneNames[i]))
#        ranks['SVMlinearindvGenes'] = dict(zip(geneNames, indGeneScores ))
#        print('SVMlinearindvGenes complete.')
    
    xgb = XGBClassifier()
    xgb.fit(X, y)
    ranks['XGBoostCLF'] = OD.Ranks2Dict(xgb.feature_importances_, geneNames)
    print('XGBoostCLF complete.')
    
    lda =  LinearDiscriminantAnalysis()#(solver='eigen',shrinkage='auto')
    lda.fit(X, y)
    ranks['LDA'] = OD.Ranks2Dict(np.abs((lda.coef_.T).T[0]), geneNames)
    print('LDA complete.')
    
#        # Computing genes ranking based on their individual performance using the LDA model.
#        indGeneScores = []
#        for i in range(X.shape[1]):
#             score = cross_val_score(lda, X[:, i:i+1], y, scoring='accuracy',
#                                     cv=10, n_jobs=-1)#StratifiedKFold(n_splits=10, shuffle=True))#ShuffleSplit(len(X), 3, .3))
#             indGeneScores.append((round(np.mean(score), 10)))#, geneNames[i]))
#        ranks['LDAindvGenes'] = dict(zip(geneNames, indGeneScores ))
#        print('LDAindvGenes complete.')
    
    ridgeCLF = RidgeClassifier()
    ridgeCLF.fit(X, y)
    ranks['RidgeCLF'] = OD.Ranks2Dict(np.abs((ridgeCLF.coef_.T).T[0]), geneNames)
    print('RidgeCLF complete.')
    
    # Ridge and RidgeClassifier both provide same ranking as per my observation. So keeping only one.
    #ridge = Ridge()
    #ridge.fit(X, y)
    #ranks['Ridge'] = OD.Ranks2Dict(np.abs(ridge.coef_), geneNames)
    
    extc = ExtraTreesClassifier(n_estimators=num_est, random_state=RS)
    extc.fit(X,y)
    ranks['ExTreeCLF'] = OD.Ranks2Dict(extc.feature_importances_, geneNames)
    print('ExTreeCLF complete.')
    
    AdabCLF = AdaBoostClassifier(n_estimators=num_est)
    AdabCLF.fit(X,y)
    ranks['adaBoostCLF'] = OD.Ranks2Dict(AdabCLF.feature_importances_, geneNames)
    print('AdaBoostCLF complete.')
    
    rlasso = RandomizedLasso(alpha='bic')
    rlasso.fit(X, y)
    ranks['StabilityRandLasso'] = OD.Ranks2Dict(np.abs(rlasso.scores_), geneNames)
    print('StabilityRandLasso complete.')
    
    r = {}
    for name in geneNames:
        r[name] = round(np.mean([ranks[method][name] for method in ranks.keys()]), 10)
     
    #methods = sorted(ranks.keys())
    ranks['Average'] = r
    #methods.append('Average')
    
#    rfeSVM = RFE(svm.SVC(kernel='linear'), n_features_to_select=1)
#    rfeSVM.fit(X,y)
#    ranks['rfeSVM'] = dict(zip(geneNames, rfeSVM.ranking_ ))
#    print('rfeSVM complete.')
    #methods.append('SVMrfe')
    #
    #rfe = RFE(extc, n_features_to_select=1)
    #rfe.fit(X,y)
    #ranks['rfeExtraTree'] = dict(zip(geneNames, rfe.ranking_ ))
    #methods.append('rfeExtraTree')
    
#    rfeRFC = RFE(RandomForestClassifier(n_estimators=200, random_state=RS), n_features_to_select=1)
#    rfeRFC.fit(X,y)
#    ranks['rfeRFC'] = dict(zip(geneNames, rfeRFC.ranking_ ))#OD.Ranks2Dict(np.array(rfeRFC.ranking_, dtype=float), geneNames)#, order=1)
#    print('rfeRFC complete.')    

    dfRanks = pd.DataFrame.from_dict(ranks)
    dfRanks.reset_index(inplace=True)
    dfRanks.rename(columns={'index':'GeneNames'}, inplace=True)
    dfRanks.sort_values(by='Average', inplace=True, ascending=False)
    
    print('\nDone!\n')
    
    
    models = [ExtraTreesClassifier(n_estimators=num_est, random_state=RS), # 0
              RandomForestClassifier(n_estimators=num_est, random_state=RS), # 1
              AdaBoostClassifier(n_estimators=num_est), # 2
              XGBClassifier(), # 3
              LinearDiscriminantAnalysis(), # 4
              RidgeClassifier(), # 5
              KNeighborsClassifier(20), # 6
              svm.SVC(kernel='linear'), # 7
              #svm.SVC(kernel='rbf', C=5) # 8  
             ]
    # Run models for crossvalidated average score of accuracy
    CV = 'Validation: SKF'
    shuffle = True
    scoring = 'accuracy'
    folds = 10
    #dfCVscores = OD.CVScorer([svm.NuSVC()], CV, X, y, scoring, shuffle)
    print('Performing models CV analysis using accuracy.')
    dfCVscoresAUC = OD.CVScorer(models, CV, X, y, scoring, shuffle, folds)
    #dfCVscores = OD.CVScorer([models[i] for i in [3, 6]], CV, X, y, scoring, shuffle)
    print('Done!')
    #print('Cancer Type: ', CancerType)
    #dfCVscoresAUC
    
    if len(VarLevelsToKeep) == 2:
        scoring = 'roc_auc'
        print('Performing models CV analysis using area under the ROC curve.')
        dfCVscoresROC = OD.CVScorer(models, CV, X, y, scoring, shuffle, folds)
        print('Done!')
    else:
        print('Skipping CV analysis using area under the ROC curve. ' \
              'This is possible for binary problems only.')
    
    
        
    print('Writing dataset, genees ranking and CV analysis results to a' \
          'directory named:{0}' \
          .format(CancerType))
    os.makedirs(CancerType , exist_ok=True)
    
    ClassVarName = ClassVar.split(sep='_')[1]
    
#    writer = pd.ExcelWriter(CancerType + '/' + CancerType + '-' + ClassVarName + 'GenesRanking.xlsx', engine='xlsxwriter')
    
    # Convert the dataframe to an XlsxWriter Excel object.
#    dfRanks.to_excel(writer, sheet_name='GenesRanking', index=False, float_format='%.10f')
#    dfCVscoresAUC.to_excel(writer, sheet_name='CVscoresAccuracy', index=False, float_format='%.10f')
    if len(VarLevelsToKeep) == 2:
#        dfCVscoresROC.to_excel(writer, sheet_name='CVscoresAreaUnderROC', index=False, float_format='%.10f')
        dfCVscoresROC.to_csv(CancerType + '/' + CancerType + '-' + ClassVarName + 'CVscoresAreaUnderROC.csv', index=False)
    
#    writer.save()
    
    dfRanks.to_csv(CancerType + '/' + CancerType + '-' + ClassVarName + 'GenesRanking.csv', index=False)
    dfCVscoresAUC.to_csv(CancerType + '/' + CancerType + '-' + ClassVarName + 'CVscoresAccuracy.csv', index=False)
#        dfAnalysis_fl_cd.to_csv(CancerType + '/' + CancerType + '-' + ClassVarName + 'Data.csv', index=False)
    
    print('Done!')
