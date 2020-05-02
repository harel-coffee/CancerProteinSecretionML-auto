#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: azams
"""

#%%
import Omics.OmicsData as OD
#import os
#import docx


import numpy as np
import pandas as pd

from sklearn.cross_validation import cross_val_score
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.linear_model import RidgeClassifier, RandomizedLasso #, RandomizedLogisticRegression,  Ridge, Lasso, LinearRegression
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier #, AdaBoostClassifier
#from sklearn.tree  import DecisionTreeClassifier
from xgboost import XGBClassifier
#from sklearn.neighbors import KNeighborsClassifier
from sklearn import svm
from sklearn.feature_selection import RFE, f_regression #, VarianceThreshold
#from sklearn.model_selection import KFold, StratifiedKFold, cross_val_score

from minepy import MINE
RS = 20170628

#%%
CancerType = 'TCGA-THCA'
ClassVar = 'CancerStatus'

# two possibilities for data, OriginalData or CleanedData
data = 'CleanedData'

dfAnalysis= pd.read_csv(CancerType + '/' + CancerType + '-' + ClassVar + data +'.csv')


print("Performing ranking of the genes.")

geneNames = dfAnalysis.columns[1:].tolist()
array = dfAnalysis.values
X = array[:,1:len(dfAnalysis.columns)]
y = array[:,0]
ranks = {}

#lr = LinearRegression(normalize=True)
#lr.fit(X, y)
#ranks["LinearReg"] = OD.Ranks2Dict(np.abs(lr.coef_), geneNames)

f, pval  = f_regression(X, y, center=True)
ranks["LinearCorr"] = OD.Ranks2Dict(f, geneNames)

print("Ranking using Linear Correlation completed...")

mine = MINE()
mic_scores = []
for i in range(X.shape[1]):
    mine.compute_score(X[:,i], y)
    m = mine.mic()
    mic_scores.append(m)
ranks["MaxInfoCont"] =  OD.Ranks2Dict(mic_scores, geneNames) 

print("Ranking using MaxInfoCont completed...")


rfc = RandomForestClassifier(n_estimators=200, random_state=RS) 
rfc.fit(X,y)
ranks["RandomForest"] = OD.Ranks2Dict(rfc.feature_importances_, geneNames)

print("Ranking using RandomForest completed...")


svmSVC = svm.SVC(kernel='linear')
svmSVC.fit(X,y)
ranks["SVMlinear"] = OD.Ranks2Dict(np.abs((svmSVC.coef_.T).T[0]), geneNames)

print("Ranking using SVMlinear completed...")


# Computing genes ranking based on their individual performance using the SVM model.
indGeneScores = []
for i in range(X.shape[1]):
     score = cross_val_score(svmSVC, X[:, i:i+1], y, scoring='accuracy', cv=10)#, n_jobs=-1)#StratifiedKFold(n_splits=10, shuffle=True))#ShuffleSplit(len(X), 3, .3))
     indGeneScores.append((round(np.mean(score), 10)))#, geneNames[i]))
ranks["SVMlinearindvGenes"] = dict(zip(geneNames, indGeneScores ))

print("Ranking using SVMlinearindvGenes completed...")


xgb = XGBClassifier()
xgb.fit(X, y)
ranks["XGBoostCLF"] = OD.Ranks2Dict(xgb.feature_importances_, geneNames)    

print("Ranking using XGBoostCLF completed...")


lda =  LinearDiscriminantAnalysis()
lda.fit(X, y)
ranks["LDA"] = OD.Ranks2Dict(np.abs((lda.coef_.T).T[0]), geneNames)

print("Ranking using LDA completed...")


# Computing genes ranking based on their individual performance using the LDA model.
indGeneScores = []
for i in range(X.shape[1]):
     score = cross_val_score(lda, X[:, i:i+1], y, scoring='accuracy', cv=10)#, n_jobs=-1)#StratifiedKFold(n_splits=10, shuffle=True))#ShuffleSplit(len(X), 3, .3))
     indGeneScores.append((round(np.mean(score), 10)))#, geneNames[i]))
ranks["LDAindvGenes"] = dict(zip(geneNames, indGeneScores ))

print("Ranking using LDAindvGenes completed...")

    
ridgeCLF = RidgeClassifier()
ridgeCLF.fit(X, y)
ranks["RidgeCLF"] = OD.Ranks2Dict(np.abs((ridgeCLF.coef_.T).T[0]), geneNames)

print("Ranking using RidgeCLF completed...")

# Ridge and RidgeClassifier both provide same ranking as per my observation. So keeping only one.
#ridge = Ridge()
#ridge.fit(X, y)
#ranks["Ridge"] = OD.Ranks2Dict(np.abs(ridge.coef_), geneNames)

extc = ExtraTreesClassifier(n_estimators=200, random_state=RS)
extc.fit(X,y)
ranks["ExTreeCLF"] = OD.Ranks2Dict(extc.feature_importances_, geneNames)

print("Ranking using RidgeCLF completed...")

#AdabCLF = AdaBoostClassifier(n_estimators=200)
#AdabCLF.fit(X,y)
#ranks["adaBoostCLF"] = OD.Ranks2Dict(AdabCLF.feature_importances_, geneNames)

rlasso = RandomizedLasso(alpha='bic')
rlasso.fit(X, y)
ranks["StabilityRandLasso"] = OD.Ranks2Dict(np.abs(rlasso.scores_), geneNames)

print("Ranking using StabilityRandLasso completed...")


r = {}
for name in geneNames:
    r[name] = round(np.mean([ranks[method][name] for method in ranks.keys()]), 2)
 
#methods = sorted(ranks.keys())
ranks["Average"] = r

print("Ranking average completed...")
#methods.append("Average")

rfeSVM = RFE(svm.SVC(kernel='linear'), n_features_to_select=1)
rfeSVM.fit(X,y)
ranks["rfeSVM"] = dict(zip(geneNames, rfeSVM.ranking_ ))
print("Ranking using rfeSVM completed...")
#methods.append("SVMrfe")
#
#rfe = RFE(extc, n_features_to_select=1)
#rfe.fit(X,y)
#ranks["rfeExtraTree"] = dict(zip(geneNames, rfe.ranking_ ))
#methods.append("rfeExtraTree")

rfeRFC = RFE(RandomForestClassifier(n_estimators=200, random_state=RS), n_features_to_select=1)
rfeRFC.fit(X,y)
ranks["rfeRFC"] = dict(zip(geneNames, rfeRFC.ranking_ ))#OD.Ranks2Dict(np.array(rfeRFC.ranking_, dtype=float), geneNames)#, order=1)
print("Ranking using rfeRFC completed...")


dfRanks = pd.DataFrame.from_dict(ranks)
dfRanks.reset_index(inplace=True)
dfRanks.rename(columns={'index':'GeneNames'}, inplace=True)
dfRanks.sort_values(by='Average', inplace=True, ascending=False)

print("Done!")
#%%
    
print("Writing genes ranking results to the directory named:{0}".format(CancerType))

writer = pd.ExcelWriter(CancerType + '/' + CancerType + '-' + ClassVar + '-' + data + 'GenesRanking.xlsx', engine='xlsxwriter')
dfRanks.to_excel(writer, sheet_name='GenesRanking', index=False, float_format='%.10f')
writer.save()

dfRanks.to_csv(CancerType + '/' + CancerType + '-' + ClassVar + '-' + data + 'GenesRanking.csv', index=False)

print("Done!")
