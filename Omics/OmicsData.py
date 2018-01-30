#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 09:50:10 2017

@author: azams
Updated by Jonathan L. Robinson, 2018-01-30
    - Generalized to use data beyond PSN transcriptomics.
    - Removed/commented older code no longer relevant or used.
    

code to autoreload 
%load_ext autoreload
%autoreload 2

OR run the following
import importlib
importlib.reload(OD)
"""

#print(__doc__)
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt
from sklearn.preprocessing import LabelEncoder, StandardScaler
import math as mt
from scipy import interp
from sklearn.metrics import roc_curve, auc
from itertools import cycle

from sklearn.model_selection import KFold, StratifiedKFold, cross_val_score
#from sklearn.linear_model import LogisticRegression, RidgeClassifier, SGDClassifier
#from sklearn.ensemble import RandomForestClassifier
#from sklearn.tree import DecisionTreeClassifier
#from sklearn.naive_bayes import MultinomialNB
#from sklearn.neighbors import KNeighborsClassifier
#dfPSN = pd.read_csv('BigDataProject/PSN_gene_data.csv')
#dfSEC = pd.read_csv('BigDataProject/SEC_gene_data.csv')

###############################################################################

###############################################################################

def ReadOMICSdataCSV(fName):

    """
    Reads in the omics data CSV file,
    removes the redundant column 'Unnamed: 0'
    and returns the data as a dataframe.
    
    """   
    df = pd.read_csv(fName + ".csv")
    # There is a column with name Unnamed: 0. 
    # Dropping it here.
    if 'Unnamed: 0' in df.columns:
        df.drop('Unnamed: 0', axis=1, inplace=True)

    return df
###############################################################################
###############################################################################


#def prepareDFCancerType(dfSpecificCancerData):
#
#    """
#    Given the cancer data which is already filtered for a specific 
#    cancer type (optional, NOT necessary!) is passed as a dataframe, it
#     i) separates the genes (features) from the class variables (particular to cancer types only), 
#        such as: TumorStage, Race, and Gender
#    ii) asks the user to choose one of the class variables
#    iii)removes the rest of the class variables and adds the chosen variable as the first column(s), followed by
#    the data for all the genes (features) and returns as a dataframe ready to work on.
#    """
#    
#    # Determine the number of genes (features) and class variables in the dataframe
#    # Note that the dataframe is arranged such that "CancerStatus" is the first
#    # class variable, so all columns before "CancerStatus" are genes.
#    numFeatures = dfSpecificCancerData.columns.get_loc('CancerStatus')
#    numClassVars = len(dfSpecificCancerData.columns) - numFeatures
#    
#    features = dfSpecificCancerData.iloc[:, 0:numFeatures]
##    f_names = list(features.columns)
#    
#    targets = dfSpecificCancerData.iloc[:, 577:len(dfSpecificCancerData.columns)]
#    t_names = list(targets.columns)
#    print("\n*******************************************************************************")
#    while True:
#        ClassVar = input("Choose a class variable (" + ' '.join(t_names) + "): ")
#        if ClassVar in t_names:
#            break
#        else:
#            print("Please splell correctly!")
#            
#    print("\n*******************************************************************************")
#    target = targets[ClassVar]
#    df = features
#    df[ClassVar] = target
#    
#    # Class variable is the last column, bringing it to the first place.
#    cols = df.columns.tolist()
#    cols = cols[-1:] + cols[:-1]
#    df = df[cols]
#
#    return df
###############################################################################

###############################################################################


def prepareDFgeneral(dfAllCancerData):

    """
    Given that the entire PSN cancer data is passed as a dataframe, it
    i) separates the genes (features) from the target variables, 
    that is: CancerStatus, Project, TumorStage, Race, and Gender
    ii) asks the user to chose one of them as the class variable
    iii)removes the rest and adds the chosen variable as the first column, followed by
    the entire genes (features) and returns as a dataframe ready to work on.    
    """      
    
    # Determine the number of genes (features) in the dataframe
    # Note that the dataframe is arranged such that "CancerStatus" is the first
    # class variable, so all columns before "CancerStatus" are genes.
    numFeatures = dfAllCancerData.columns.get_loc('CancerStatus')
    
    features = dfAllCancerData.iloc[:, 0:numFeatures]
#    f_names = list(features.columns)
    
    targets = dfAllCancerData.iloc[:, numFeatures:len(dfAllCancerData.columns)]
    t_names = list(targets.columns)
    print("\n*******************************************************************************")
    while True:
        ClassVar = input("Choose a class variable (" + ' '.join(t_names) + "): ")
        if ClassVar in t_names:
            break
        else:
            print("Please splell correctly!")
            
    print("\n*******************************************************************************")
    target = targets[ClassVar]
    df = features
    df[ClassVar] = target
    
    # Class variable is the last column, bringing it to the first place.
    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df = df[cols]

    return df
###############################################################################
def prepareDF(dfAllCancerData, ClassVar):

    """
    Given that the entire PSN cancer data is passed as a dataframe, it
    i) separates the genes (features) from the target variables, 
    that is: CancerStatus, Project, TumorStage, Race, and Gender
    ii) keeps the column corresponding to ClassVar and removes the rest
    iii) and moves it to be the first column, followed by
    the entire genes (features) and returns as a dataframe ready to work on.    
    """
    
    # Determine the number of genes (features) in the dataframe
    # Note that the dataframe is arranged such that "CancerStatus" is the first
    # class variable, so all columns before "CancerStatus" are genes.
    numFeatures = dfAllCancerData.columns.get_loc('CancerStatus')
    
    features = dfAllCancerData.iloc[:, 0:numFeatures]
    target = dfAllCancerData[ClassVar]
    df = features
    df[ClassVar] = target
    
    # Class variable is the last column, bringing it to the first place.
    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df = df[cols]

    return df

###############################################################################

#from Pramod
def prepareDF_Mod(dfAllCancerData, TargetVariable):

    """
    Given that the PSN cancer data is passed as a dataframe, it
    i) separates the Genes (features) from the target variables,
    ii) asks the user to chose one of them as the class variable
    iii) adds it as the first column of features and returns a dataframe
        ready to work on.
    """
    
    # Determine the number of genes (features) in the dataframe
    # Note that the dataframe is arranged such that "CancerStatus" is the first
    # class variable, so all columns before "CancerStatus" are genes.
    numFeatures = dfAllCancerData.columns.get_loc('CancerStatus')
    
    features = dfAllCancerData.iloc[:, 0:numFeatures]   
    CancerStatuses = dfAllCancerData[TargetVariable]
    df = features    
    df[TargetVariable] = CancerStatuses
    
    # Class variable is the last column, bringing it to the first place.
    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df = df[cols]

    return df
###############################################################################

###############################################################################

def printNaNs(df):

    """
    Given that the PSN cancer data is passed as a dataframe, it
    i) prints the number of missing values (if any) in each of the columns.
    ii)reports if no missing data is found.
    """   
    print("\n********************************************************************")
    print("Number of samples in the dataset: {0}".format(df.shape[0]))
    print("********************************************************************")
    print("Printing missing values count (if any) in each of the columns. ")
    flag = True
    for c in df.columns:
        if df[c].isnull().sum():
            print("{:_<12} : {:5d}".format(c,df[c].isnull().sum()))
            flag = False
    if flag:
        print('Hurrah, No missing data right now!')
    print("********************************************************************")
###############################################################################

###############################################################################
def dropNaNs(df, ClassVar='none'):

    """
    Given the omics data passed as a dataframe, and (optionally) the name of a
    class variable, it
    
    i) prints the total number of samples in the dataset.
    ii) if none of the samples have any missing values, it returns the same dataframe, else:
        a) number of samples having missing values are reported
        b) these samples are removed from the dataset
        c) number of samples remained in the dataset after removing with missing values, are reported
        d) returns the updated dataframe
    """   
    print("\n********************************************************************")
    print("Number of samples in the dataset: {0}".format(df.shape[0]))
    print("********************************************************************")

    if ClassVar == 'none':
        dfdna = df.dropna()
    else:
        dfdna = df.dropna(subset=[ClassVar])

    if df.shape[0] > dfdna.shape[0]:
        print("Number of samples having missing values: {0}".format(df.shape[0]- dfdna.shape[0]))
        print("Number of samples remained after dropping samples with missing data: {0}".format(dfdna.shape[0]))
    else:
        print("There are no samples with missing values!")
        
    return dfdna        
###############################################################################

###############################################################################
def printClassVarValCounts(df, ClassVar):

    """
    Given that the PSN cancer data is passed as a dataframe, and Class variable as string, it
    i) prints the total number of samples in the dataset
    ii) Prints all distinct values and the number of samples corresponding to these values
    iii) It also displays the total number of missing values(if any) as NaN
    """        
    print("\n********************************************************************")
    print("Number of samples in the dataset: {0}".format(df.shape[0]))
    print("Target variable, {0}, has {1} unique values,".format(ClassVar, len(df[ClassVar].unique())))
    print("with the following distribution of the data.")
    print(df[ClassVar].value_counts(dropna=False))
    print("********************************************************************")

###############################################################################

###############################################################################
# If there are some levels of the Class Variable that we want to exclude,
# we can use this method.
###############################################################################
def RemoveExtraLevels(df, ClassVar, toRemove):

    """
    Given that the PSN cancer data is passed as a dataframe, Class variable as string, 
    and a list of values of class variable which need to be removed
    i) prints the total number of samples in the dataset
    ii) Prints all distinct values and the number of samples corresponding to these values
    iii) It also displays the total number of missing values(if any) as NaN
    """        
    for x in toRemove:
        df.drop(df.index[df[ClassVar] == x], inplace= True)
    printClassVarValCounts(df, ClassVar)
    return df
###############################################################################

###############################################################################
def FilterLevels(df, ClassVar, toKeep, printStats='yes'):

    """
    Given the cancer data as a dataframe, Class variable as string, 
    and a list of values of that class variable which should be kept:
    i) prints the total number of samples in the dataset
    ii) Prints all distinct class variable values and the number of samples
        corresponding to these values
    iii) It also displays the total number of missing values(if any) as NaN
    """
    df_new = pd.DataFrame()
    for x in toKeep:
        df_temp = df[df[ClassVar] == x]
        df_new = df_new.append(df_temp)
    if printStats=='yes':
        printClassVarValCounts(df_new, ClassVar)
    return df_new
###############################################################################

###############################################################################
def returnVarLevels(df, var):

    """
    Returns the unique values/levels of the given variable.
    """
    
    return df[var].unique().tolist()
    
###############################################################################

###############################################################################
def mapClassVar(dfdna, ClassVar):

    """
    Pass it a dataframe, dfdna, after removing NaNs using dropNaNs(df), 
    and class variable, it will
    i) map the levels of string levels of the variable to integers
    ii) apply this mapping to dfdna 
    iii) return the new df
    iv) print the mapping. 
    """
    
    le = LabelEncoder()
    df_le = dfdna.copy()
    df_le[ClassVar] = le.fit_transform(df_le[ClassVar])
    ClassVarLvlsEnc = df_le[ClassVar].unique().tolist()
    ClassVarLvlsEnc.sort()
    ClassVarLvlsIDs = le.inverse_transform(ClassVarLvlsEnc).tolist()
    print("\n***************************************************************")
    print('Following label encoding has been assigned to the values of {0}.'.format(ClassVar))
    dictionary = dict(zip(ClassVarLvlsEnc,ClassVarLvlsIDs))
    print(dictionary)
    print("\n***************************************************************")
#    printClassVarValCounts(df, ClassVar)
    return df_le, ClassVarLvlsIDs
###############################################################################

############################################################################### 
def fitScalarTransform(df):

    """
    Standardize the data so that variance is 1 and mean is zero.
    Returns X_scaled and y. 
    y is the column corresponding to the class variable.
    X_scaled contains are all other variables on which scaling is applied.
    contains
    """
    array = df.values
    X = array[:,1:len(df.columns)]
#    y = array[:,0]
# above way of getting y changes y to floats, whereas y is simply the label 
# encoded class variable. The line below is to remedy this. 
# Doing so, we get y as int64 which is required at some places e.g., when plotting t-SNE results.     
    y = np.asarray(df.iloc[:,0])
    scaler = StandardScaler().fit(X)
#    scaler
#    scaler.mean_                                      
#    scaler.scale_                                       
    X_scaled = scaler.transform(X)                               
    return X_scaled, y
    
###############################################################################
from sklearn.preprocessing import FunctionTransformer

###############################################################################
def fitLogTransform(df):
    array = df.values
    X = array[:,1:len(df.columns)]
#    y = array[:,0]
# above way of getting y changes y to floats, whereas y is simply the label 
# encoded class variable. The line below is to remedy this. 
# Doing so, we get y as int64 which is required at some places e.g., when plotting t-SNE results.     
    y = np.asarray(df.iloc[:,0])

    logScaler = FunctionTransformer(np.log1p)
    X_scaled = logScaler.transform(X)
    return X_scaled, y
###############################################################################

#from Pramod
###############################################################################
def dffitLogTransform(df):
    
    """
    Takes a dataframe with the first column as the classificatin variable and 
    gene expression levels as the rest of the columns and returns a new dataframe
    with log transformed gene expression levels
    """   
    gene_names = df.ix[:,0:].columns.values    
    df_new = pd.DataFrame(index=range(len(df)))
    logScaler = FunctionTransformer(np.log1p)
        
    for gene in gene_names:
        X = df[gene]
        X_scaled = logScaler.transform(X)
        df_new[gene] = X_scaled.reshape(-1,1)
          
    return df_new
    
###############################################################################
def PrepareLogitResults(df, ClassVar):
 
    """
    Pass it a dataframe, it fits a logit model using ClassVar and then returns 
    the following model parameters:
    'Beta', 'p-Value', 'OR', 'CI (2.5%)', 'CI (97.5%)'   
    """
    df['intercept']=1.0
    train_cols=df.columns[1:]
    res = sm.Logit(df[ClassVar], df[train_cols]).fit(maxiter=10000, method='ncg')#'ncg') #bfgs
    params = res.params
    conf = res.conf_int()
    conf['OR'] = params
    conf.columns = ['CI (2.5%)', 'CI (97.5%)', 'OR']
    conf = np.exp(conf)
    conf['p-Value'] = res.pvalues
    conf['Beta'] = res.params.values
    cols_order = ['Beta', 'p-Value', 'OR', 'CI (2.5%)', 'CI (97.5%)']
    conf = conf[cols_order] 
    conf.reset_index(level=0, inplace=True)
    conf = conf.rename(columns={'index':'Variable'})
    return conf
###############################################################################
from sklearn.preprocessing import MinMaxScaler
###############################################################################
def Ranks2Dict(ranks, names, order=1):
    minmax = MinMaxScaler()
    ranks = minmax.fit_transform(order*np.array([ranks]).T).T[0]
    ranks = map(lambda x: round(x, 10), ranks)
    return dict(zip(names, ranks ))

##############################################################################
###############################################################################
def PrepareCorrResults(df):

    """
    Pass it a dataframe, and it returns Pairwise Pearson correlation coefficient values 
    for the entire variables of the datadframe.
    The first columns of the returned dfCORR contains correlation values of the classvariable
    versus all other variables. 
    """
    dfCORR = df.corr()    
    dfCORR.reset_index(level=0, inplace=True)
    dfCORR = dfCORR.rename(columns={'index':'Variable'})
    return dfCORR
###############################################################################
def CVScorer(models, CV, X, y, scoring, shuffle, folds=10):
    
    if CV == 'Validation: SKF':
        cv = StratifiedKFold(n_splits=folds, shuffle=shuffle)
    elif CV == 'Validation: KF':
        cv = KFold(n_splits=folds, shuffle=shuffle)
    
    dfCVscores = pd.DataFrame(columns=['Model', 'Scoring', 'Score', 'CI-lower', 'CI-high'])    
    for model in models:
        modelName = str(model).partition('(')[0]
        scores = cross_val_score(model, X, y, scoring=scoring, cv=cv) #, n_jobs=-1) this was causing problem for LDA and RidgeClassifier
        dfCVscores = dfCVscores.append(pd.Series([modelName, scoring, scores.mean(),(scores.mean() - 2*scores.std()), (scores.mean() + 2*scores.std())],
                                                 index=dfCVscores.columns), ignore_index=True)
        #print("{3} [-/+]: {0:.2f} [{1:.2f}, {2:.2f}]".format(scores.mean(), 
        #    (scores.mean() - 2*scores.std()), (scores.mean() + 2*scores.std()), 'Model: ' + modelName + ', Cross validated average score of ' + scoring))
        
    return dfCVscores
###############################################################################
def ROCanalysis(mod_name, CV, classifier, X, y, shuffle, folds=10):
    
  """
  Plot ROC curve generated using 10-fold cross validation for the given model.
  mod_name:: Name of the classifier.
  CV:: chose one of these: 'Validation: SKF', 'Validation: KF'
  classifier:: e.g., LogisticRegression()
  X:: Featuers/variables
  y:: the class variable
  shuffle:: True
  """
# Classification and ROC analysis

# Run classifier with cross-validation and plot ROC curves
#  array = df.values
#  X = array[:,1:len(df.columns)]
#  y = array[:,0]

  if CV == 'Validation: SKF':
      cv = StratifiedKFold(n_splits=folds, shuffle=shuffle)
  elif CV == 'Validation: KF':
      cv = KFold(n_splits=folds, shuffle=shuffle)

  mean_tpr = 0.0
  mean_fpr = np.linspace(0, 1, 101)
  mean_acc = []
  tprs = []
  
  colors = cycle(['darkcyan', 'indigo', 'darkgreen', 'darkgoldenrod', 'darkblue'
    , 'darkorange', 'mediumvioletred', 'crimson', 'darksalmon', 'darkred'])
  lw = 2
  plt.figure(figsize=(8, 8))
  i = 0
  for (train, test), color in zip(cv.split(X, y), colors):
      
      if mod_name.startswith('Ridge'):
          classifier.fit(X[train], y[train])
          confScores = classifier.decision_function(X[test])
          fpr, tpr, thresholds = roc_curve(y[test], confScores, pos_label=1)
      else:
          probas_ = classifier.fit(X[train], y[train]).predict_proba(X[test])
          fpr, tpr, thresholds = roc_curve(y[test], probas_[:, 1], pos_label=1)
          
#      model = classifier.fit(X[train], y[train])
                                                    #[:, 1]
#      fpr, tpr, thresholds = roc_curve(y[test], probas, pos_label=1)
      
      mean_acc.append(classifier.score(X[test],y[test]))
      # Compute ROC curve and area the curve

      mean_tpr += interp(mean_fpr, fpr, tpr)
      mean_tpr[0] = 0.0
      tp = interp(mean_fpr, fpr, tpr)
      tp[0]=0.0
      tprs.append(tp)
             
      roc_auc = auc(fpr, tpr)
      plt.plot(fpr, tpr, alpha=0.55, lw=lw, color=color,
             label='ROC fold %d (area = %0.2f)' % (i+1, roc_auc))

      i += 1
  plt.plot([0, 1], [0, 1], linestyle='--', lw=lw, color='k')#,
         #label='Luck')
#
  tprs = np.array(tprs)
  mean_tprs = tprs.mean(axis=0)
  std = tprs.std(axis=0)
  tprs_upper = np.minimum(mean_tprs + std, 1)
  tprs_lower = mean_tprs - std

#  mean_auc_test = auc(mean_fpr, mean_tprs)
#  print(mean_auc_test)
#  plt.plot(base_fpr, mean_tprs, 'b',label="Mean ROC (area = %0.2f)" % (mean_auc), lw=lw)
  plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=0.4)
#
  mean_tpr /= cv.get_n_splits(X, y)
  mean_tpr[-1] = 1.0
  mean_auc = auc(mean_fpr, mean_tpr)
#  print(mean_auc)
  plt.plot(mean_fpr, mean_tpr, color='b', linestyle=':', 
         label='Mean ROC (area = %0.2f)' % mean_auc, lw=4)

  plt.xlim([-0.01, 1.01])
  plt.ylim([-0.01, 1.01])
  plt.xlabel("False Positive Rate (1 - Specificity) \n Cross-Validation Average"
        + " Score of Accuracy: %0.3f%%" % (np.mean(mean_acc)*100), size=12)
  plt.ylabel('True Positive Rate (Sensitivity)', size=12)
  plt.title("Receiver Operating Characteristic Curve (%s) \n Model: %s" 
            % (CV, mod_name), size=13)
  plt.legend(loc="lower right")
  plt.grid(True)
  plt.show()
#  plt.savefig("auc.png")
##
#  doc.add_picture("auc.png", width=docx.shared.Cm(20), height=docx.shared.Cm(20))
#  doc.add_paragraph(" Model: " + mod_name + "\n" + CV)
#  return X, y
###############################################################################




def GeneExpression(df,Level):
    """
    This function takes in a data frame with only the values of the Gene expression
    and provides the list of genes whose gene expression values lie in the lower 2 percentile of the population.
    """
    df_Gene=df.iloc[:,1:]
    data_stats = pd.DataFrame()
    median = []
    mean = []
    for i in df_Gene.columns.values:
        median.append(mt.log10(df_Gene[i].median()+1))
        mean.append(mt.log10(df_Gene[i].mean()+1))
    data_stats['Gene Name']=df_Gene.columns.values
    data_stats['Median'] = median
    data_stats['Mean'] = mean
    percentile = np.percentile(data_stats['Median'],Level)
    LowCountGene = []
    for i in data_stats.index.values:
        if data_stats['Median'][i] <= percentile:
            LowCountGene.append(data_stats['Gene Name'][i])           
    
    return data_stats, np.array(LowCountGene)

def CleanData (df, Level):
    data_stats, LowCountGene = GeneExpression(df, Level)
    df_clean = df.drop(LowCountGene,1)
    return df_clean
    
###############################################################################
def prepCancerTypeDict(hdfStore=False):
    """
    This function loads the entire PSN cancer dataset from  'allcancerdata.csv' 
    and returns a disctionaly of dataframes, where each dataframe corresponds 
    to a cancer types. 
    """ 
    
    # Import data from csv to a data frame
    df = ReadOMICSdataCSV("data/allcancerdata")
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
        cancerTypesDic[project]= FilterLevels(df, ClassVar, toKeep, printStats='no')
    
    # For hdfStore=True, we write the dictionay to a hdfStore.
    if hdfStore:
        CancerDataStore = pd.HDFStore('data/CancerDataStore.h5')
        for (key, value) in cancerTypesDic.items():
            # keys are names of cancers, e.g., TCGA-BRCA. Using split to ignore the TCGA- part and use
            # the rest as the name. With prefix TCGA-, it is not a valid Python identifier.
            CancerDataStore.put(key.split('-')[1], value)
            print("{0} successfully saved in store!".format(key))
        
        print(CancerDataStore)
        CancerDataStore.close()
    
    return cancerTypesDic
    
###############################################################################

from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis #, QuadraticDiscriminantAnalysis
###############################################################################

def plotPCA(X, y, nComp, target_names, save=False):
    """
    =======================================================
    Comparison of LDA and PCA 2D projection of PSN dataset
    =======================================================
    
    Principal Component Analysis (PCA) applied to this data identifies the
    combination of attributes (principal components, or directions in the
    feature space) that account for the most variance in the data. Here we
    plot the different samples on the possible pairs of principal components.
    
    Linear Discriminant Analysis (LDA) tries to identify attributes that
    account for the most variance *between classes*. In particular,
    LDA, in contrast to PCA, is a supervised method, using known class labels.
    """
    pca = PCA(n_components=nComp)
    X_r = pca.fit(X).transform(X)
    # Percentage of variance explained for each components
    #print('explained variance ratio (first two components): %s'
    #      % str(pca.explained_variance_ratio_))

    l = len(target_names)
    colors = ['darkcyan', 'indigo', 'darkgreen', 'darkgoldenrod', 'darkblue'
        , 'darkorange', 'mediumvioletred', 'crimson', 'darksalmon', 'darkred', 'cyan', 'orange','green']   
    colors = colors[0:l]
    target_codes = list(range(0,l))
    
    
    plt.figure(figsize=(8, 8))
    lw = 2
    for xComp in range(1,nComp+1):
        for yComp in range(xComp+1,nComp+1):
            for color, i, target_name in zip(colors, target_codes, target_names):
                plt.scatter(X_r[y == i, xComp-1], X_r[y == i, yComp-1], color=color, alpha=.8, lw=lw,
                            label=target_name)
            plt.legend(loc='best', shadow=False, scatterpoints=1)
            plt.title('PCA applied to PSN dataset')
            plt.xlabel('PCA ' + str(xComp))
            plt.ylabel('PCA ' + str(yComp))
            if save:
                plt.savefig('PCA component ' + str(xComp) + ' by ' + str(yComp) + '.png') 
            
            plt.show()
###############################################################################
###############################################################################
def plotLDA(X, y, nComp, target_names, save=False):
    """
    =======================================================
    Comparison of LDA and PCA 2D projection of PSN dataset
    =======================================================
    
    Principal Component Analysis (PCA) applied to this data identifies the
    combination of attributes (principal components, or directions in the
    feature space) that account for the most variance in the data. Here we
    plot the different samples on the possible pairs of principal components.
    
    Linear Discriminant Analysis (LDA) tries to identify attributes that
    account for the most variance *between classes*. In particular,
    LDA, in contrast to PCA, is a supervised method, using known class labels.
    """
    lda = LinearDiscriminantAnalysis(n_components=nComp)
    X_r2 = lda.fit(X, y).transform(X)
    
    l = len(target_names)
    colors = ['darkcyan', 'indigo', 'darkgreen', 'darkgoldenrod', 'darkblue'
        , 'darkorange', 'mediumvioletred', 'crimson', 'darksalmon', 'darkred', 'cyan', 'orange','green']   
    colors = colors[0:l]
    target_codes = list(range(0,l))
    
    
    plt.figure(figsize=(8, 8))
    lw = 2
    for xComp in range(1,nComp+1):
        for yComp in range(xComp+1,nComp+1):
            for color, i, target_name in zip(colors, target_codes, target_names):
                plt.scatter(X_r2[y == i, xComp-1], X_r2[y == i, yComp-1], alpha=.8, color=color, lw=lw,
                            label=target_name)
            plt.legend(loc='best', shadow=False, scatterpoints=1)
            plt.title('LDA applied to PSN dataset')
            plt.xlabel('LDA ' + str(xComp))
            plt.ylabel('LDA ' + str(yComp))
            if save:
                plt.savefig('LDA component ' + str(xComp) + ' by ' + str(yComp) + '.png') 
            
            plt.show()
###############################################################################


###############################################################################
def plotPCAvsLDA(X, y, nComp, target_names, save=False):
    """
    =======================================================
    Comparison of LDA and PCA 2D projection of PSN dataset
    =======================================================
    
    Principal Component Analysis (PCA) applied to this data identifies the
    combination of attributes (principal components, or directions in the
    feature space) that account for the most variance in the data. Here we
    plot the different samples on the possible pairs of principal components.
    
    Linear Discriminant Analysis (LDA) tries to identify attributes that
    account for the most variance *between classes*. In particular,
    LDA, in contrast to PCA, is a supervised method, using known class labels.
    """
    pca = PCA(n_components=nComp)
    X_r = pca.fit(X).transform(X)
    
#    qda = QuadraticDiscriminantAnalysis(n_components=nComp)
#    X_r = qda.fit(X, y).transform(X)
    
    lda = LinearDiscriminantAnalysis(n_components=nComp)
    X_r2 = lda.fit(X, y).transform(X)
    
    # Percentage of variance explained for each components
    #print('explained variance ratio (first two components): %s'
    #      % str(pca.explained_variance_ratio_))

    l = len(target_names)
    colors = ['darkcyan', 'indigo', 'darkgreen', 'darkgoldenrod', 'darkblue'
        , 'darkorange', 'mediumvioletred', 'crimson', 'darksalmon', 'darkred', 'cyan', 'orange','green']   
    colors = colors[0:l]
    target_codes = list(range(0,l))
    
    plt.figure(figsize=(8, 8))
    lw = 2
    for xComp in range(1,nComp+1):
        for yComp in range(xComp+1,nComp+1):
            for color, i, target_name in zip(colors, target_codes, target_names):
                plt.scatter(X_r[y == i, xComp-1], X_r[y == i, yComp-1], color=color, alpha=.8, lw=lw,
                            label=target_name)
            plt.legend(loc='best', shadow=False, scatterpoints=1)
            plt.title('PCA applied to PSN dataset')
            plt.xlabel('PCA ' + str(xComp))
            plt.ylabel('PCA ' + str(yComp))
            if save:
                plt.savefig('PCA component ' + str(xComp) + ' by ' + str(yComp) + '.png') 
            
            plt.figure()
            for color, i, target_name in zip(colors, target_codes, target_names):
                plt.scatter(X_r2[y == i, xComp-1], X_r2[y == i, yComp-1], alpha=.8, color=color, lw=lw,
                            label=target_name)
            plt.legend(loc='best', shadow=False, scatterpoints=1)
            plt.title('LDA applied to PSN dataset')
            plt.xlabel('LDA ' + str(xComp))
            plt.ylabel('LDA ' + str(yComp))
            if save:
                plt.savefig('LDA component ' + str(xComp) + ' by ' + str(yComp) + '.png') 
            
            plt.show()
###############################################################################
#def CancerTypesDiscAnalysis(dfAllOD, CancerTypes, nComp = 2, save=False):
#    """    
#    We want to analyze how well different cancer types are separated from each other.
#    We filter out samples where 'CancerStatus' = 'Primary solid Tumor' and ClassVar = 'Project'.
#    Then we chose which CancerTypes to compare against each other and draw plots using PCA and LDA
#    for the analysis purposes. 
#    dfAllOD is dataframe of all PSN data
#    CancerTypes is a list of the cancer types that we want to compare against each other. 
#    To be able to see LDA plots, compare a min of 3 cancer types at a time.
#    """    
#    # from CancerStatus keep only 'Primary solid Tumor'
#    ClassVar = 'CancerStatus'
#    toKeep = ['Primary solid Tumor']
#    df_pst = FilterLevels(dfAllOD, ClassVar, toKeep)
#    
#    # Now remove extra variables, we keep only Project 
#    df_pst.drop(['CancerStatus', 'TumorStage', 'Race', 'Gender'], axis=1, inplace=True)
#    
##    # Print counts for missing values.
##    OD.printNaNs(df_pst)
#
#    # drop all the rows where there is any missing data
#    dfdna_pst = dropNaNs(df_pst)
#    
#    # Class variable is the last column, bringing it to the first place.
#    cols = dfdna_pst.columns.tolist()
#    cols = cols[-1:] + cols[:-1]
#    dfdna_pst = dfdna_pst[cols]
#    
#    ClassVar = 'Project'
##    OD.printClassVarValCounts(dfdna_pst,ClassVar)
##    ProjectIDS = OD.returnVarLevels(dfdna_pst, ClassVar)
#    
#    dfdna_pst_fl = FilterLevels(dfdna_pst, ClassVar, CancerTypes)
#    
#    dfdna_pst_fl, ClassVarEncOrder = mapClassVar(dfdna_pst_fl,ClassVar)
#    
#    dfdna_pst_fl_cd = CleanData(dfdna_pst_fl,2)
#        
#    X_scaled_lg, y_lg = fitLogTransform(dfdna_pst_fl_cd)
#    
##    target_names = ClassVarEncOrder
#    plotPCAvsLDA(X_scaled_lg, y_lg, nComp, ClassVarEncOrder, save=save)
#
##    return ClassVarEncOrder
###############################################################################
def dfCancerTypesOrdered(dfAllOD):
    """    
    We want to analyze how well different cancer types are separated from each other.
    We filter out samples where 'CancerStatus' = 'Primary solid Tumor' and ClassVar = 'Project'.
    Then we chose which CancerTypes to compare against each other and draw plots using PCA and LDA
    for the analysis purposes. 
    dfAllOD is dataframe of all PSN data
    CancerTypes is a list of the cancer types that we want to compare against each other. 
    To be able to see LDA plots, compare a min of 3 cancer types at a time.
    """    
    ClassVar = 'CancerStatus'
    #toKeep = ['Primary solid Tumor']
    toKeep = ['Solid Tissue Normal', 'Primary solid Tumor']
    df_pst = FilterLevels(dfAllOD, ClassVar, toKeep)
    
    # Determine the number of genes (features) in the dataframe
    # Note that the dataframe is arranged such that "CancerStatus" is the first
    # class variable, so all columns before "CancerStatus" are genes.
    numFeatures = dfAllOD.columns.get_loc('CancerStatus')
    
    # Now remove extra variables, we keep only Project
    remVars = df_pst.columns[numFeatures:].tolist()
    remVars.remove('Project')
    df_pst.drop(remVars, axis=1, inplace=True)
    
    # drop all the rows where there is any missing data
    dfdna_pst = dropNaNs(df_pst)
    
    # Class variable is the last column, bringing it to the first place.
    cols = dfdna_pst.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    dfdna_pst = dfdna_pst[cols]

    # create a data frame of unique values of all cancer types (Project IDs) sorted 
    # with respect to descending frequency.
    ClassVar = 'Project'
    VarLevels = pd.DataFrame(dfdna_pst[ClassVar].value_counts())
    VarLevels.reset_index(inplace=True)
    VarLevels.rename(columns={'index':ClassVar,ClassVar:'Frequency'}, inplace=True)
    VarLevels.sort_values(by='Frequency', inplace=True, ascending=False)

    # Here we get a list of all uniques values of Project sorted by descending frequency 
    CancerTypesSorted = VarLevels[ClassVar].tolist()

    return dfdna_pst, CancerTypesSorted

###############################################################################
###############################################################################
def returnVarLevelsSorted(dfdna, ClassVar):

    """
    Returns the unique values/levels of the given variable.
    """
    # create a data frame of unique values of all cancer types (Project IDs) sorted 
    # with respect to descending frequency.

    VarLevels = pd.DataFrame(dfdna[ClassVar].value_counts())
    VarLevels.reset_index(inplace=True)
    VarLevels.rename(columns={'index':ClassVar,ClassVar:'Frequency'}, inplace=True)
    VarLevels.sort_values(by='Frequency', inplace=True, ascending=False)

    # Here we get a list of all uniques values of Project sorted by descending frequency 
    VarLevelsSorted = VarLevels[ClassVar].tolist()

    return VarLevels, VarLevelsSorted

    
###############################################################################
###############################################################################


import seaborn as sns
import matplotlib.patheffects as PathEffects

def tSNEscatter(x, colors, ClassVarEncOrder, nClasses):
    # We choose a color palette with seaborn.
    palette = np.array(sns.color_palette("husl", n_colors=nClasses))

    # We create a scatter plot.
    f = plt.figure(figsize=(8, 8))
    ax = plt.subplot(aspect='equal')
    sc = ax.scatter(x[:,0], x[:,1], lw=0, s=40, c=palette[colors.astype(np.int)])
 #   plt.xlim(-25, 25)
 #   plt.ylim(-25, 25)
    ax.axis('off')
    ax.axis('tight')
    plt.title('TSNE 2D projection applied to dataset')

    # We add the labels for each class.
    txts = []
    for i in range(nClasses):
        # Position of each label.
        xtext, ytext = np.median(x[colors == i, :], axis=0)
        name = ClassVarEncOrder[i]
        txt = ax.text(xtext, ytext, name, fontsize=12) #name[5:]
        txt.set_path_effects([
            PathEffects.Stroke(linewidth=5, foreground="w"),
            PathEffects.Normal()])
        txts.append(txt)
    
    plt.show()    
  #  return f, ax, sc, txts
###############################################################################
#from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.feature_selection import RFECV
from itertools import compress

###############################################################################

def RecursiceFeatureElimCV(mod_name, CV, classifier, data,n_splits, scoring):

  col_names = list(data)
  feature_names = col_names[1:]  
  array = data.values
#  X = array[:,1:len(data.columns)]
  y = array[:,0]
  
  X, _ = fitLogTransform(data)

  if CV == 'SKF':
      cv = StratifiedKFold(n_splits=n_splits, shuffle=True)
  elif CV == 'KF':
      cv = KFold(n_splits=n_splits, shuffle=True)


# Create the RFE object and compute a cross-validated score.
#svc = SVC(kernel="linear")
  rfecv = RFECV(estimator=classifier, step=1, cv=cv,
              scoring=scoring)
  rfecv.fit(X, y)

  print("Optimal number of Genes selected: %d" % rfecv.n_features_)

#
#print("Num Features:", fit.n_features_)
  print("Selected Genes:") #, rfecv.support_)
  fil = list(rfecv.support_)
  selected_genes = list(compress(feature_names, fil))
  print(selected_genes)
  #np.invert
  print("\nGenes not selected {0}:".format(len(feature_names)- rfecv.n_features_))
  notselected_genes = list(compress(feature_names, np.invert(fil)))
  print(notselected_genes)
  
# print("Feature Ranking:", rfecv.ranking_)
#
# Plot number of features VS. cross-validation scores
  plt.figure()
  plt.xlabel("Number of genes selected")
  plt.ylabel("Cross validation score")
  plt.title("Selection of Most Important Genes using RFECV (%s) \n Model: %s" % (CV, mod_name)) #, size=13)
  
  plt.plot(range(1, len(rfecv.grid_scores_) + 1), rfecv.grid_scores_, 'bo-')
  plt.grid(True)
  plt.show()
  return selected_genes, notselected_genes#X, y, rfecv.grid_scores_
###############################################################################

#def BRCA_TumorStageMapping(x):
#    if x in ['stage ia','stage ib']:
#        return 'stage i'
#    elif x in ['stage iia','stage iib']:
#        return 'stage ii'
#    elif x in ['stage iiia','stage iiib','stage iiic']:
#        return 'stage iii'
#    else:
#        return x
###############################################################################
        
def BeegleSearchCommonGenes(beegleSearchResults, localGeneSet=False):
    if localGeneSet is False:
        dfGeneNamesMappingPSN = pd.read_csv("dfGeneNamesMappingPSN", sep=",")
        localGeneSet = dfGeneNamesMappingPSN['GeneName'].tolist()
    
    dfBeegleResults = pd.read_table(beegleSearchResults + ".tsv")
    beegleGeneSet = dfBeegleResults['Gene Symbol'].tolist()
    #return the intersection of two lists 
    return list(set(localGeneSet) & set(beegleGeneSet))    
