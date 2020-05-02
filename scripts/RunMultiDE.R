# script to run DE analysis on all mutation class variables and cancer types


source('~/Documents/PostDoc/CancerOmicsDataExploration/CancerDiffExpAnalysis.R')
setwd("~/Documents/PostDoc/CancerOmicsDataExploration")

# get list of all available cancer types
allCancerTypes <- CancerDiffExpAnalysis()

# get list of all available mutation class variables
allClassVars <- CancerDiffExpAnalysis(cancerType='all')
mutClassVars <- allClassVars[substr(allClassVars,1,3) == 'mut']

# loop through all cancer types and mutation variables
for (cancer in allCancerTypes) {
  for (mutVar in mutClassVars) {
    # check that the analysis has not already been performed
    if ( any(grepl(paste(cancer, mutVar, 'DEresults', sep='_'), dir(paste('results', cancer, sep='/')))) ) {
      message('Already analyzed ', mutVar, ' in cancer type ', cancer, ' - skipping.')
      next
    }
    # run analysis
    message('Analyzing ', mutVar, ' in cancer type ', cancer, '.')
    CancerDiffExpAnalysis(cancerType=cancer, classVar=mutVar, classVarLevels=c('FALSE', 'TRUE'))
  }
}




