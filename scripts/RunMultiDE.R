# script to run DE analysis on all mutation class variables and cancer types, as well as all tumor stage combinations.


source('/Users/jonrob/Documents/PostDoc/CancerProteinSecretionML/scripts/CancerDiffExpAnalysis.R')
setwd("/Users/jonrob/Documents/PostDoc/CancerProteinSecretionML")

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


# loop through all TumorStageMerged combinations
stage_pairs <- list(c('stage i', 'stage ii'), c('stage i','stage iii'), c('stage i','stage iv'), c('stage i','stage x'),
                    c('stage ii','stage iii'), c('stage ii','stage iv'), c('stage ii','stage x'), c('stage iii','stage iv'),
                    c('stage iii','stage x'), c('stage iv','stage x'))
for (cancer in allCancerTypes) {
  for (sp in stage_pairs) {
    # check that the analysis has not already been performed
    if ( any(grepl(paste(cancer, 'TumorStage', gsub(' ','',paste(sp, collapse='_')), 'DEresults', sep='_'), dir(paste('results', cancer, sep='/')))) ) {
      message('Already analyzed ', paste(sp, collapse=' vs. '), ' in cancer type ', cancer, ' - skipping.')
      next
    }
    # run analysis
    message('Analyzing ', paste(sp, collapse=' vs. '), ' in cancer type ', cancer, '.')
    CancerDiffExpAnalysis(cancerType=cancer, classVar='TumorStageMerged', classVarLevels=sp)
  }
}


# evaluate TumorStageBinary for all cancer types
CancerDiffExpAnalysis(cancerType='all', classVar='TumorStageBinary', classVarLevels=c('stage i-iii','stage iv'))






