# script to run DE analysis on all available mutation class variables and
# cancer types, as well as all available tumor stage combinations.


##### SPECIFY LOCATION OF THE CancerProteinSecretionML DIRECTORY #####
main_dir <- "/usr/yourname/Documents/CancerProteinSecretionML"
######################################################################


# source DE function
source(file.path(main_dir, 'scripts', 'CancerDiffExpAnalysis.R'))

# get list of all available cancer types
allCancerTypes <- CancerDiffExpAnalysis(main_dir=main_dir)

# get list of all available mutation class variables
allClassVars <- CancerDiffExpAnalysis(cancerType='all', main_dir=main_dir)
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
    CancerDiffExpAnalysis(cancerType=cancer, classVar=mutVar, classVarLevels=c('FALSE', 'TRUE'), main_dir=main_dir)
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
    CancerDiffExpAnalysis(cancerType=cancer, classVar='TumorStageMerged', classVarLevels=sp, main_dir=main_dir)
  }
}


# evaluate TumorStageBinary for all cancer types
CancerDiffExpAnalysis(cancerType='all', classVar='TumorStageBinary', classVarLevels=c('stage i-iii','stage iv'), main_dir=main_dir)






