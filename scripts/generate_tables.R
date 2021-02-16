# Script to generate information for supplementary table information

library(dplyr)
library(tidyr)
library(tibble)

# specify project and output directory
proj_dir <- '/Users/jonrob/Documents/PostDoc/CancerProteinSecretionML'
supp_dir <- file.path(proj_dir, 'doc', 'manuscript', 'REVISION', 'supporting_information')
results_folder <- 'results'


##############################
### Load and organize data ###
##############################

# load PSP expression data
expdata <- readRDS(paste0(proj_dir, '/data/allcancerdata.rds'))
expdata$Project <- sub('TCGA-', '', expdata$Project)  # remove "TCGA-" prefix on cancer type abbrevs

# define scoring function for DE FDR-adjusted p-values
score_pvals <- function(x) {
  y <- -log10(x)
  (y - min(y)) / (max(y) - min(y))
}

# load gene annotations and ML gene scores 
gene_data <- read.delim(paste0(proj_dir, '/data/PSPgenes.txt'), row.names=1)
genes <- rownames(gene_data)
cancers <- dir(file.path(proj_dir, results_folder))
stages <- c('stagei_stageii', 'stagei_stageiii', 'stagei_stageiv', 'stageii_stageiii', 'stageii_stageiv', 'stageiii_stageiv')
stageNames <- c('stage1v2', 'stage1v3', 'stage1v4', 'stage2v3', 'stage2v4', 'stage3v4')
cancer_stages <- unlist(lapply(cancers, function(x) paste(x, stageNames, sep='_')))

score_vars <- c('scores_cancerStatus', 'scores_mutTP53', 'scores_tumorStage', 'scores_stageRegress')
# DE_vars <- c('DE_cancerStatus', 'DE_mutTP53', 'DE_tumorStage')

model_names <- c('ExtraTreesClassifier', 'RandomForestClassifier', 'AdaBoostClassifier', 'XGBClassifier',
                 'LinearDiscriminantAnalysis', 'SVC', 'LassoRegression', 'RidgeRegression', 'Average')
reg_model_names <- c('ExtraTreesRegressor', 'RandomForestRegressor', 'AdaBoostRegressor', 'XGBRegressor',
                     'SVR', 'LassoRegression', 'RidgeRegression', 'Average')

new_model_names <- c('Extra Trees', 'Random Forest', 'AdaBoost', 'XGBoost', 'LDA', 'SVM', 'Lasso', 'Ridge', 'Average')
new_reg_model_names <- c('Extra Trees', 'Random Forest', 'AdaBoost', 'XGBoost', 'SVM', 'Lasso', 'Ridge', 'Average')

DE_vars <- c('DE_log2FC', 'DE_FDR', 'DE_FDRscore')

for (x in score_vars) assign(x, list())
for (model in c(model_names, DE_vars)) {
  scores_cancerStatus[[model]] <- matrix(NA, nrow=length(genes), ncol=length(cancers), dimnames=list(genes, cancers))
  scores_mutTP53[[model]]      <- matrix(NA, nrow=length(genes), ncol=length(cancers), dimnames=list(genes, cancers))
  scores_tumorStage[[model]]   <- matrix(NA, nrow=length(genes), ncol=length(cancer_stages), dimnames=list(genes, cancer_stages))
}
for (reg_model in c(reg_model_names)) {
  scores_stageRegress[[reg_model]] <- matrix(NA, nrow=length(genes), ncol=length(cancers), dimnames=list(genes, cancers))
}

scores_cancerStatus$model_score <- matrix(NA, nrow=length(model_names)-1, ncol=length(cancers), dimnames=list(model_names[1:(length(model_names)-1)], cancers))
scores_mutTP53$model_score      <- matrix(NA, nrow=length(model_names)-1, ncol=length(cancers), dimnames=list(model_names[1:(length(model_names)-1)], cancers))
scores_tumorStage$model_score   <- matrix(NA, nrow=length(model_names)-1, ncol=length(cancer_stages), dimnames=list(model_names[1:(length(model_names)-1)], cancer_stages))
scores_stageRegress$model_score <- matrix(NA, nrow=length(reg_model_names)-1, ncol=length(cancers), dimnames=list(reg_model_names[1:(length(reg_model_names)-1)], cancers))

for (cancer in cancers) {
  cancer_path <- file.path(proj_dir, results_folder, cancer)
  files_cancerStatus <- dir(cancer_path, 'CancerStatus')
  files_mutTP53 <- dir(cancer_path, 'mutTP53')
  files_stageRegress <- dir(cancer_path, 'TumorStage_regression')
  files_tumorStage <- setdiff(dir(cancer_path, 'TumorStage'), files_stageRegress)
  
  if (length(files_cancerStatus) > 0) {
    scores <- read.csv(paste0(cancer_path, '/', files_cancerStatus[grepl('GenesRanking', files_cancerStatus)]), row.names=1)
    for (model in model_names) {
      scores_cancerStatus[[model]][, cancer] <- scores[match(genes, rownames(scores)), model]
    }
    suppressWarnings(try({
      de_res <- read.delim(paste0(cancer_path, '/', files_cancerStatus[grepl('DEresults', files_cancerStatus)]), row.names=1)
      scores_cancerStatus$DE_log2FC[, cancer] <- de_res$logFC[match(genes, rownames(de_res))]
      scores_cancerStatus$DE_FDR[, cancer] <- de_res$FDR[match(genes, rownames(de_res))]
      scores_cancerStatus$DE_FDRscore[, cancer] <- score_pvals(de_res$FDR[match(genes, rownames(de_res))])
    }, silent=T))
    
    model_score <- read.csv(paste0(cancer_path, '/', files_cancerStatus[grepl('CVscores', files_cancerStatus)]), row.names=1)
    scores_cancerStatus$model_score[, cancer] <- model_score$Score[match(model_names[1:(length(model_names)-1)], rownames(model_score))]
  }
  
  if (length(files_mutTP53) > 0) {
    scores <- read.csv(paste0(cancer_path, '/', files_mutTP53[grepl('GenesRanking', files_mutTP53)]), row.names=1)
    for (model in model_names) {
      scores_mutTP53[[model]][, cancer] <- scores[match(genes, rownames(scores)), model]
    }
    suppressWarnings(try({
      de_res <- read.delim(paste0(cancer_path, '/', files_mutTP53[grepl('DEresults', files_mutTP53)]), row.names=1)
      scores_mutTP53$DE_log2FC[, cancer] <- de_res$logFC[match(genes, rownames(de_res))]
      scores_mutTP53$DE_FDR[, cancer] <- de_res$FDR[match(genes, rownames(de_res))]
      scores_mutTP53$DE_FDRscore[, cancer] <- score_pvals(de_res$FDR[match(genes, rownames(de_res))])
    }, silent=T))
    
    model_score <- read.csv(paste0(cancer_path, '/', files_mutTP53[grepl('CVscores', files_mutTP53)]), row.names=1)
    scores_mutTP53$model_score[, cancer] <- model_score$Score[match(model_names[1:(length(model_names)-1)], rownames(model_score))]
  }
  
  if (length(files_tumorStage) > 0) {
    files_tumorStage_geneScores <- files_tumorStage[grepl('GenesRanking', files_tumorStage)]
    for (f in files_tumorStage_geneScores) {
      scores <- read.csv(paste0(cancer_path, '/', f), row.names=1)
      f_stage <- paste(unlist(strsplit(f, '_'))[3:4], collapse='_')
      if (!(f_stage %in% stages)) {
        next
      }
      f_name <- paste(cancer, stageNames[stages %in% f_stage], sep='_')
      for (model in model_names) {
        scores_tumorStage[[model]][, f_name] <- scores[match(genes, rownames(scores)), model]
      }
      suppressWarnings(try({
        de_res <- read.delim(paste0(cancer_path, '/', sub('GenesRanking.csv', 'DEresults.txt', f)), row.names=1)
        scores_tumorStage$DE_log2FC[, f_name] <- de_res$logFC[match(genes, rownames(de_res))]
        scores_tumorStage$DE_FDR[, f_name] <- de_res$FDR[match(genes, rownames(de_res))]
        scores_tumorStage$DE_FDRscore[, f_name] <- score_pvals(de_res$FDR[match(genes, rownames(de_res))])
      }, silent=T))
      
      model_score <- read.csv(paste0(cancer_path, '/', sub('GenesRanking', 'CVscores', f)), row.names=1)
      scores_tumorStage$model_score[, f_name] <- model_score$Score[match(model_names[1:(length(model_names)-1)], rownames(model_score))]
    }
  }
  
  if (length(files_stageRegress) > 0) {
    scores <- read.csv(paste0(cancer_path, '/', files_stageRegress[grepl('GenesRanking', files_stageRegress)]), row.names=1)
    for (reg_model in reg_model_names) {
      scores_stageRegress[[reg_model]][, cancer] <- scores[match(genes, rownames(scores)), reg_model]
    }
    model_score <- read.csv(paste0(cancer_path, '/', files_stageRegress[grepl('CVscores', files_stageRegress)]), row.names=1)
    scores_stageRegress$model_score[, cancer] <- model_score$Score[match(reg_model_names[1:(length(reg_model_names)-1)], rownames(model_score))]
  }
  
}

# remove cancer types that were not included in each analysis
keep_cancers <- colnames(scores_cancerStatus$Average)[colSums(is.na(scores_cancerStatus$Average)) < nrow(scores_cancerStatus$Average)]
for (item in names(scores_cancerStatus)) {
  scores_cancerStatus[[item]] <- scores_cancerStatus[[item]][, keep_cancers] %>% replace_na(0)
}
keep_cancers <- colnames(scores_mutTP53$Average)[colSums(is.na(scores_mutTP53$Average)) < nrow(scores_mutTP53$Average)]
for (item in names(scores_mutTP53)) {
  scores_mutTP53[[item]] <- scores_mutTP53[[item]][, keep_cancers] %>% replace_na(0)
}
keep_cancers <- colnames(scores_tumorStage$Average)[colSums(is.na(scores_tumorStage$Average)) < nrow(scores_tumorStage$Average)]
for (item in names(scores_tumorStage)) {
  scores_tumorStage[[item]] <- scores_tumorStage[[item]][, keep_cancers] %>% replace_na(0)
}
keep_cancers <- colnames(scores_stageRegress$Average)[colSums(is.na(scores_stageRegress$Average)) < nrow(scores_stageRegress$Average)]
for (item in names(scores_stageRegress)) {
  scores_stageRegress[[item]] <- scores_stageRegress[[item]][, keep_cancers] %>% replace_na(0)
}

# rename models in model_score slots
rownames(scores_cancerStatus$model_score) <- new_model_names[match(rownames(scores_cancerStatus$model_score), model_names)]
rownames(scores_mutTP53$model_score) <- new_model_names[match(rownames(scores_mutTP53$model_score), model_names)]
rownames(scores_tumorStage$model_score) <- new_model_names[match(rownames(scores_tumorStage$model_score), model_names)]
rownames(scores_stageRegress$model_score) <- new_reg_model_names[match(rownames(scores_stageRegress$model_score), reg_model_names)]

# combine scores into list and remove intermediate variables
allscores <- list(cancerStatus=scores_cancerStatus, mutTP53=scores_mutTP53, tumorStage=scores_tumorStage, stageRegress=scores_stageRegress)
rm(list=setdiff(ls(), c('proj_dir', 'expdata', 'allscores', 'supp_dir', 'gene_data')))
invisible(gc())



#############################################
### Collect details of expression dataset ###
#############################################
keep <- !is.na(expdata$Project)
expdata <- expdata[keep, ]
cancers <- unique(expdata$Project)
cancer_name_data <- list(abbrev=c('ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC',
                                  'KIRP','LAML','LGG','LIHC','LUAD','LUSC','MESO','OV','PAAD','PCPG','PRAD','READ',
                                  'SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC','UCS','UVM'),
                         name=c('Adrenocortical Carcinoma','Bladder Urothelial Carcinoma','Breast Invasive Carcinoma',
                                'Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma','Cholangiocarcinoma',
                                'Colon Adenocarcinoma','Lymphoid Neoplasm Diffuse Large B-cell Lymphoma',
                                'Esophageal Carcinoma','Glioblastoma Multiforme','Head and Neck Squamous Cell Carcinoma',
                                'Kidney Chromophobe','Kidney Renal Clear Cell Carcinoma',
                                'Kidney Renal Papillary Cell Carcinoma','Acute Myeloid Leukemia','Brain Lower Grade Glioma',
                                'Liver Hepatocellular Carcinoma','Lung Adenocarcinoma','Lung Squamous Cell Carcinoma',
                                'Mesothelioma','Ovarian Serous Cystadenocarcinoma','Pancreatic Adenocarcinoma',
                                'Pheochromocytoma and Paraganglioma','Prostate Adenocarcinoma','Rectum Adenocarcinoma',
                                'Sarcoma','Skin Cutaneous Melanoma','Stomach Adenocarcinoma','Testicular Germ Cell Tumors',
                                'Thyroid Carcinoma','Thymoma','Uterine Corpus Endometrial Carcinoma','Uterine Carcinosarcoma',
                                'Uveal Melanoma'))

class_vars <- c('CancerStatus', 'mutTP53', 'TumorStageMerged')
class_levels <- list(CancerStatus=c('Solid Tissue Normal', 'Primary solid Tumor'),
                     mutTP53=c('FALSE', 'TRUE'),
                     TumorStageMerged=c('stage i', 'stage ii', 'stage iii', 'stage iv'))

metadat <- as.data.frame(matrix(0, nrow=length(cancers), ncol=length(unlist(class_levels)), dimnames=list(cancers, unlist(class_levels))))
for (cvar in class_vars) {
  dat <- expdata[, c('Project', cvar)]
  for (clev in class_levels[[cvar]]) {
    metadat[, clev] <- unlist(lapply(cancers, function(x) sum(dat[dat$Project %in% x, cvar] == clev)))
  }
}

cancer_names <- cancer_name_data$name[match(cancers, cancer_name_data$abbrev)]
metadat <- cbind(cancer_names, metadat)

write.table(metadat %>% rownames_to_column('Abbrev.'), file=file.path(supp_dir, 'sample_metadata.txt'), sep='\t', row.names=F)


###########################################################
### Export data to text file (for supplementary tables) ###
###########################################################
write.table(as.data.frame(allscores$cancerStatus$Average) %>% rownames_to_column('Gene'),
            file=file.path(supp_dir, 'ConsensusMLGeneScores_CancerStatus.txt'), sep='\t', row.names=F)
write.table(as.data.frame(allscores$mutTP53$Average) %>% rownames_to_column('Gene'),
            file=file.path(supp_dir, 'ConsensusMLGeneScores_mutTP53.txt'), sep='\t', row.names=F)
write.table(as.data.frame(allscores$tumorStage$Average) %>% rownames_to_column('Gene'),
            file=file.path(supp_dir, 'ConsensusMLGeneScores_tumorStage.txt'), sep='\t', row.names=F)
write.table(as.data.frame(allscores$stageRegress$Average) %>% rownames_to_column('Gene'),
            file=file.path(supp_dir, 'ConsensusMLGeneScores_StageRegress.txt'), sep='\t', row.names=F)

write.table(as.data.frame(allscores$cancerStatus$DE_log2FC) %>% rownames_to_column('Gene'),
            file=file.path(supp_dir, 'DE_log2FC_CancerStatus.txt'), sep='\t', row.names=F)
write.table(as.data.frame(allscores$mutTP53$DE_log2FC) %>% rownames_to_column('Gene'),
            file=file.path(supp_dir, 'DE_log2FC_mutTP53.txt'), sep='\t', row.names=F)
write.table(as.data.frame(allscores$tumorStage$DE_log2FC) %>% rownames_to_column('Gene'),
            file=file.path(supp_dir, 'DE_log2FC_tumorStage.txt'), sep='\t', row.names=F)

write.table(as.data.frame(allscores$cancerStatus$DE_FDR) %>% rownames_to_column('Gene'),
            file=file.path(supp_dir, 'DE_FDR_CancerStatus.txt'), sep='\t', row.names=F)
write.table(as.data.frame(allscores$mutTP53$DE_FDR) %>% rownames_to_column('Gene'),
            file=file.path(supp_dir, 'DE_FDR_mutTP53.txt'), sep='\t', row.names=F)
write.table(as.data.frame(allscores$tumorStage$DE_FDR) %>% rownames_to_column('Gene'),
            file=file.path(supp_dir, 'DE_FDR_tumorStage.txt'), sep='\t', row.names=F)

write.table(as.data.frame(allscores$cancerStatus$DE_FDRscore) %>% rownames_to_column('Gene'),
            file=file.path(supp_dir, 'DE_FDRscore_CancerStatus.txt'), sep='\t', row.names=F)
write.table(as.data.frame(allscores$mutTP53$DE_FDRscore) %>% rownames_to_column('Gene'),
            file=file.path(supp_dir, 'DE_FDRscore_mutTP53.txt'), sep='\t', row.names=F)
write.table(as.data.frame(allscores$tumorStage$DE_FDRscore) %>% rownames_to_column('Gene'),
            file=file.path(supp_dir, 'DE_FDRscore_tumorStage.txt'), sep='\t', row.names=F)

write.table(as.data.frame(allscores$cancerStatus$model_score) %>% rownames_to_column('ML algorithm'),
            file=file.path(supp_dir, 'ROCAUC_CancerStatus.txt'), sep='\t', row.names=F)
write.table(as.data.frame(allscores$mutTP53$model_score) %>% rownames_to_column('ML algorithm'),
            file=file.path(supp_dir, 'ROCAUC_mutTP53.txt'), sep='\t', row.names=F)
write.table(as.data.frame(allscores$tumorStage$model_score) %>% rownames_to_column('ML algorithm'),
            file=file.path(supp_dir, 'ROCAUC_tumorStage.txt'), sep='\t', row.names=F)
write.table(as.data.frame(allscores$stageRegress$model_score) %>% rownames_to_column('ML algorithm'),
            file=file.path(supp_dir, 'NegMSE_stageRegress.txt'), sep='\t', row.names=F)

