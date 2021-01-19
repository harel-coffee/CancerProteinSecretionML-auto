
# load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(viridisLite)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(cowplot)


# specify main project directory and the directory where to save figures
proj_dir <- '/Users/jonrob/Documents/PostDoc/CancerProteinSecretionML'
fig_dir <- paste0(proj_dir, '/doc/manuscript/figures/fig_pieces')


##############################
### Load and organize data ###
##############################

# load PSP expression data
expdata <- readRDS(paste0(proj_dir, '/data/allcancerdata_psp.rds'))
expdata$Project <- sub('TCGA-', '', expdata$Project)  # remove "TCGA-" prefix on cancer type abbrevs

# define scoring function for DE FDR-adjusted p-values
score_pvals <- function(x) {
  y <- -log10(x)
  (y - min(y)) / (max(y) - min(y))
}

# load gene annotations and ML gene scores 
gene_data <- read.delim(paste0(proj_dir, '/data/PSPgenes.txt'), row.names=1)
genes <- rownames(gene_data)
cancers <- dir(paste0(proj_dir, '/results'))
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

scores_cancerStatus$roc_auc <- matrix(NA, nrow=length(model_names)-1, ncol=length(cancers), dimnames=list(model_names[1:(length(model_names)-1)], cancers))
scores_mutTP53$roc_auc      <- matrix(NA, nrow=length(model_names)-1, ncol=length(cancers), dimnames=list(model_names[1:(length(model_names)-1)], cancers))
scores_tumorStage$roc_auc   <- matrix(NA, nrow=length(model_names)-1, ncol=length(cancer_stages), dimnames=list(model_names[1:(length(model_names)-1)], cancer_stages))
scores_stageRegress$neg_mse <- matrix(NA, nrow=length(reg_model_names)-1, ncol=length(cancers), dimnames=list(reg_model_names[1:(length(reg_model_names)-1)], cancers))

for (cancer in cancers) {
  cancer_path <- paste0(proj_dir, '/results/', cancer)
  files_cancerStatus <- dir(cancer_path, 'CancerStatus')
  files_mutTP53 <- dir(cancer_path, 'mutTP53')
  files_stageRegress <- dir(cancer_path, 'TumorStage_regression')
  files_tumorStage <- setdiff(dir(cancer_path, 'TumorStage'), files_stageRegress)
  
  if (length(files_cancerStatus) > 0) {
    scores <- read.csv(paste0(cancer_path, '/', files_cancerStatus[grepl('GenesRanking', files_cancerStatus)]), row.names=1)
    for (model in model_names) {
      scores_cancerStatus[[model]][, cancer] <- scores[match(genes, rownames(scores)), model]
    }
    de_res <- read.delim(paste0(cancer_path, '/', files_cancerStatus[grepl('DEresults', files_cancerStatus)]), row.names=1)
    scores_cancerStatus$DE_log2FC[, cancer] <- de_res$logFC[match(genes, rownames(de_res))]
    scores_cancerStatus$DE_FDR[, cancer] <- de_res$FDR[match(genes, rownames(de_res))]
    scores_cancerStatus$DE_FDRscore[, cancer] <- score_pvals(de_res$FDR[match(genes, rownames(de_res))])
    
    roc_auc <- read.csv(paste0(cancer_path, '/', files_cancerStatus[grepl('CVscores', files_cancerStatus)]), row.names=1)
    scores_cancerStatus$roc_auc[, cancer] <- roc_auc$Score[match(model_names[1:(length(model_names)-1)], rownames(roc_auc))]
  }
  
  if (length(files_mutTP53) > 0) {
    scores <- read.csv(paste0(cancer_path, '/', files_mutTP53[grepl('GenesRanking', files_mutTP53)]), row.names=1)
    for (model in model_names) {
      scores_mutTP53[[model]][, cancer] <- scores[match(genes, rownames(scores)), model]
    }
    de_res <- read.delim(paste0(cancer_path, '/', files_mutTP53[grepl('DEresults', files_mutTP53)]), row.names=1)
    scores_mutTP53$DE_log2FC[, cancer] <- de_res$logFC[match(genes, rownames(de_res))]
    scores_mutTP53$DE_FDR[, cancer] <- de_res$FDR[match(genes, rownames(de_res))]
    scores_mutTP53$DE_FDRscore[, cancer] <- score_pvals(de_res$FDR[match(genes, rownames(de_res))])
    
    roc_auc <- read.csv(paste0(cancer_path, '/', files_mutTP53[grepl('CVscores', files_mutTP53)]), row.names=1)
    scores_mutTP53$roc_auc[, cancer] <- roc_auc$Score[match(model_names[1:(length(model_names)-1)], rownames(roc_auc))]
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
      de_res <- read.delim(paste0(cancer_path, '/', sub('GenesRanking.csv', 'DEresults.txt', f)), row.names=1)
      scores_tumorStage$DE_log2FC[, f_name] <- de_res$logFC[match(genes, rownames(de_res))]
      scores_tumorStage$DE_FDR[, f_name] <- de_res$FDR[match(genes, rownames(de_res))]
      scores_tumorStage$DE_FDRscore[, f_name] <- score_pvals(de_res$FDR[match(genes, rownames(de_res))])
      
      roc_auc <- read.csv(paste0(cancer_path, '/', sub('GenesRanking', 'CVscores', f)), row.names=1)
      scores_tumorStage$roc_auc[, f_name] <- roc_auc$Score[match(model_names[1:(length(model_names)-1)], rownames(roc_auc))]
    }
  }
  
  if (length(files_stageRegress) > 0) {
    scores <- read.csv(paste0(cancer_path, '/', files_stageRegress[grepl('GenesRanking', files_stageRegress)]), row.names=1)
    for (reg_model in reg_model_names) {
      scores_stageRegress[[reg_model]][, cancer] <- scores[match(genes, rownames(scores)), reg_model]
    }
    neg_mse <- read.csv(paste0(cancer_path, '/', files_stageRegress[grepl('CVscores', files_stageRegress)]), row.names=1)
    scores_stageRegress$neg_mse[, cancer] <- neg_mse$Score[match(reg_model_names[1:(length(reg_model_names)-1)], rownames(neg_mse))]
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

# rename models in roc_auc and neg_mse slots
rownames(scores_cancerStatus$roc_auc) <- new_model_names[match(rownames(scores_cancerStatus$roc_auc), model_names)]
rownames(scores_mutTP53$roc_auc) <- new_model_names[match(rownames(scores_mutTP53$roc_auc), model_names)]
rownames(scores_tumorStage$roc_auc) <- new_model_names[match(rownames(scores_tumorStage$roc_auc), model_names)]
rownames(scores_stageRegress$neg_mse) <- new_reg_model_names[match(rownames(scores_stageRegress$neg_mse), reg_model_names)]

# combine scores into list and remove intermediate variables
allscores <- list(cancerStatus=scores_cancerStatus, mutTP53=scores_mutTP53, tumorStage=scores_tumorStage, stageRegress=scores_stageRegress)
rm(list=setdiff(ls(), c('proj_dir', 'expdata', 'allscores', 'fig_dir', 'gene_data')))
invisible(gc())



####################################################################
### Fig. 1, toy heatmap demonstrating consensus scoring approach ###
####################################################################

dat <- as.data.frame(matrix(runif(8*8)^2, nrow=8, ncol=8))
dat$Average <- apply(dat, 1, mean)
dat <- dat[order(dat$Average, decreasing=T), ]
pdf(file=paste0(fig_dir, '/toy_heatmap.pdf'), width=1.8, height=1.5, onefile=F)
pheatmap(dat,
         scale='none',
         color=viridis(100),
         cluster_rows=F,
         cluster_cols=F,
         breaks=seq(0,1,len=100))
grid_linewidth <- 1.5
grid_color <- 'white'
grid.ls(grid.force())
grid.gedit('matrix::GRID.rect', gp=gpar(col=grid_color, lwd=grid_linewidth))
invisible(dev.off())


#####################################################
### Fig. X, Panel X: Boxplot of top-scoring genes ###
#####################################################

# specify parameters
classVar <- 'stageRegress'  # 'mutTP53', 'cancerStatus', 'tumorStage', or 'stageRegress'
n_genes <- 10  # specify number of genes to include
sort_by <- 'mean'  # 'mean' or 'top each'
model <- 'Average'  # e.g., 'Average' or 'DE_FDRscore'

# prepare data
scores <- allscores[[classVar]]
dat <- as.data.frame(scores[[model]])
if (sort_by == 'mean') {
  o <- order(apply(dat, 1, mean), decreasing=T)
  top_genes <- rownames(dat)[o[1:n_genes]]
} else if (sort_by == 'top each') {
  o <- apply(dat, 2, function(x) order(x, decreasing=T))
  top_genes <- rownames(dat)[unique(as.vector(o[1:n_genes, ]))]
}
if (classVar == 'mutTP53') {
  dat$targets <- as.factor(rownames(dat) %in% c('BAX', 'HSPA4L', 'KIF23'))
} else {
  dat$targets <- 'NA'
}
dat <- dat[top_genes, ] %>% rownames_to_column('Gene')
dat$Gene <- factor(dat$Gene, levels=rev(top_genes), ordered=T)
dat <- pivot_longer(dat, cols=colnames(dat)[2:(ncol(dat)-1)], values_to='Score', names_to='Cancer')

# generate plot
pdf(file=paste0(fig_dir, '/' , classVar, '_', model, '_boxplot.pdf'), width=4, height=3)
ggplot(dat, aes(x=Score, y=Gene, fill=targets)) +
  geom_boxplot() +
  # geom_violin(trim=T, scale='area', draw_quantiles=0.5) + 
  # geom_jitter(height = 0, width = 0.1) +
  theme_minimal() +
  theme(legend.position='none',
        axis.text=element_text(color='black', size=12),
        axis.title=element_text(size=12)) + 
  xlab('Consensus gene score') +
  scale_fill_manual(values=brewer.pal(3,'Paired')[1:2]) +
  coord_cartesian(xlim=c(0,1))
invisible(dev.off())


#####################################################
### Fig. X, Panel X: Heatmap of top-scoring genes ###
#####################################################

# specify parameters
classVar <- 'stageRegress'  # 'mutTP53', 'cancerStatus', 'tumorStage', or 'stageRegress'
group_cancers <- T  # if 'tumorStage', stages will be grouped (averaged) together by cancer type
n_genes <- 10  # specify number of genes to include
sort_by <- 'mean'  # 'mean' or 'top each'
model <- 'Average'  # e.g., 'Average' or 'DE_FDRscore'
cancer_types <- NULL #c('ACC','KIRP','KIRC','THCA','TGCT')  # use NULL to keep all available cancer types

# prepare data
scores <- allscores[[classVar]]
dat <- as.data.frame(scores[[model]])
if (!is.null(cancer_types)) {
  keep <- unlist(lapply(colnames(dat), function(x) any(startsWith(x, cancer_types))))
  dat <- dat[, keep]
}
if (classVar == 'tumorStage' && group_cancers == T) {
  # merge (average) cancer types
  cancers <- unlist(lapply(colnames(dat), function(x) head(unlist(strsplit(x, '_')), 1)))
  uniq_cancers <- unique(cancers)
  for (cancer in uniq_cancers) {
    dat[, cancer] <- apply(dat[, cancers %in% cancer, drop=F], 1, mean)
  }
  dat <- dat[, uniq_cancers]
}
if (sort_by == 'mean') {
  o <- order(apply(dat, 1, mean), decreasing=T)
  top_genes <- rownames(dat)[o[1:n_genes]]
} else if (sort_by == 'top each') {
  o <- apply(dat, 2, function(x) order(x, decreasing=T))
  top_genes <- rownames(dat)[unique(as.vector(o[1:n_genes, ]))]
}
dat <- dat[top_genes, ]

if (classVar == 'tumorStage' && group_cancers == F) {
  colnames(dat) <- sub('_', ' (', colnames(dat))
  colnames(dat) <- paste0(sub('stage', 'stage ', colnames(dat)), ')')
  plot_height=3.5
} else {
  plot_height=3
}

# generate heatmap
pdf(file=paste0(fig_dir, '/', classVar, '_', model, '_heatmap.pdf'), width=3+2.3*ncol(dat)/22, height=plot_height, onefile=F)
pheatmap(dat,
         scale='none',
         color=magma(100),
         clustering_distance_rows='euclidean',
         clustering_distance_cols='euclidean',
         clustering_method='complete',  # ward.D, ward.D2, single, complete, average, mcquitty, median, centroid. Use "complete" for tumor stages, and "ward.D2" otherwise.
         breaks=seq(0,0.8,len=100),
         angle_col=90,
         border_color='black')
invisible(dev.off())


###################################################################################
### Fig. X, Panel X: Heatmap of top-scoring genes, custom row & column ordering ###
###################################################################################

# specify parameters
classVar <- 'stageRegress'  # 'mutTP53', 'cancerStatus', 'tumorStage', or 'stageRegress'
n_genes <- 5  # specify number of genes to include
model <- 'Average'  # e.g., 'Average' or 'DE_FDRscore'
cancers <- NULL #c('STAD', 'READ', 'COAD', 'KICH', 'THCA')

# prepare data
scores <- allscores[[classVar]]
dat <- as.data.frame(scores[[model]])
if (!is.null(cancers)) {
  dat <- dat[, colnames(dat) %in% cancers]
}
o <- apply(dat, 2, function(x) order(x, decreasing=T))
top_genes <- rownames(dat)[unique(as.vector(o[1:n_genes, ]))]
dat <- dat[top_genes, ]

# cluster columns based on correlation distance
# col_order <- hclust(as.dist(cor(dat, method='pearson')), method='average')$order
col_order <- hclust(dist(t(dat), method='euclidean'), method='ward.D2')$order
dat <- dat[, col_order]
row_order <- NULL
for (i in seq(ncol(dat))) {
  ro <- order(dat[, i], decreasing=T)
  row_order <- c(row_order, setdiff(ro[1:n_genes], row_order))
}
dat <- dat[row_order, ]

# specify annotation colors
# annColors <- list(module = brewer.pal(4, 'Set1') %>% setNames(levels(gene_data$module)),
#                   subsystem = c(brewer.pal(12, 'Set3'), '#969696') %>% setNames(levels(gene_data$subsystem)))
unique_modules <- unique(gene_data$module[row_order])
annColors <- list(Function = viridis(4) %>%
                    setNames(c('Capacity control', 'Folding', 'Trafficking', 'Glycosylation')))
annColors$Function <- annColors$Function[intersect(names(annColors$Function), unique_modules)]
annData <- gene_data[, c('module'), drop=F] %>% setNames('Function')

# generate heatmap
pdf(file=paste0(fig_dir, '/', classVar, '_', model, '_topEach_heatmap.pdf'), width=3+2.3*ncol(dat)/22, height=4, onefile=F)
pheatmap(dat,
         scale='none',
         color=magma(100),
         cluster_rows=F,
         cluster_cols=F,
         breaks=seq(0,0.6,len=100),
         angle_col=90,
         annotation_row=annData,
         annotation_colors=annColors)
grid_linewidth <- 1
grid_color <- 'black'
grid.ls(grid.force())
grid.gedit('matrix::GRID.rect', gp=gpar(col=grid_color, lwd=grid_linewidth))
grid.gedit('row_annotation', gp=gpar(col=grid_color, lwd=grid_linewidth))
grid.gedit('annotation_legend::GRID.rect', gp=gpar(col=grid_color, lwd=grid_linewidth))
invisible(dev.off())


#################################################
### Fig. X, Panel X: Histogram of gene scores ###
#################################################

# specify parameters
classVar <- 'stageRegress'  # 'mutTP53', 'cancerStatus', 'tumorStage', or 'stageRegress'
model <- 'Average'  # e.g., 'Average' or 'DE_FDRscore'

# prepare data
scores <- allscores[[classVar]]
dat <- as.data.frame(apply(scores[[model]], 1, mean)) %>% setNames('Mean score')

# generate plot
pdf(file=paste0(fig_dir, '/', classVar, '_', model, '_histogram.pdf'), width=4, height=3)
ggplot(dat, aes(x=`Mean score`)) +
  geom_histogram(color='black', fill=brewer.pal(3,'Paired')[1], breaks=seq(0, 0.5, length.out=50)) +
  theme_minimal() +
  theme(axis.text=element_text(color='black', size=12),
        axis.title=element_text(size=12)) + 
  ylab('Count') +
  coord_cartesian(xlim=c(0,0.45))
invisible(dev.off())
 

######################################################################
### Fig. X, Panel X: Combined histogram and boxplot of gene scores ###
######################################################################

# specify parameters
classVar <- 'tumorStage'  # 'mutTP53', 'cancerStatus', or 'tumorStage'
model <- 'DE_FDRscore'  # e.g., 'Average' or 'DE_FDRscore'
n_genes <- 10  # specify number of genes to include
cancer_types <- c('ACC','KIRP','KIRC','THCA','TGCT')  # use NULL to keep all available cancer types

# prepare data
scores <- allscores[[classVar]]
dat <- scores[[model]]
if (!is.null(cancer_types)) {
  keep <- unlist(lapply(colnames(dat), function(x) any(startsWith(x, cancer_types))))
  dat <- dat[, keep]
}
dat <- as.data.frame(apply(dat, 1, mean)) %>% setNames('Mean score')

if (model == 'Average') {
  hist_xlab <- 'Mean gene ML score'
  box_xlab <- 'Gene ML score'
} else {
  hist_xlab <- 'Mean gene DE score'
  box_xlab <- 'Gene DE score'
}

# generate histogram
p_hist <- ggplot(dat, aes(x=`Mean score`)) +
  geom_histogram(color='black', fill=brewer.pal(3,'Paired')[1], breaks=seq(0, 0.5, length.out=50)) +
  theme_minimal() +
  theme(axis.text=element_text(color='black', size=12),
        axis.title=element_text(size=12)) + 
  coord_cartesian(xlim=c(0,0.45)) +
  ylab('Genes') +
  xlab(hist_xlab)

# prepare data
dat <- as.data.frame(scores[[model]])
if (!is.null(cancer_types)) {
  keep <- unlist(lapply(colnames(dat), function(x) any(startsWith(x, cancer_types))))
  dat <- dat[, keep]
}
o <- order(apply(dat, 1, mean), decreasing=T)
top_genes <- rownames(dat)[o[1:n_genes]]
if (classVar == 'mutTP53') {
  dat$targets <- as.factor(rownames(dat) %in% c('BAX', 'HSPA4L', 'KIF23'))
} else {
  dat$targets <- 'NA'
}
dat <- dat[top_genes, ] %>% rownames_to_column('Gene')
dat$Gene <- factor(dat$Gene, levels=rev(top_genes), ordered=T)
dat <- pivot_longer(dat, cols=colnames(dat)[2:(ncol(dat)-1)], values_to='Score', names_to='Cancer')

# generate boxplots
p_box <- ggplot(dat, aes(x=Score, y=Gene, fill=targets)) +
  geom_boxplot(color='black') +
  theme_minimal() +
  theme(legend.position='none',
        axis.text=element_text(color='black', size=12),
        axis.title=element_text(size=12)) + 
  xlab(box_xlab) +
  scale_fill_manual(values=brewer.pal(3,'Paired')[1:2]) +
  coord_cartesian(xlim=c(0,1))

# generate combined plot
pdf(file=paste0(fig_dir, '/', classVar, '_', model, '_combined_HistBox.pdf'), width=7, height=3.5)
plot_grid(p_hist, p_box, ncol=2, rel_widths=c(1,1), align='h', axis='tb', scale = 0.98)
invisible(dev.off())


##################################################
### Fig. X, Panel X: Boxplot of ROC AUC values ###
##################################################

# specify parameters
classVar <- 'tumorStage'  # 'mutTP53', 'cancerStatus', or 'tumorStage'
groupby <- 'Cancer'  # 'Cancer' or 'Model'

# prepare data
scores <- allscores[[classVar]]
dat <- as.data.frame(scores$roc_auc)

if (classVar == 'tumorStage') {
  colnames(dat) <- sub('_', ' (', colnames(dat))
  colnames(dat) <- paste0(sub('stage', 'stage ', colnames(dat)), ')')
}

o_model <- rownames(dat)[order(apply(dat, 1, median), decreasing=T)]
o_cancer <- colnames(dat)[order(apply(dat, 2, median), decreasing=T)]
dat <- dat[o_model, o_cancer] %>% rownames_to_column('Model')
dat <- pivot_longer(dat, cols=colnames(dat)[-1], values_to='ROC AUC', names_to='Cancer')
dat$Model <- factor(dat$Model, levels=o_model, ordered=T)
dat$Cancer <- factor(dat$Cancer, levels=o_cancer, ordered=T)

# generate plot
if (groupby == 'Cancer') {
  if (classVar == 'tumorStage') {
    plot_width <- 15
  } else {
    plot_width <- 6
  }
  plot_height <- 4
  hues <- c(0,100)
  text_angle <- 90
  vjust_val <- 0.5
} else if (groupby == 'Model') {
  plot_width <- 3
  plot_height <- 4.5
  hues <- c(230,330)
  text_angle <- 45
  vjust_val <- 1
}
pdf(file=paste0(fig_dir, '/', classVar, '_ROCAUC_', groupby, '_boxplot.pdf'), width=plot_width, height=plot_height)
ggplot(dat, aes_string(x=groupby, y="`ROC AUC`", fill=groupby)) +
  geom_boxplot() +
  # geom_violin(draw_quantiles=0.5, scale='width') +
  theme_minimal() +
  theme(legend.position='none',
        axis.text=element_text(color='black', size=12),
        axis.text.x=element_text(angle=text_angle, hjust=1, vjust=vjust_val),
        axis.title=element_text(size=12)) + 
  xlab(groupby) +
  scale_fill_hue(h=hues, c=75)
  # coord_cartesian(ylim=c(0,1))
# geom_violin(trim=F) + 
invisible(dev.off())


##################################################
### Fig. X, Panel X: Boxplot of NEG MSE values ###
##################################################

# specify parameters
classVar <- 'stageRegress'  # only relevant for 'stageRegress'
groupby <- 'Model'  # 'Cancer' or 'Model'

# prepare data
scores <- allscores[[classVar]]
dat <- as.data.frame(scores$neg_mse)

o_model <- rownames(dat)[order(apply(dat, 1, median), decreasing=T)]
o_cancer <- colnames(dat)[order(apply(dat, 2, median), decreasing=T)]
dat <- dat[o_model, o_cancer] %>% rownames_to_column('Model')
dat <- pivot_longer(dat, cols=colnames(dat)[-1], values_to='NEG MSE', names_to='Cancer')
dat$Model <- factor(dat$Model, levels=o_model, ordered=T)
dat$Cancer <- factor(dat$Cancer, levels=o_cancer, ordered=T)

# generate plot
if (groupby == 'Cancer') {
  plot_width <- 6
  plot_height <- 4
  hues <- c(0,100)
  text_angle <- 90
  vjust_val <- 0.5
} else if (groupby == 'Model') {
  plot_width <- 3
  plot_height <- 4.5
  hues <- c(230,330)
  text_angle <- 45
  vjust_val <- 1
}

pdf(file=paste0(fig_dir, '/', classVar, '_NEGMSE_', groupby, '_boxplot.pdf'), width=plot_width, height=plot_height)
ggplot(dat, aes_string(x=groupby, y="`NEG MSE`", fill=groupby)) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position='none',
        axis.text=element_text(color='black', size=12),
        axis.text.x=element_text(angle=text_angle, hjust=1, vjust=vjust_val),
        axis.title=element_text(size=12)) + 
  xlab(groupby) +
  scale_fill_hue(h=hues, c=75)
invisible(dev.off())


#####################################################################################
### Fig. X, Panel X: Boxplot of tumorStage ROC AUC values, grouped by cancer type ###
#####################################################################################

# specify parameters
classVar <- 'tumorStage'
groupby <- 'Cancer'  # 'Cancer' or 'Model'
highlight_cancers <- c('THCA', 'TGCT', 'KIRP', 'KIRC', 'ACC')   # NULL to not highlight any cancers

# prepare data
scores <- allscores[[classVar]]
dat <- as.data.frame(scores$roc_auc)

# merge (average) cancer types
cancers <- unlist(lapply(colnames(dat), function(x) head(unlist(strsplit(x, '_')), 1)))
uniq_cancers <- unique(cancers[duplicated(cancers)])  # only keep those with >1 stage comparisons
for (cancer in uniq_cancers) {
  dat[, cancer] <- apply(dat[, cancers %in% cancer, drop=F], 1, mean)
}
dat <- dat[, uniq_cancers]

# order data by decreasing cancer ROC AUC
o_cancer <- colnames(dat)[order(apply(dat, 2, median), decreasing=T)]
dat <- dat[, o_cancer] %>% rownames_to_column('Model')
dat <- pivot_longer(dat, cols=colnames(dat)[-1], values_to='ROC AUC', names_to='Cancer')
dat$Cancer <- factor(dat$Cancer, levels=o_cancer, ordered=T)
dat$Highlight <- factor(dat$Cancer %in% highlight_cancers)

# generate plot
pdf(file=paste0(fig_dir, '/', classVar, '_ROCAUC_CancerGrouped_boxplot.pdf'), width=3.75, height=4)
ggplot(dat, aes(x=Cancer, y=`ROC AUC`, fill=Highlight)) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position='none',
        axis.text=element_text(color='black', size=12),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        axis.title=element_text(size=12)) + 
  xlab('Cancer') +
  scale_fill_manual(values=brewer.pal(3, 'Paired'))
# coord_cartesian(ylim=c(0,1))
invisible(dev.off())


####################################################
### Fig. X, Panel X: Histogram of ROC AUC values ###
####################################################

# specify parameters
classVar <- 'mutTP53'  # 'mutTP53', 'cancerStatus', or 'tumorStage'

# prepare data
scores <- allscores[[classVar]]
dat <- as.data.frame(scores$roc_auc)
dat <- pivot_longer(dat, cols=colnames(dat), values_to='ROC AUC', names_to='Cancer')

# generate plot
pdf(file=paste0(fig_dir, '/', classVar, '_ROCAUC_histogram.pdf'), width=1.2, height=3.5)
ggplot(dat, aes(y=`ROC AUC`)) +
  # geom_histogram(aes(y=..density..), color='black', fill='white', bins=50) +
  geom_density(alpha=0.4, fill='steelblue3') +
  theme_classic() +
  theme(axis.text=element_text(color='black', size=12),
        axis.title=element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank()) +
  xlab('Density')
  # coord_cartesian(xlim=c(0,1))
invisible(dev.off())


##############################################################################
### Fig. X, Panel X: Combined histogram and two boxplots of ROC AUC values ###
##############################################################################

# specify parameters
classVar <- 'tumorStage'  # 'mutTP53', 'cancerStatus', or 'tumorStage'

for (groupby in c('Cancer', 'Model')) {
  
  # prepare boxplot data
  scores <- allscores[[classVar]]
  dat <- as.data.frame(scores$roc_auc)
  
  if (classVar == 'tumorStage') {
    # merge (average) cancer types
    cancers <- unlist(lapply(colnames(dat), function(x) head(unlist(strsplit(x, '_')), 1)))
    uniq_cancers <- unique(cancers)
    for (cancer in uniq_cancers) {
      dat[, cancer] <- apply(dat[, cancers %in% cancer, drop=F], 1, mean)
    }
    dat <- dat[, uniq_cancers]
  }
  
  o_model <- rownames(dat)[order(apply(dat, 1, median), decreasing=T)]
  o_cancer <- colnames(dat)[order(apply(dat, 2, median), decreasing=T)]
  dat <- dat[o_model, o_cancer] %>% rownames_to_column('Model')
  dat <- pivot_longer(dat, cols=colnames(dat)[-1], values_to='ROC AUC', names_to='Cancer')
  dat$Model <- factor(dat$Model, levels=o_model, ordered=T)
  dat$Cancer <- factor(dat$Cancer, levels=o_cancer, ordered=T)
  
  # generate plot
  if (groupby == 'Cancer') {
    hues <- c(0,100)
    text_angle <- 90
    vjust_val <- 0.5
  } else if (groupby == 'Model') {
    hues <- c(230,330)
    text_angle <- 45
    vjust_val <- 1
  }
  p <- ggplot(dat, aes_string(x=groupby, y="`ROC AUC`", fill=groupby)) +
    geom_boxplot() +
    theme_minimal() +
    theme(legend.position='none',
          axis.text=element_text(color='black', size=12),
          axis.text.x=element_text(angle=text_angle, hjust=1, vjust=vjust_val),
          axis.title=element_text(size=12),
          axis.title.x=element_blank()) + 
    scale_fill_hue(h=hues, c=75) +
    coord_cartesian(ylim=c(0,1))
  
  assign(paste0('p_box_', groupby), p)
}

# prepare histogram data
scores <- allscores[[classVar]]
dat <- as.data.frame(scores$roc_auc)

if (classVar == 'tumorStage') {
  # merge (average) cancer types
  cancers <- unlist(lapply(colnames(dat), function(x) head(unlist(strsplit(x, '_')), 1)))
  uniq_cancers <- unique(cancers)
  for (cancer in uniq_cancers) {
    dat[, cancer] <- apply(dat[, cancers %in% cancer, drop=F], 1, mean)
  }
  dat <- dat[, uniq_cancers]
}

dat <- pivot_longer(dat, cols=colnames(dat), values_to='ROC AUC', names_to='Cancer')

# generate plot
p_hist <- ggplot(dat, aes(y=`ROC AUC`)) +
  geom_density(alpha=0.4, fill='steelblue3') +
  theme_classic() +
  theme(axis.text=element_text(color='black', size=12),
        axis.title=element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank()) +
  xlab('Density') +
  coord_cartesian(ylim=c(0,1))

# generate combined plot
pdf(file=paste0(fig_dir, '/', classVar, '_ROCAUC_CombinedPlots.pdf'), width=10, height=3.5)
plot_grid(p_hist, p_box_Cancer, p_box_Model, ncol=3, rel_widths=c(4,18,8), align='h', axis='tb', scale = 0.98)
# plot_grid(p_hist, p_box_Cancer, ncol=2, rel_widths=c(2,9), align='h', axis='tb', scale = 0.98)
invisible(dev.off())



#######################################################################
### Fig. X, Panel X: Combined boxplots of top ML and DE gene scores ###
#######################################################################

# specify parameters
classVar <- 'cancerStatus'  # 'mutTP53', 'cancerStatus', or 'tumorStage'
n_genes <- 10  # specify number of genes to include
sort_by <- 'mean'  # 'mean' or 'top each'

for (model in c('Average', 'DE_FDRscore')) {
  
  # prepare data
  scores <- allscores[[classVar]]
  dat <- as.data.frame(scores[[model]])
  if (sort_by == 'mean') {
    o <- order(apply(dat, 1, mean), decreasing=T)
    top_genes <- rownames(dat)[o[1:n_genes]]
  } else if (sort_by == 'top each') {
    o <- apply(dat, 2, function(x) order(x, decreasing=T))
    top_genes <- rownames(dat)[unique(o[1:n_genes, ])]
  }
  if (classVar == 'mutTP53') {
    dat$targets <- as.factor(rownames(dat) %in% c('BAX', 'HSPA4L', 'KIF23'))
  } else {
    dat$targets <- 'NA'
  }
  dat <- dat[top_genes, ] %>% rownames_to_column('Gene')
  dat$Gene <- factor(dat$Gene, levels=rev(top_genes), ordered=T)
  dat <- pivot_longer(dat, cols=colnames(dat)[2:(ncol(dat)-1)], values_to='Score', names_to='Cancer')

  # specify color palette
  if (model == 'Average') {
    pal <- brewer.pal(4, 'Paired')[1:2]
    xlabel <- 'Consensus ML gene score'
  } else {
    pal <- brewer.pal(4, 'Paired')[3:4]
    xlabel <- 'Normalized DE gene score'
  }
  
  # generate plot
  p <- ggplot(dat, aes(x=Score, y=Gene, fill=targets)) +
    geom_boxplot() +
    theme_minimal() +
    theme(legend.position='none',
          axis.text=element_text(color='black', size=12),
          axis.title=element_text(size=12)) + 
    xlab(xlabel) +
    scale_fill_manual(values=pal) +
    coord_cartesian(xlim=c(0,1))
  
  assign(paste0('p_', model), p)
}

# generate combined boxplot
pdf(file=paste0(fig_dir, '/', classVar, '_combined_boxplots.pdf'), width=8.5, height=4.5)
plot_grid(p_Average, p_DE_FDRscore, ncol=2, rel_widths=c(1,1), align='h', axis='tb', scale = 0.98)
invisible(dev.off())



########################################################################
### Fig. X, Panel X: Bubble plot of gene DE fold-change and p-values ###
########################################################################

# specify parameters
classVar <- 'cancerStatus'  # 'mutTP53', 'cancerStatus', or 'tumorStage'
gene_start <- 'KIF'  # e.g., 'KIF', 'VAMP', 'STX', etc.
dist_metric <- 'euclidean'

# prepare data
genes <- rownames(gene_data)[startsWith(rownames(gene_data), gene_start)]
scores <- allscores[[classVar]]
DE_logFC <- as.data.frame(scores$DE_log2FC[rownames(scores$DE_log2FC) %in% genes, ])
DE_logPval <- as.data.frame(-log10(scores$DE_FDR[rownames(scores$DE_FDR) %in% genes, ]))

# cluster rows and columns
clust_method <- 'complete'
if (dist_metric == 'correlation') {
  gene_order <- rownames(DE_logFC)[hclust(as.dist(cor(t(DE_logFC))), method=clust_method)$order]
  cancer_order <- colnames(DE_logFC)[hclust(as.dist(cor(DE_logFC)), method=clust_method)$order]
} else {
  gene_order <- rownames(DE_logFC)[hclust(dist(DE_logFC, method=dist_metric), method=clust_method)$order]
  cancer_order <- colnames(DE_logFC)[hclust(dist(t(DE_logFC), method=dist_metric), method=clust_method)$order]
}

# finish preparing data
DE_logFC <- DE_logFC %>% rownames_to_column('Gene')
DE_logFC <- pivot_longer(DE_logFC, cols=colnames(DE_logFC)[2:ncol(DE_logFC)], values_to='logFC', names_to='Cancer')
DE_logPval <- DE_logPval %>% rownames_to_column('Gene')
DE_logPval <- pivot_longer(DE_logPval, cols=colnames(DE_logPval)[2:ncol(DE_logPval)], values_to='logPval', names_to='Cancer')
dat <- merge(DE_logFC, DE_logPval, by=c('Gene', 'Cancer'), )

dat$Gene <- ordered(dat$Gene, levels=gene_order)
dat$Cancer <- ordered(dat$Cancer, levels=cancer_order)

dat$logPval[dat$logPval > 5] <- 5  # trim max pval to scale
dat$logFC[dat$logFC > 2] <- 2  # trim min logFC to scale
dat$logFC[dat$logFC < -2] <- -2  # trim max logFC to scale

if (nlevels(dat$Gene) > 12) {
  legend_orientation <- 'vertical'
  plot_width <- 4.2+0.1*nlevels(dat$Cancer)
} else {
  legend_orientation <- 'horizontal'
  plot_width <- 5+0.1*nlevels(dat$Cancer)
}

# generate bubble plot
pdf(file=paste0(fig_dir, '/', classVar, '_', gene_start, '_bubble.pdf'),
    width=plot_width,
    height=0.9+1.7/7*nlevels(dat$Gene))
ggplot(dat, aes(x=Cancer, y=Gene, size=logPval, fill=logFC)) +
  geom_point(shape=21, stroke=0.5) +
  scale_size(range=c(1, 7), limits=c(0,5), name='-log10(p)') +
  scale_fill_gradientn(colors=rev(brewer.pal(9,'RdBu')), limits=c(-2,2), name='log2FC') +
  theme_classic() +
  theme(axis.text=element_text(color='black', size=12),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        axis.title=element_text(size=12),
        legend.text=element_text(color='black', size=12),
        legend.box=legend_orientation)
invisible(dev.off())



##################################################################################
### Plot expression of gene(s), grouped by cancer and subgrouped by a classVar ###
##################################################################################

# specify parameters
gene <- 'KIF20A' #c('KIF20A', 'KIF23') # c('AGR2', 'NET1', 'GALNT6', 'GALNT14', 'GALNT18')
classVar <- 'CancerStatus'  # 'mutTP53'
classLevels <- c('Primary solid Tumor', 'Solid Tissue Normal')  # c('TRUE', 'FALSE')
cancers <- 'all'  # use 'all' to include all possible cancers
useFacets <- T

# filter data
dat <- expdata %>% select(all_of(gene), all_of(classVar), 'Project')
if (classVar != 'CancerStatus') {
  dat <- dat[expdata$CancerStatus == 'Primary solid Tumor', ]
}
dat <- dat[dat[, classVar] %in% classLevels, ]

# log-transform counts
dat[, gene] <- log2(dat[, gene] + 1)

# remove cancer types without at least 10 samples in each class level
if (length(cancers) == 1 && casefold(cancers) == 'all') {
  cancers <- unique(na.omit(dat$Project))
}
tab <- table(dat[, c(classVar, 'Project')])[classLevels, ]
cancers <- intersect(cancers, colnames(tab)[colSums(tab >= 10) == 2])
if (length(cancers) == 0) {
  stop('No cancers remain after removing those with fewer than 10 samples in each class!')
}
dat <- dat[dat$Project %in% cancers, ]

# rename CancerStatus class levels if used
if (classVar == 'CancerStatus') {
  indx <- dat$CancerStatus == 'Primary solid Tumor'
  dat$CancerStatus[indx] = 'Tumor'
  dat$CancerStatus[!indx] = 'Normal'
}

# convert matrix to long form where gene names are in a new column
dat <- pivot_longer(dat, cols=all_of(gene), names_to='Gene', values_to='log2(TPM)')


# generate plot
if (length(cancers) == 1) {
  
  if (useFacets) {
    
    ggplot(dat, aes_string(x="Gene", y="`log2(TPM)`", fill=classVar)) +
      geom_boxplot() +
      facet_grid(. ~ Gene, scales='free') + 
      scale_x_discrete(position = "top") +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
    
  } else {
    
    ggplot(dat, aes_string(x="Gene", y="`log2(TPM)`", fill=classVar)) +
      geom_boxplot()
    
  }
  
} else if (length(gene) > 1) {
  
  ggplot(dat, aes_string(x="Project", y="`log2(TPM)`", fill=classVar)) +
    geom_boxplot() +
    facet_grid(Gene ~ ., scales='free')
  
  ggplot(dat, aes_string(x="Project", y="`log2(TPM)`", fill=classVar)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(Gene ~ ., scales='free')# +
    #coord_cartesian(ylim=c(0,3.1))
  
} else {
  
  ggplot(dat, aes_string(x="Project", y="`log2(TPM)`", fill=classVar)) +
    geom_boxplot()
}


###############################################################################
### Plot expression of HAS genes, grouped by cancer and subgrouped by stage ###
###############################################################################

# specify parameters
gene <- c('HAS1','HAS2','HAS3')
classVar <- 'TumorStageMerged'  # 'CancerStatus', 'mutTP53', 'TumorStageMerged'
classLevels <- c('stage i', 'stage ii', 'stage iii', 'stage iv')   # c('Primary solid Tumor', 'Solid Tissue Normal')  # c('TRUE', 'FALSE')
cancers <- c('TGCT', 'THCA')  # use 'all' to include all possible cancers
pal <- brewer.pal(9, 'YlOrRd')[c(3,5,7,9)]

# filter data
dat <- expdata %>% select(all_of(gene), all_of(classVar), 'Project')
dat <- dat[expdata$CancerStatus == 'Primary solid Tumor', ]
dat <- dat[dat[, classVar] %in% classLevels, ]
dat <- dat[dat$Project %in% cancers, ]

# rename levels
new_levels <- c('I', 'II', 'III', 'IV')
dat$TumorStageMerged <- new_levels[match(dat$TumorStageMerged, classLevels)]

# log-transform counts
dat[, gene] <- log2(dat[, gene] + 1)

# convert matrix to long form where gene names are in a new column
dat <- pivot_longer(dat, cols=all_of(gene), names_to='Gene', values_to='log2(TPM)')

# generate plots
p1 <- dat %>% filter(Gene == 'HAS1') %>% 
  ggplot(aes_string(x="Project", y="`log2(TPM)`", fill=classVar)) +
  geom_boxplot(outlier.shape = NA, color='black') +
  theme_minimal() +
  theme(legend.position='top',
        axis.text=element_text(size=12, color='black'),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank()) +
  scale_fill_manual(values=pal) +
  labs(fill='Stage') +
  coord_cartesian(ylim=c(0,3.1)) +
  ylab('HAS1 log2(TPM)')

p2 <- dat %>% filter(Gene == 'HAS2') %>% 
  ggplot(aes_string(x="Project", y="`log2(TPM)`", fill=classVar)) +
  geom_boxplot(outlier.shape = NA, color='black') +
  theme_minimal() +
  theme(legend.position='none',
        axis.text=element_text(size=12, color='black'),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank()) +
  scale_fill_manual(values=pal) +
  coord_cartesian(ylim=c(0,7.1)) +
  ylab('HAS2 log2(TPM)')

p3 <- dat %>% filter(Gene == 'HAS3') %>% 
  ggplot(aes_string(x="Project", y="`log2(TPM)`", fill=classVar)) +
  geom_boxplot(outlier.shape = NA, color='black') +
  theme_minimal() +
  theme(legend.position='none',
        axis.text=element_text(size=12, color='black')) +
  scale_fill_manual(values=pal) +
  coord_cartesian(ylim=c(0,6.7)) +
  xlab('Cancer type') +
  ylab('HAS3 log2(TPM)')

pdf(file=paste0(fig_dir, '/stage_expression_HAS.pdf'), width=3.5, height=5.5)
plot_grid(p1, p2, p3, align='h', axis='lr', ncol=1, rel_heights=c(8.5,7,8))
invisible(dev.off())



########################################################################################################
### Plot expression of CRYAB, grouped by cancer, subgrouped by classVar, and faceted by DE direction ###
########################################################################################################


# specify parameters
gene <- 'CRYAB'  # THIS SECTION IS HARD-CODED FOR CRYAB
classVar <- 'CancerStatus'  # 'mutTP53'
classLevels <- c('Primary solid Tumor', 'Solid Tissue Normal')  # c('TRUE', 'FALSE')
cancers <- c('BRCA','LUAD','COAD','BLCA','HNSC','STAD',
             'THCA','PRAD','READ','UCEC','ESCA','KICH',
             'LUSC','LIHC','KIRC','KIRP')  # ordered by DE direction and p-value
cancers_down <- cancers[1:10]
cancers_ns <- cancers[11:13]
cancers_up <- cancers[14:16]


# filter data
dat <- expdata %>% select(all_of(gene), all_of(classVar), 'Project')
if (classVar != 'CancerStatus') {
  dat <- dat[expdata$CancerStatus == 'Primary solid Tumor', ]
}
dat <- dat[dat[, classVar] %in% classLevels, ]

# log-transform counts
dat[, gene] <- log2(dat[, gene] + 1)

# order cancer types by DE significance and direction
dat$Project <- factor(dat$Project, levels=cancers, ordered=T)

# remove cancer types without at least 10 samples in each class level
if (length(cancers) == 1 && casefold(cancers) == 'all') {
  cancers <- unique(na.omit(dat$Project))
}
tab <- table(dat[, c(classVar, 'Project')])[classLevels, ]
cancers <- intersect(cancers, colnames(tab)[colSums(tab >= 10) == 2])
if (length(cancers) == 0) {
  stop('No cancers remain after removing those with fewer than 10 samples in each class!')
}
dat <- dat[dat$Project %in% cancers, ]

# rename CancerStatus class levels if used
if (classVar == 'CancerStatus') {
  indx <- dat$CancerStatus == 'Primary solid Tumor'
  dat$CancerStatus[indx] = 'Tumor'
  dat$CancerStatus[!indx] = 'Normal'
}

# add column specifying DE direction
DE_labels <- c('Decreased', 'Not Significant', 'Increased')
dat$DE <- NA
dat$DE[dat$Project %in% cancers_down] <- DE_labels[1]
dat$DE[dat$Project %in% cancers_ns] <- DE_labels[2]
dat$DE[dat$Project %in% cancers_up] <- DE_labels[3]
dat$DE <- factor(dat$DE, levels=DE_labels, ordered=T)

# convert matrix to long form where gene names are in a new column
dat <- pivot_longer(dat, cols=all_of(gene), names_to='Gene', values_to='CRYAB log2(TPM)')

# generate plot
pal <- brewer.pal(9, 'YlOrRd')[c(1,8)]
pdf(paste0(fig_dir, '/CRYAB_expression_boxplot.pdf'), height=3, width=10)
ggplot(dat, aes_string(x="Project", y="`CRYAB log2(TPM)`", fill=classVar)) +
  geom_boxplot(color='black') + 
  facet_grid(. ~ DE, scales='free_x', space='free_x') +
  scale_fill_manual(values=pal) +
  theme_bw() +
  theme(axis.text=element_text(size=10, color='black'),
        strip.text=element_text(size=11, color='black')) +
  xlab('Cancer type')
invisible(dev.off())


#####################################################################
### Plot expression of gene(s) grouped by normal and tumor stages ###
#####################################################################

# specify parameters
gene <- 'KIF20A' #c('KIF20A', 'KIF23') # c('AGR2', 'NET1', 'GALNT6', 'GALNT14', 'GALNT18')
classLevels <- c('normal', 'stage i', 'stage ii', 'stage iii', 'stage iv')
classVar <- 'TumorStageMerged'
cancers <- 'all'  # use 'all' to include all possible cancers
pal <- brewer.pal(9, 'YlOrRd')[c(1,3,5,7,9)]

# filter data
dat <- expdata
dat$TumorStageMerged[dat$CancerStatus == 'Solid Tissue Normal'] <- 'normal'
dat <- dat %>% select(all_of(gene), all_of(classVar), 'Project')
dat <- dat[dat[, classVar] %in% classLevels, ]

# log-transform counts
dat[, gene] <- log2(dat[, gene] + 1)

# remove cancer types without tumor stage data
if (length(cancers) == 1 && casefold(cancers) == 'all') {
  cancers <- unique(na.omit(dat$Project))
}
keepCancers <- NULL
for (cancer in cancers) {
  if ((sum(grepl('stage', dat$TumorStageMerged[dat$Project == cancer])) > 10) &&
      (sum(grepl('normal', dat$TumorStageMerged[dat$Project == cancer])) > 1)) {
    keepCancers <- c(keepCancers, cancer)
  }
}
cancers <- keepCancers
dat <- dat[dat$Project %in% cancers, ]

# rename levels
new_levels <- c('Normal', 'Stage I', 'Stage II', 'Stage III', 'Stage IV')
dat$TumorStageMerged <- new_levels[match(dat$TumorStageMerged, classLevels)]

# add column specifying cancer group (renal or non-renal)
dat$renal <- 'Other'
dat$renal[dat$Project %in% c('KIRC','KIRP','KICH')] <- 'Renal Carcinomas'
dat$renal <- factor(dat$renal, levels=c('Renal Carcinomas', 'Other'), ordered=T)


# convert matrix to long form where gene names are in a new column
dat <- pivot_longer(dat, cols=all_of(gene), names_to='Gene', values_to='log2(TPM)')


# generate plot
pdf(file=paste0(fig_dir, '/stage_expression_KIF20A.pdf'), width=12, height=4)
ggplot(dat, aes_string(x="Project", y="`log2(TPM)`", fill=classVar)) +
  geom_boxplot(outlier.shape=NA, color='black') +
  theme_bw() +
  theme(axis.text=element_text(size=12, color='black'),
        axis.text.x=element_text(size=10),
        strip.text=element_text(size=11, color='black'),
        legend.position='bottom',
        legend.justification='right') +
  scale_fill_manual(values=pal) +
  facet_grid(. ~ renal, scales='free_x', space='free_x') +
  labs(fill='Tumor stage:') +
  xlab('Cancer type') +
  ylab(paste(gene, 'log2(TPM)'))
invisible(dev.off())
  
  
  

