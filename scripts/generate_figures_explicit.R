
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


# specify main project directory
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

score_vars <- c('scores_cancerStatus', 'scores_mutTP53', 'scores_tumorStage')
# DE_vars <- c('DE_cancerStatus', 'DE_mutTP53', 'DE_tumorStage')
model_names <- c('ExtraTreesClassifier', 'RandomForestClassifier', 'AdaBoostClassifier', 'XGBClassifier',
                 'LinearDiscriminantAnalysis', 'SVC', 'LassoRegression', 'RidgeRegression', 'Average')
new_model_names <- c('Extra Trees', 'Random Forest', 'AdaBoost', 'XGBoost', 'LDA', 'SVM', 'Lasso', 'Ridge', 'Average')
DE_vars <- c('DE_log2FC', 'DE_FDR', 'DE_FDRscore')

for (x in score_vars) assign(x, list())
for (model in c(model_names, DE_vars)) {
  scores_cancerStatus[[model]] <- matrix(NA, nrow=length(genes), ncol=length(cancers), dimnames=list(genes, cancers))
  scores_mutTP53[[model]]      <- matrix(NA, nrow=length(genes), ncol=length(cancers), dimnames=list(genes, cancers))
  scores_tumorStage[[model]]   <- matrix(NA, nrow=length(genes), ncol=length(cancer_stages), dimnames=list(genes, cancer_stages))
}
scores_cancerStatus$roc_auc <- matrix(NA, nrow=length(model_names)-1, ncol=length(cancers), dimnames=list(model_names[1:(length(model_names)-1)], cancers))
scores_mutTP53$roc_auc      <- matrix(NA, nrow=length(model_names)-1, ncol=length(cancers), dimnames=list(model_names[1:(length(model_names)-1)], cancers))
scores_tumorStage$roc_auc   <- matrix(NA, nrow=length(model_names)-1, ncol=length(cancer_stages), dimnames=list(model_names[1:(length(model_names)-1)], cancer_stages))

for (cancer in cancers) {
  cancer_path <- paste0(proj_dir, '/results/', cancer)
  files_cancerStatus <- dir(cancer_path, 'CancerStatus')
  files_mutTP53 <- dir(cancer_path, 'mutTP53')
  files_tumorStage <- dir(cancer_path, 'TumorStage')
  
  if (length(files_cancerStatus) > 0) {
    scores <- read.csv(paste0(cancer_path, '/', files_cancerStatus[grepl('GenesRanking', files_cancerStatus)]), row.names=1)
    for (model in model_names) {
      scores_cancerStatus[[model]][, cancer] <- scores[match(genes, rownames(scores)), model]
    }
    de_res <- read.delim(paste0(cancer_path, '/', files_cancerStatus[grepl('DEresults', files_cancerStatus)]), row.names=1)
    scores_cancerStatus$DE_log2FC[, cancer] <- de_res$logFC[match(genes, rownames(de_res))]
    scores_cancerStatus$DE_FDR[, cancer] <- de_res$FDR[match(genes, rownames(de_res))]
    scores_cancerStatus$DE_FDRscore[, cancer] <- score_pvals(de_res$FDR[match(genes, rownames(de_res))])
    
    roc_auc <- read.csv(paste0(cancer_path, '/', files_cancerStatus[grepl('ROCAUC', files_cancerStatus)]), row.names=1)
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
    
    roc_auc <- read.csv(paste0(cancer_path, '/', files_mutTP53[grepl('ROCAUC', files_mutTP53)]), row.names=1)
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
      
      roc_auc <- read.csv(paste0(cancer_path, '/', sub('GenesRanking', 'CVscores_ROCAUC', f)), row.names=1)
      scores_tumorStage$roc_auc[, f_name] <- roc_auc$Score[match(model_names[1:(length(model_names)-1)], rownames(roc_auc))]
    }
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

# rename models in roc_auc slot
rownames(scores_cancerStatus$roc_auc) <- new_model_names[match(rownames(scores_cancerStatus$roc_auc), model_names)]
rownames(scores_mutTP53$roc_auc) <- new_model_names[match(rownames(scores_mutTP53$roc_auc), model_names)]
rownames(scores_tumorStage$roc_auc) <- new_model_names[match(rownames(scores_tumorStage$roc_auc), model_names)]

# combine scores into list and remove intermediate variables
allscores <- list(cancerStatus=scores_cancerStatus, mutTP53=scores_mutTP53, tumorStage=scores_tumorStage)
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


#########################################################
### Fig. 2, Panel A: Histogram of mutTP53 gene scores ###
#########################################################

# specify parameters
classVar <- 'mutTP53'  # 'mutTP53', 'cancerStatus', or 'tumorStage'
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



#############################################################
### Fig. 2, Panel B: Heatmap of top-scoring mutTP53 genes ###
#############################################################

# specify parameters
classVar <- 'mutTP53'  # 'mutTP53', 'cancerStatus', or 'tumorStage'
n_genes <- 10  # specify number of genes to include
model <- 'Average'  # e.g., 'Average' or 'DE_FDRscore'

# prepare data
scores <- allscores[[classVar]]
dat <- as.data.frame(scores[[model]])
o <- order(apply(dat, 1, mean), decreasing=T)
top_genes <- rownames(dat)[o[1:n_genes]]
dat <- dat[top_genes, ]

# generate heatmap
pdf(file=paste0(fig_dir, '/', classVar, '_', model, '_heatmap.pdf'), width=3+2.3*ncol(dat)/22, height=3, onefile=F)
pheatmap(dat,
         scale='none',
         color=magma(100),
         clustering_distance_rows='euclidean',
         clustering_distance_cols='euclidean',
         clustering_method='ward.D2',
         breaks=seq(0,0.8,len=100),
         angle_col=90,
         border_color='black')
invisible(dev.off())



#############################################################################################
### Fig. 2, Panels C-E: Combined histogram and two boxplots of ROC AUC values for mutTP53 ###
#############################################################################################

# specify parameters
classVar <- 'mutTP53'  # 'mutTP53', 'cancerStatus', or 'tumorStage'

for (groupby in c('Cancer', 'Model')) {
  
  # prepare boxplot data
  scores <- allscores[[classVar]]
  dat <- as.data.frame(scores$roc_auc)
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
dat <- pivot_longer(dat, cols=colnames(dat), values_to='ROC AUC', names_to='Cancer')

# generate plot
p_hist <- ggplot(dat, aes(y=`ROC AUC`)) +
  geom_density(fill=brewer.pal(3, 'Paired')[1]) +
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
invisible(dev.off())



##############################################################################################
### Fig. 4: Heatmap of top-5 genes for each cancer type, CancerStatus, subset of 5 cancers ###
##############################################################################################

# specify parameters
classVar <- 'cancerStatus'
n_genes <- 5
model <- 'Average'
cancers <- c('STAD', 'READ', 'COAD', 'KICH', 'THCA')

# prepare data
scores <- allscores[[classVar]]
dat <- as.data.frame(scores[[model]])
dat <- dat[, colnames(dat) %in% cancers]
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
annColors <- list(Function = viridis(3, begin=0.4) %>% setNames(unique(gene_data$module[row_order])))
annData <- gene_data[, c('module'), drop=F] %>% setNames('Function')

# generate heatmap
pdf(file=paste0(fig_dir, '/', classVar, '_', model, '_topEach_heatmap.pdf'), width=3+2.3*ncol(dat)/22, height=4, onefile=F)
pheatmap(dat,
         scale='none',
         color=magma(100),
         cluster_rows=F,
         cluster_cols=F,
         breaks=seq(0,0.8,len=100),
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








###########################################################################
### Fig. SX: Heatmap of gene scores for top 5 genes of each cancer type ###
###########################################################################

# specify parameters
classVar <- 'cancerStatus'  # 'mutTP53', 'cancerStatus', or 'tumorStage'
n_genes <- 5  # specify number of genes to include
model <- 'Average'  # e.g., 'Average' or 'DE_FDRscore'

# prepare data
scores <- allscores[[classVar]]
dat <- as.data.frame(scores[[model]])
o <- apply(dat, 2, function(x) order(x, decreasing=T))
top_genes <- rownames(dat)[unique(as.vector(o[1:n_genes, ]))]
dat <- dat[top_genes, ]

# cluster columns based on correlation distance
col_order <- hclust(as.dist(cor(dat, method='pearson')), method='average')$order
dat <- dat[, col_order]
row_order <- NULL
for (i in seq(ncol(dat))) {
  ro <- order(dat[, i], decreasing=T)
  row_order <- c(row_order, setdiff(ro[1:n_genes], row_order))
}
dat <- dat[row_order, ]

# generate heatmap
pdf(file=paste0(fig_dir, '/', classVar, '_', model, '_topEach_allcancers_heatmap.pdf'), width=4, height=10)
pheatmap(dat,
         scale='none',
         color=magma(100),
         cluster_rows=F,
         cluster_cols=F,
         breaks=seq(0,0.8,len=100),
         angle_col=90,
         border_color='black')
invisible(dev.off())



####################################################################
### Fig. SX: Bubble plot of KIF gene DE fold-change and p-values ###
####################################################################

# specify parameters
classVar <- 'cancerStatus'  # 'mutTP53', 'cancerStatus', or 'tumorStage'
gene_start <- 'KIF'
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

# generate bubble plot
pdf(file=paste0(fig_dir, '/', classVar, '_', gene_start, '_bubble.pdf'), width=4.9, height=2.8)
ggplot(dat, aes(x=Cancer, y=Gene, size=logPval, fill=logFC)) +
  geom_point(shape=21, stroke=0.5) +
  scale_size(range=c(1, 7), limits=c(0,5), name='-log10(p)', guide=guide_legend(nrow=1)) +
  scale_fill_gradientn(colors=rev(brewer.pal(9,'RdBu')), limits=c(-2,2), name='log2FC') +
  theme_classic() +
  theme(axis.text=element_text(color='black', size=12),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        axis.title=element_text(size=12),
        legend.text=element_text(color='black', size=12),
        legend.box='vertical',
        legend.position='bottom',
        legend.box.just='left')
invisible(dev.off())



