
# load packages (installed in "psp-cancer-r" conda environment)
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
fig_dir <- paste0(proj_dir, '/doc/manuscript/REVISION/figures/fig_pieces')


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

scores_cancerStatus$model_score <- matrix(NA, nrow=length(model_names)-1, ncol=length(cancers), dimnames=list(model_names[1:(length(model_names)-1)], cancers))
scores_mutTP53$model_score      <- matrix(NA, nrow=length(model_names)-1, ncol=length(cancers), dimnames=list(model_names[1:(length(model_names)-1)], cancers))
scores_tumorStage$model_score   <- matrix(NA, nrow=length(model_names)-1, ncol=length(cancer_stages), dimnames=list(model_names[1:(length(model_names)-1)], cancer_stages))
scores_stageRegress$model_score <- matrix(NA, nrow=length(reg_model_names)-1, ncol=length(cancers), dimnames=list(reg_model_names[1:(length(reg_model_names)-1)], cancers))

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
    
    model_score <- read.csv(paste0(cancer_path, '/', files_cancerStatus[grepl('CVscores', files_cancerStatus)]), row.names=1)
    scores_cancerStatus$model_score[, cancer] <- model_score$Score[match(model_names[1:(length(model_names)-1)], rownames(model_score))]
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
      de_res <- read.delim(paste0(cancer_path, '/', sub('GenesRanking.csv', 'DEresults.txt', f)), row.names=1)
      scores_tumorStage$DE_log2FC[, f_name] <- de_res$logFC[match(genes, rownames(de_res))]
      scores_tumorStage$DE_FDR[, f_name] <- de_res$FDR[match(genes, rownames(de_res))]
      scores_tumorStage$DE_FDRscore[, f_name] <- score_pvals(de_res$FDR[match(genes, rownames(de_res))])
      
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

# rename models in model_score slot
rownames(scores_cancerStatus$model_score) <- new_model_names[match(rownames(scores_cancerStatus$model_score), model_names)]
rownames(scores_mutTP53$model_score) <- new_model_names[match(rownames(scores_mutTP53$model_score), model_names)]
rownames(scores_tumorStage$model_score) <- new_model_names[match(rownames(scores_tumorStage$model_score), model_names)]
rownames(scores_stageRegress$model_score) <- new_reg_model_names[match(rownames(scores_stageRegress$model_score), reg_model_names)]

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



#################################################################################################
### Fig. 2, Panels C-E: Combined histogram and two boxplots of model score values for mutTP53 ###
#################################################################################################

# specify parameters
classVar <- 'mutTP53'  # 'mutTP53', 'cancerStatus', or 'tumorStage'
score_name <- 'ROC AUC'  # name of model scoring metric for use in labeling plot

# remove problematic chacters from score_name
score_var <- gsub(' |[.]', '', score_name)

for (groupby in c('Cancer', 'Model')) {
  
  # prepare boxplot data
  scores <- allscores[[classVar]]
  dat <- as.data.frame(scores$model_score)
  o_model <- rownames(dat)[order(apply(dat, 1, median), decreasing=T)]
  o_cancer <- colnames(dat)[order(apply(dat, 2, median), decreasing=T)]
  dat <- dat[o_model, o_cancer] %>% rownames_to_column('Model')
  dat <- pivot_longer(dat, cols=colnames(dat)[-1], values_to=score_var, names_to='Cancer')
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
  p <- ggplot(dat, aes_string(x=groupby, y=score_var, fill=groupby)) +
    geom_boxplot() +
    theme_minimal() +
    theme(legend.position='none',
          axis.text=element_text(color='black', size=12),
          axis.text.x=element_text(angle=text_angle, hjust=1, vjust=vjust_val),
          axis.title=element_text(size=12),
          axis.title.x=element_blank()) + 
    scale_fill_hue(h=hues, c=75) +
    ylab(score_name) +
    coord_cartesian(ylim=c(0,1))
  
  assign(paste0('p_box_', groupby), p)
}

# prepare histogram data
scores <- allscores[[classVar]]
dat <- as.data.frame(scores$model_score)
dat <- pivot_longer(dat, cols=colnames(dat), values_to=score_var, names_to='Cancer')

# generate plot
p_hist <- ggplot(dat, aes_string(y=score_var)) +
  geom_density(fill=brewer.pal(3, 'Paired')[1]) +
  theme_classic() +
  theme(axis.text=element_text(color='black', size=12),
        axis.title=element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank()) +
  xlab('Density') +
  ylab(score_name) +
  coord_cartesian(ylim=c(0,1))

# generate combined plot
pdf(file=paste0(fig_dir, '/', classVar, '_', score_var, '_CombinedPlots.pdf'), width=10, height=3.5)
plot_grid(p_hist, p_box_Cancer, p_box_Model, ncol=3, rel_widths=c(4,18,8), align='h', axis='tb', scale = 0.98)
invisible(dev.off())


##################################################################################
### Fig. 3, Panels A-B: Boxplots of top ML and DE gene scores for CancerStatus ###
##################################################################################

# specify parameters
classVar <- 'cancerStatus'
n_genes <- 10  # number of genes to include
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


##############################################################################
### Fig. 3, Panel C: Bubble plot of STX genes DE fold-changes and p-values ###
##############################################################################

# specify parameters
classVar <- 'cancerStatus'
gene_start <- 'STX'
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
dat <- merge(DE_logFC, DE_logPval, by=c('Gene', 'Cancer'))

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


#########################################################################################
### Fig. 5, Panel A: Boxplot of tumorStage model score values, grouped by cancer type ###
#########################################################################################

# specify parameters
classVar <- 'tumorStage'
groupby <- 'Cancer'  # 'Cancer' or 'Model'
highlight_cancers <- c('THCA', 'TGCT', 'KIRP', 'KIRC', 'ACC')   # NULL to not highlight any cancers
score_name <- 'ROC AUC'  # name of model scoring metric for use in labeling plot

# remove problematic chacters from score_name
score_var <- gsub(' |[.]', '', score_name)

# prepare data
scores <- allscores[[classVar]]
dat <- as.data.frame(scores$model_score)

# merge (average) cancer types
cancers <- unlist(lapply(colnames(dat), function(x) head(unlist(strsplit(x, '_')), 1)))
uniq_cancers <- unique(cancers[duplicated(cancers)])  # only keep those with >1 stage comparisons
# uniq_cancers <- unique(cancers)  # keep all cancers
for (cancer in uniq_cancers) {
  dat[, cancer] <- apply(dat[, cancers %in% cancer, drop=F], 1, mean)
}
dat <- dat[, uniq_cancers]

# order data by decreasing cancer model score
o_cancer <- colnames(dat)[order(apply(dat, 2, median), decreasing=T)]
dat <- dat[, o_cancer] %>% rownames_to_column('Model')
dat <- pivot_longer(dat, cols=colnames(dat)[-1], values_to=score_var, names_to='Cancer')
dat$Cancer <- factor(dat$Cancer, levels=o_cancer, ordered=T)
dat$Highlight <- factor(dat$Cancer %in% highlight_cancers)

# generate plot
pdf(file=paste0(fig_dir, '/', classVar, '_', score_var, '_CancerGrouped_boxplot.pdf'), width=3.75, height=4)
ggplot(dat, aes_string(x='Cancer', y=score_var, fill='Highlight')) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position='none',
        axis.text=element_text(color='black', size=12),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        axis.title=element_text(size=12)) + 
  xlab('Cancer') +
  ylab(score_name) +
  scale_fill_manual(values=brewer.pal(3, 'Paired'))
# coord_cartesian(ylim=c(0,1))
invisible(dev.off())


####################################################################
### Fig. 5, Panel B: Heatmap of top-scoring genes for tumorStage ###
####################################################################

# specify parameters
classVar <- 'tumorStage'
group_cancers <- F  # group (average) tumor stages together by cancer type
n_genes <- 10  # number of genes to include
sort_by <- 'mean'  # 'mean' or 'top each'
model <- 'Average'  # e.g., 'Average' or 'DE_FDRscore'
cancer_types <- c('ACC','KIRP','KIRC','THCA','TGCT')

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



##############################################################################################
### Fig. 5, Panel C: Plot KIF20A TPM for all cancer types grouped by normal & tumor stages ###
##############################################################################################

# specify parameters
gene <- 'KIF20A' 
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



######################################################################
#_____________________________________________________________________
#
#                        SUPPLEMENTARY FIGURES
#_____________________________________________________________________


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
dat <- merge(DE_logFC, DE_logPval, by=c('Gene', 'Cancer'))

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


###########################################################################################
### Fig. SX, Panels XXX: Combined histogram and boxplot of gene scores for stageRegress ###
###########################################################################################

# specify parameters
classVar <- 'stageRegress'
model <- 'Average'
n_genes <- 10

# prepare data
scores <- allscores[[classVar]]
dat <- scores[[model]]
dat <- as.data.frame(apply(dat, 1, mean)) %>% setNames('Mean score')

hist_xlab <- 'Mean gene ML score'
box_xlab <- 'Gene ML score'

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
o <- order(apply(dat, 1, mean), decreasing=T)
top_genes <- rownames(dat)[o[1:n_genes]]

dat <- dat[top_genes, ] %>% rownames_to_column('Gene')
dat$Gene <- factor(dat$Gene, levels=rev(top_genes), ordered=T)
dat <- pivot_longer(dat, cols=colnames(dat)[2:(ncol(dat)-1)], values_to='Score', names_to='Cancer')

# generate boxplots
p_box <- ggplot(dat, aes(x=Score, y=Gene)) +
  geom_boxplot(color='black', fill='#A6CEE3') +
  theme_minimal() +
  theme(legend.position='none',
        axis.text=element_text(color='black', size=12),
        axis.title=element_text(size=12)) + 
  xlab(box_xlab) +
  coord_cartesian(xlim=c(0,1))

# generate combined plot
pdf(file=paste0(fig_dir, '/', classVar, '_', model, '_combined_HistBox.pdf'), width=7, height=3.5)
plot_grid(p_hist, p_box, ncol=2, rel_widths=c(1,1), align='h', axis='tb', scale = 0.98)
invisible(dev.off())



#######################################################################
### Fig. SX, Panel X: Heatmap of top-scoring genes for stageRegress ###
#######################################################################

# specify parameters
classVar <- 'stageRegress'
n_genes <- 10
sort_by <- 'mean'
model <- 'Average'

# prepare data
scores <- allscores[[classVar]]
dat <- as.data.frame(scores[[model]])
o <- order(apply(dat, 1, mean), decreasing=T)
top_genes <- rownames(dat)[o[1:n_genes]]
dat <- dat[top_genes, ]

# generate heatmap
pdf(file=paste0(fig_dir, '/', classVar, '_', model, '_heatmap.pdf'), width=2.8+2*ncol(dat)/22, height=3, onefile=F)
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



#######################################################################################################
### Fig. SX, Panels XXX: Combined histogram and two boxplots of model score values for stageRegress ###
#######################################################################################################

# specify parameters
classVar <- 'stageRegress'
score_name <- 'Neg. MSE'  # name of model scoring metric for use in labeling plot

# remove problematic chacters from score_name
score_var <- gsub(' |[.]', '', score_name)

for (groupby in c('Cancer', 'Model')) {
  
  # prepare boxplot data
  scores <- allscores[[classVar]]
  dat <- as.data.frame(scores$model_score)
  o_model <- rownames(dat)[order(apply(dat, 1, median), decreasing=T)]
  o_cancer <- colnames(dat)[order(apply(dat, 2, median), decreasing=T)]
  dat <- dat[o_model, o_cancer] %>% rownames_to_column('Model')
  dat <- pivot_longer(dat, cols=colnames(dat)[-1], values_to=score_var, names_to='Cancer')
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
  p <- ggplot(dat, aes_string(x=groupby, y=score_var, fill=groupby)) +
    geom_boxplot() +
    theme_minimal() +
    theme(legend.position='none',
          axis.text=element_text(color='black', size=12),
          axis.text.x=element_text(angle=text_angle, hjust=1, vjust=vjust_val),
          axis.title=element_text(size=12),
          axis.title.x=element_blank()) + 
    ylab(score_name) +
    scale_fill_hue(h=hues, c=75)
  
  assign(paste0('p_box_', groupby), p)
}

# prepare histogram data
scores <- allscores[[classVar]]
dat <- as.data.frame(scores$model_score)
dat <- pivot_longer(dat, cols=colnames(dat), values_to=score_var, names_to='Cancer')

# generate plot
p_hist <- ggplot(dat, aes_string(y=score_var)) +
  geom_density(fill=brewer.pal(3, 'Paired')[1]) +
  theme_classic() +
  theme(axis.text=element_text(color='black', size=12),
        axis.title=element_text(size=12),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank()) +
  xlab('Density') +
  ylab(score_name)

# generate combined plot
pdf(file=paste0(fig_dir, '/', classVar, '_', score_var, '_CombinedPlots.pdf'), width=7, height=3.5)
plot_grid(p_hist, p_box_Cancer, p_box_Model, ncol=3, rel_widths=c(4,10,7), align='h', axis='tb', scale = 0.98)
invisible(dev.off())


###############################################################################
### Fig. SX, Panel X: Heatmap of top cancer-specific genes for stageRegress ###
###############################################################################

# specify parameters
classVar <- 'stageRegress'
n_genes <- 3
model <- 'Average'

# prepare data
scores <- allscores[[classVar]]
dat <- as.data.frame(scores[[model]])
o <- apply(dat, 2, function(x) order(x, decreasing=T))
top_genes <- rownames(dat)[unique(as.vector(o[1:n_genes, ]))]
dat <- dat[top_genes, ]

# cluster columns based on correlation distance
col_order <- hclust(dist(t(dat), method='euclidean'), method='ward.D2')$order
dat <- dat[, col_order]
row_order <- NULL
for (i in seq(ncol(dat))) {
  ro <- order(dat[, i], decreasing=T)
  row_order <- c(row_order, setdiff(ro[1:n_genes], row_order))
}
dat <- dat[row_order, ]

# specify annotation colors
unique_modules <- unique(gene_data$module[row_order])
annColors <- list(Function = viridis(4) %>%
                    setNames(c('Capacity control', 'Folding', 'Trafficking', 'Glycosylation')))
annColors$Function <- annColors$Function[intersect(names(annColors$Function), unique_modules)]
annData <- gene_data[, c('module'), drop=F] %>% setNames('Function')

# generate heatmap
pdf(file=paste0(fig_dir, '/', classVar, '_', model, '_topEach_heatmap.pdf'), width=4.5, height=5.5, onefile=F)
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






