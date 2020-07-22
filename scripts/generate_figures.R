

library(dplyr)
library(tidyr)
library(ggplot2)



# Load PSP expression data
alldata <- readRDS('/Users/jonrob/Documents/PostDoc/CancerProteinSecretionML/data/allcancerdata_psp.rds')
alldata$Project <- sub('TCGA-', '', alldata$Project)  # remove "TCGA-" prefix on cancer type abbrevs







# Plot expression of one gene, grouped by cancer and subgrouped by a classVar
# ===========================================================================

# specify parameters
gene <- 'CRYAB'
classVar <- 'CancerStatus'  # 'mutTP53'
classLevels <- c('Primary solid Tumor', 'Solid Tissue Normal')  # c('TRUE', 'FALSE')
cancers <- 'all'  # use 'all' to include all possible cancers

# filter data
dat <- alldata %>% select(all_of(gene), all_of(classVar), 'Project')
dat <- dat[dat[, classVar] %in% classLevels, ]

# log-transform counts
dat[, gene] <- log2(dat[, gene] + 1)

# remove cancer types without at least 10 samples in each class level
if (length(cancers) == 1 && casefold(cancers) == 'all') { cancers <- unique(na.omit(dat$Project))}
tab <- table(dat[, c(classVar, 'Project')])[classLevels, ]
cancers <- colnames(tab)[colSums(tab >= 10) == 2]
dat <- dat[dat$Project %in% cancers, ]

# generate plot
ggplot(dat, aes_string(x="Project", y=gene, fill=classVar)) +
  geom_boxplot()











