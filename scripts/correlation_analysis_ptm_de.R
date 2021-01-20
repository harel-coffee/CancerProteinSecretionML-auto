
# session info -----------------------------------------------------------
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Big Sur 10.16
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] pheatmap_1.0.12 dplyr_1.0.2     readxl_1.3.1    here_1.0.1     
# 
# loaded via a namespace (and not attached):
# [1] Rcpp_1.0.5         rstudioapi_0.11    magrittr_1.5       tidyselect_1.1.0   munsell_0.5.0      colorspace_1.4-1  
# [7] R6_2.4.1           rlang_0.4.8        tools_4.0.3        grid_4.0.3         gtable_0.3.0       ellipsis_0.3.1    
# [13] rprojroot_2.0.2    tibble_3.0.3       lifecycle_0.2.0    crayon_1.3.4       purrr_0.3.4        RColorBrewer_1.1-2
# [19] vctrs_0.3.4        glue_1.4.2         compiler_4.0.3     pillar_1.4.6       cellranger_1.1.0   generics_0.0.2    
# [25] scales_1.1.1       pkgconfig_2.0.3   

# loaded packages ---------------------------------------------------------

library(here)
library(readxl)
library(dplyr)
library(pheatmap)

# functions ---------------------------------------------------------------
# function for loading all sheets in excel file
read_excel_allsheets <- function(filename, tibble = TRUE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
 # function for saving figure as PFD
save_pheatmap_pdf <- function(x, filename, width = 2.2, height = 10) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# data -------------------------------------------------------------------
# read table of DE values as a list
tbls2 <- read_excel_allsheets(here("data", 'Table_S2.xlsx'))
tbls2$Legend <- NULL
# extract list of DE values as a data frame of logarithmic fold changes for each sample
nam <- names(tbls2)
genes <- tbls2[["ACC"]][["ensembl"]]
tbl_fc_1 <- as.data.frame(lapply(tbls2, "[", "FC_1"))
tbl_fc_2 <- as.data.frame(lapply(tbls2, "[", "FC_2"))
tbl_fc <- cbind(tbl_fc_1, tbl_fc_2)
colnames(tbl_fc) <- paste(nam, rep(c('FC1', 'FC2'), each=32), sep='_')
# load table of post translational modifications
ptm <- read.delim(here("data", "uniprot", "uniprot_reviewed_info.txt"), header = TRUE)
ptm <- ptm[which(ptm$ensembl_gene_id %in% genes),]
ptm <- ptm[match(genes, ptm$ensembl_gene_id), ]
ptm$ensembl_gene_id <- genes
remove(genes, tbls2, tbl_fc_1, tbl_fc_2)

# correlation analysis ---------------------------------------------------
# pearson
corel_data_pr <- cor(tbl_fc, ptm[, 4:6],  use = "complete.obs", method = "pearson")
corel_data_pr <- corel_data_pr[complete.cases(corel_data_pr), ]
P1 <- pheatmap(corel_data_pr, cluster_rows=FALSE, show_rownames=TRUE, show_colnames=TRUE, cluster_cols=FALSE, main = 'Pearson')
save_pheatmap_pdf(P1, "secAI_correlation_analysis_ptm_de_Pearson.pdf")

# spearman
corel_data_sp <- cor(tbl_fc, ptm[, 4:6],  use = "complete.obs", method = "spearman")
corel_data_sp <- corel_data_sp[complete.cases(corel_data_sp), ]
P2 <- pheatmap(corel_data_sp, cluster_rows=FALSE, show_rownames=TRUE, show_colnames=TRUE, cluster_cols=FALSE, main = 'Spearman')
save_pheatmap_pdf(P2, "secAI_correlation_analysis_ptm_de_Spearman.pdf")

