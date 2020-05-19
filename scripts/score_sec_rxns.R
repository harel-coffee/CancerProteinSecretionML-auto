generate_rxnScoreData <- function(geneRxnAssocFile=NULL, scoreMethod=median) {
  
  # load gene-reaction association matrix
  geneRxnAssoc <- read.delim(geneRxnAssocFile, row.names=1)
  
  # generate gene expression .csv file containing genes in the gene-reaction association matrix
  genes <- colnames(geneRxnAssoc)
  data_path <- dirname(geneRxnAssocFile)
  data <- generate_allcancerdata_csv(gene_list=genes, file_name=NULL)
  meta_data <- data[, which(colnames(data) == 'CancerStatus'):ncol(data)]
  tpm_data <- data[, 1:(which(colnames(data) == 'CancerStatus')-1)]
  
  # remove genes (columns) of geneRxnAssoc that are not in tpm_data
  geneRxnAssoc <- geneRxnAssoc[, (genes %in% colnames(tpm_data))]
  genes <- genes[genes %in% colnames(tpm_data)]
  if (any(duplicated(geneRxnAssoc))) {
    stop('Need to remove duplicated reaction(s)!')
  }
  
  # score reactions
  rxnScores <- data.frame(matrix(0, nrow=nrow(tpm_data), ncol=nrow(geneRxnAssoc))) %>% setNames(rownames(geneRxnAssoc))
  if (!is.function(scoreMethod)) {
    if (casefold(scoreMethod) == 'geomean') {
      scoreMethod = function(x) exp(sum(log(x[x > 0]))/length(x))
    }
  }
  for (i in seq(nrow(geneRxnAssoc))) {
    rxnScores[, i] <- apply(tpm_data[genes[geneRxnAssoc[i, ] > 0]], 1, scoreMethod)
  }
  
}