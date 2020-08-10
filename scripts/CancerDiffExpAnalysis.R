CancerDiffExpAnalysis <- function(cancerType=NULL, classVar=NULL, classVarLevels=NULL, gene=NULL, main_dir=NULL) {
  # Run a differential expression analysis using edgeR.
  #
  # Input:
  #
  #   cancerType      Character of one cancer type, or list of multiple cancer types to analyze.
  #                   Specifying 'all' will analyze all available cancer types in the data.
  #                   If NULL, the function will return a list of all available cancer types.
  #
  #   classVar        Class variable by which to separate the samples. For example:
  #                   'CancerStatus', 'TumorStage', 'TumorStageMerged', 'TumorStageBinary',
  #                   'Race', 'Gender', or 'mutTP53'
  #                   If NULL, the function will return a list of all available class variables.
  #
  #   classVarLevels  The TWO values of the class variable by which the samples will be grouped.
  #                   Note that the fold-changes will be calculated as L2/L1.
  #                   For example, if classVar is 'CancerStatus', classVarLevels would likely be
  #                   c('Solid Tissue Normal', 'Primary solid Tumor').
  #                   If NULL, the function will return a list of all available levels for the 
  #                   chosen class variable.
  #
  #   gene            Leave empty (NULL) to perform a differential expression analysis.
  #                   If a gene name is specified, then the expression of that gene in the
  #                   two groups will be visualized with a plot.
  #
  #   main_dir        Path to the CancerProteinSecretionML directory.
  #                   e.g., 'usr/yourname/Documents/CancerProteinSecretionML'
  #

  if (is.null(main_dir)) {
    stop('You must specify the path to the main CancerProteinSecretionML directory!')
  }
  
  library(edgeR)
  library(ggplot2)
  library(SummarizedExperiment)
  
  # Load annotation data
  annotData <- readRDS(file.path(main_dir, 'data', 'allcancerdata_psp.rds'))
  
  # obtain/verify list of cancer types to analyze
  annotData$Project <- sub('TCGA-', '', annotData$Project)  # remove 'TCGA-' prefix from cancer codes
  allCancerTypes <- unique(annotData$Project)
  allCancerTypes <- allCancerTypes[!is.na(allCancerTypes)]
  if (is.null(cancerType)) {
    message('Returning list of available cancer types.')
    return(allCancerTypes)
  } else if ((length(cancerType) == 1) && (tolower(cancerType) == 'all')) {
    cancerType <- allCancerTypes
  } else if (!all(cancerType %in% allCancerTypes)) {
    invalidTypes <- setdiff(cancerType, allCancerTypes)
    stop('Invalid cancerType entries: ', paste(invalidTypes, collapse = ', '))
  }
  
  # obtain/verify class variables by which to analyze the data
  ind <- which(colnames(annotData) == 'CancerStatus')  # CancerStatus is the first classVar
  allClassVars <- colnames(annotData)[ind:ncol(annotData)]
  pspGenes <- colnames(annotData)[1:(ind - 1)]  # also retrieve PSP gene names
  if (is.null(classVar)) {
    message('Returning list of available class variables.')
    return(allClassVars)
  } else if (!(classVar == 'CancerStatus')) {
    # if CancerStatus is not the class variable, only include Primary solid Tumor samples
    keepSamples <- annotData$CancerStatus == 'Primary solid Tumor'
    annotData <- annotData[keepSamples, ]
  } else if (!any(allClassVars == classVar)) {
    stop('"', classVar, '" is not a valid class variable.')
  }
    
  # obtain/verify class variable levels by which to analyze the data
  # NOTE: the DE will be performed such that fold changes are level[2] / level[1]
  allClassVarLevels <- unique(annotData[, classVar])
  if (is.null(classVarLevels)) {
    message('Returning list of available class variable levels.')
    return(allClassVarLevels)
  } else if (!(length(classVarLevels) == 2)) {
    stop('ClassVarLevels must contain *exactly* 2 levels.')
  } else if (!all(classVarLevels %in% allClassVarLevels)) {
    invalidLevels <- setdiff(classVarLevels, allClassVarLevels)
    stop('Invalid level(s) in ClassVarLevels: ', paste(invalidLevels, collapse = ', '))
  }
  
  # extract barcodes
  L1_barcodes <- annotData$Barcode[annotData[, classVar] == classVarLevels[1]]
  L2_barcodes <- annotData$Barcode[annotData[, classVar] == classVarLevels[2]]
  
  
  # check if "gene" is specified
  if (!is.null(gene)) {
    
    # extract data for specified gene and class variable levels
    if (!any(colnames(annotData) == gene)) {
      stop('The specified gene "', gene, '" was not found in the dataset.')
    }
    
    # determine which cancer types have sufficient samples to be included
    keep.cancers <- NULL
    for (cType in cancerType) {
      cVals <- annotData[annotData$Project == cType, classVar]
      if ((sum(cVals %in% classVarLevels[1]) >= 10) && (sum(cVals %in% classVarLevels[2]) >= 10)) {
        keep.cancers <- c(keep.cancers, cType)
      }
    }
    
    # extract relevant subset of data
    sample.ind <- (annotData[, classVar] %in% classVarLevels) & (annotData$Project %in% keep.cancers)
    plot_data <- annotData[sample.ind, c('Project', classVar, gene)]
    plot_data[, gene] <- log(plot_data[, gene] + 1)  # log-transform TPM values
    plot_data[, classVar] <- reorder(plot_data[, classVar], plot_data[, classVar] == classVarLevels[2])
    
    # generate boxplot
    p <- ggplot(plot_data, aes_string(x=classVar, y=gene, color=classVar)) +
      geom_boxplot(outlier.size=0.5) +
      ylab(paste(gene,'log(TPM+1)')) +
      facet_grid(. ~ Project) + 
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
    
    return(p)
  }
  
  # iterate through cancer types
  for (cType in cancerType) {
    
    if (length(cancerType) > 1) {
      message('Processing cancer type: ', cType, ' (', which(cancerType == cType),
              ' of ', length(cancerType), ')')
    }
    
    # loads raw counts as SummarizedExperiment object named 'data'
    data <- readRDS(file.path(main_dir, 'data', 'TCGA_data_rawcounts', paste0(cType, '_raw_counts.rds')))
    
    # extract relevant data subset
    keepCols <- data$barcode %in% c(L1_barcodes, L2_barcodes)
    if (sum(data$barcode %in% L1_barcodes) < 10 || sum(data$barcode %in% L2_barcodes) < 10) {
      message('Insufficient samples in one or both groups - skipping cancer type.')
      next
    }
    count_data <- assay(data[, keepCols])  # count dataframe with barcodes as colnames
    geneNames <- rowData(data)$external_gene_name
    rm(data)  # no longer need this variable
    
    # create DGEList object from count_data
    y <- DGEList(counts=count_data)
    
    # specify groups such that fold-changes will be L2/L1
    groups <- factor(colnames(count_data) %in% L2_barcodes, labels=c('L1','L2'))
    design <- model.matrix(~groups)
    rownames(design) <- colnames(y)
    
    # filter low-count genes
    keep <- filterByExpr(y, design=design, min.count = 10)
    y <- y[keep, , keep.lib.sizes=F]
    geneNames <- geneNames[keep]
    
    # calculate normalization factors, and estimate dispersion
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design)
    
    # fit GLM and test for DE genes
    fit <- glmQLFit(y, design)
    qlf <- glmQLFTest(fit)
    
    # organize results table
    res.table <- qlf$table
    res.table$FDR <- p.adjust(res.table$PValue, 'BH')  # BH = Benjamini-Hochberg method
    res.table <- cbind(geneNames, res.table)
    
    # only extract results subset involving PSP genes
    res.table <- merge(data.frame(pspGenes), res.table, by.x='pspGenes', by.y='geneNames', all.x=T)
    
    # set default values for PSP genes excluded from the analysis due to low counts
    res.table$logFC[is.na(res.table$logFC)] <- 0
    res.table$PValue[is.na(res.table$PValue)] <- 1
    res.table$FDR[is.na(res.table$FDR)] <- 1
    
    # export results
    if (classVar %in% c('TumorStage', 'TumorStageMerged', 'TumorStageBinary')) {
      file_name_piece <- paste(cType, 'TumorStage', classVarLevels[1], classVarLevels[2], 'DEresults.txt', sep='_')
      file_name_piece <- gsub(' ', '', file_name_piece)  # remove spaces
    } else {
      file_name_piece <- paste(cType, classVar, 'DEresults.txt', sep='_')
    }
    res.filename <- file.path(main_dir, 'results', cType, file_name_piece)
    
    
    write.table(res.table, file=res.filename, quote=F, sep='\t', row.names=F)
    
  }
  
}



