generate_allcancerdata_csv <- function(gene_list=NULL, sample_vars=NULL, file_name='allcancerdata', count_type = 'TPM', merge_coadread = FALSE) {
  # load all cancer FPKM-UQ data, and save only a subset of the data in a .csv file.
  #
  # gene_list         List of genes (gene symbols or ensembl IDs) to include in the dataset.
  #                   DEFAULT = PSP genes.
  #
  # sample_vars       List of sample/class variables to include in the table.
  #                   DEFAULT = c('CancerStatus','Project','TumorStage',... etc. (see code)
  #
  # file_name         Name of output file to which data table will be saved. If NULL, then the
  #                   data will be returned as a dataframe without exporting to a file.
  #                   DEFAULT = 'allcancerdata_NEW.csv'
  #
  # count_type        Type of count data: 'FPKM', 'FPKM-UQ', or 'TPM'
  #                   DEFAULT = 'TPM'
  #
  # merge_coadread    Merge the TCGA-COAD and TCGA-READ projects into one (TCGA-COADREAD)
  #                   DEFAULT = FALSE
  #
  
  # load required packages
  library(SummarizedExperiment)
  
  if (is.null(gene_list)) {
    # get list of PSP genes from csv file
    psp.data <- read.csv('/Users/jonrob/Documents/PostDoc/TheCancerGenomeAtlas/RNAseq_analysis/gene_lists/secretome/human_secretory_machinery_NEW.csv',sep = ',')
    gene_list <- as.character(psp.data$gene_name)
    
    # get list of all protein-coding genes (VERY LARGE - NOT RECOMMENDED!)
    #prot.genes <- read.delim('/Users/jonrob/Documents/PostDoc/HumanProteinAtlas/data/proteinatlas_genes.txt')
    #gene_list <- as.character(prot.genes$protein_coding_genes)
    
    # get list of metabolic genes from file
    # met.genes <- read.delim('/Users/jonrob/Documents/PostDoc/TheCancerGenomeAtlas/RNAseq_analysis/gene_lists/metabolism/HMR2_genes.txt')
    # gene_list <- as.character(unique(met.genes$gene_abbrev))
  }
  
  if (is.null(sample_vars)) {
    sample_vars <- c('CancerStatus','Project','TumorStage','TumorStageMerged','TumorStageBinary','OverallSurvival',
                     'Race','Gender','Barcode','Mutations')#'HyperMut','HyperMutBinary')
    # HyperMut            Hypermutation status: 'Low', 'Hypermutant', or 'Ultrahypermutant'
    # HyperMutBin         Binary hypermutation status: 'Low', or 'Hypermutant'
  }
  
  # check that the first entry in sample_vars is "CancerStatus." This is required for the analysis scripts to know where the
  # gene names end, and where the sample variable names begin.
  CancerStatusInd <- which(sample_vars %in% 'CancerStatus')
  if (length(CancerStatusInd) == 0) {
    warning('"CancerStatus" must ALWAYS be included as the first entry in sample_vars. It has now been added to sample_vars.')
    sample_vars <- c('CancerStatus', sample_vars)
  } else if (CancerStatusInd != 1) {
    warning('"CancerStatus" must ALWAYS be the first entry in sample_vars. sample_vars has been rearranged such that this is now the case.')
    sample_vars = c('CancerStatus', sample_vars[!(sample_vars %in% 'CancerStatus')])
  }
  
  # load RNA-Seq data
  message(paste('Loading RNA-Seq (',count_type,') data... ',sep = ''), appendLF = FALSE)
  if (count_type == 'FPKM-UQ') {
    load('/Users/jonrob/Documents/PostDoc/TheCancerGenomeAtlas/RNAseqNEW_data/FPKM-UQ_RData/ALLCANCERS_fpkm_uq.RData')
  } else if (count_type == 'FPKM') {
    load('/Users/jonrob/Documents/PostDoc/TheCancerGenomeAtlas/RNAseqNEW_data/FPKM_RData/AllCancersCombined/ALLCANCERS_fpkm.RData')
  } else if (count_type == 'TPM') {
    load('/Users/jonrob/Documents/PostDoc/TheCancerGenomeAtlas/RNAseqNEW_data/FPKM_RData/AllCancersCombined/ALLCANCERS_tpm.RData')
  } else {
    stop('Invalid "count_type" input.')
  }
  message('Done.')
  
  # extract column data
  column.data <- colData(allcancer.data)
  
  # initialize sample data
  sample.data <- data.frame(matrix(nrow = nrow(column.data), ncol = 0))
  
  # loop through sample vars, and extract corresponding data from column data
  message('Processing data... ', appendLF = FALSE)
  for (s_var in sample_vars) {
    
    if (s_var == 'CancerStatus'){
      
      sample.data$CancerStatus <- column.data$definition
      
    } else if (s_var == 'Project') {
      
      sample.data$Project <- column.data$project_id
      if ( merge_coadread ) {
        sample.data$Project[sample.data$Project %in% c('TCGA-COAD','TCGA-READ')] = 'TCGA-COADREAD'
      }
      
    } else if (s_var == 'TumorStage') {
      
      sample.data$TumorStage <- column.data$tumor_stage
      
    } else if (s_var == 'TumorStageMerged') {
      
      ts = column.data$tumor_stage
      ts[ts %in% c('stage ia','stage ib')] <- 'stage i'
      ts[ts %in% c('stage iia','stage iib','stage iic')] <- 'stage ii'
      ts[ts %in% c('stage iiia','stage iiib','stage iiic')] <- 'stage iii'
      ts[ts %in% c('stage iva','stage ivb','stage ivc')] <- 'stage iv'
      sample.data$TumorStageMerged <- ts
      
    } else if (s_var == 'TumorStageBinary') {
      
      ts = column.data$tumor_stage
      # ts[ts %in% c('stage i','stage ia','stage ib','stage ii','stage iia','stage iib','stage iic')] <- 'stage i-ii'
      # ts[ts %in% c('stage iii','stage iiia','stage iiib','stage iiic','stage iv','stage iva','stage ivb','stage ivc')] <- 'stage iii-iv'
      ts[ts %in% c('stage i','stage ia','stage ib','stage ii','stage iia','stage iib','stage iic','stage iii','stage iiia','stage iiib','stage iiic')] <- 'stage i-iii'
      ts[ts %in% c('stage iv','stage iva','stage ivb','stage ivc')] <- 'stage iv'
      sample.data$TumorStageBinary <- ts
      
    } else if (s_var == 'Race') {
      
      sample.data$Race <- column.data$race
      
    } else if (s_var == 'Gender') {
      
      sample.data$Gender <- column.data$gender
      
    } else if (s_var == 'OverallSurvival') {
      
      is_alive <- column.data$vital_status %in% 'alive'
      overall_surv <- column.data$days_to_death
      overall_surv[is_alive] <- column.data$days_to_last_follow_up[is_alive]
      sample.data$OverallSurvival <- overall_surv
      
    } else if (s_var == 'Barcode') {
      
      sample.data$Barcode <- column.data$barcode
      
    # } else if (s_var == 'TP53') {  # this can be changed to any gene for which a mutation is desired
    # 
    #   mut_data <- read.delim(paste('/Users/jonrob/Documents/PostDoc/other_databases/cBioPortal/mutation_data/barcodes_mut_',s_var,'.txt',sep = ''))
    #   mut_barcodes <- as.character(mut_data$barcode)
    #   is_mut <- rep(0,length(column.data$barcode))
    #   for (bcode in mut_barcodes) {
    #     is_mut <- is_mut + grepl(bcode,column.data$barcode)
    #   }
    #   is_mut <- is_mut > 0
    #   sample.data$dummy <- is_mut
    #   names(sample.data)[length(sample.data)] <- paste('mut_', s_var, sep = '')
      
    } else if (s_var == 'Mutations') {
      
      # mut_files = list.files('/Users/jonrob/Documents/PostDoc/other_databases/cBioPortal/mutation_data/','barcodes_mut_')
      
      # load and transpose mutation data (loads dataframe 'all_mut_map')
      # load('/Users/jonrob/Documents/PostDoc/TheCancerGenomeAtlas/mutations/mutect2/tcga_mutation_map_cancergenes.RData')
      load('/Users/jonrob/Documents/PostDoc/TheCancerGenomeAtlas/mutations/mutect2/tcga_mutation_data.RData')
      
      # load list of genes that are in the top-5 most mutated for at least one cancer type
      genes <- read.delim('/Users/jonrob/Documents/PostDoc/TheCancerGenomeAtlas/mutations/mutect2/top5_mut_genes.txt',
                          stringsAsFactors=F, header=F)
      keepRows <- rownames(all_mut_map) %in% genes$V1  # only keep data for these genes
      all_mut_map <- t(all_mut_map[keepRows, ])  # orient barcodes as row names, and genes as column names
      
      # trim end of TCGA barcodes associated with mutation data (only keep patient and sample type identification portions)
      rownames(all_mut_map) <- substr(rownames(all_mut_map),1,16)
      
      # initialize mutation data matrix
      mut_data <- matrix(data = FALSE, nrow = nrow(sample.data), ncol = ncol(all_mut_map))
      colnames(mut_data) <- paste('mut', colnames(all_mut_map), sep = '')
      
      # iterate through TCGA barcodes associated with mutation data, and map to barcodes in sample.data
      bcodes <- rownames(all_mut_map)
      for (i in seq(nrow(all_mut_map))) {
        rowinds <- grep(bcodes[i],sample.data$Barcode)
        if (length(rowinds) == 1) {
          mut_data[rowinds,] <- mut_data[rowinds,] | all_mut_map[i,]
        } else if (length(rowinds) > 1) {
          for (rowind in rowinds) {
            mut_data[rowind,] <- mut_data[rowind,] | all_mut_map[i,]
          }
        }
      }
      
      # any samples that were not in the mutation dataset have an unknown mutation status
      barcode_truncated <- substr(sample.data$Barcode,1,16)
      mut_data[!(barcode_truncated %in% bcodes),] = 'unknown'
      
      # merge mutation data with sample.data
      sample.data <- cbind(sample.data,mut_data)

    } else if (grepl('HyperMut',s_var)) {
      
      # load file specifying which samples are hypermutated
      load('/Users/jonrob/Documents/PostDoc/collaborators/RaphaelFerreira/hypermutation/BarcodePatients/AllCancerHyperMutStatus.RData')
      
      # convert to binary variable if specified (change ultrahypermutant -> hypermutant)
      if (s_var == 'HyperMutBinary') {
        hMutData$burden[hMutData$burden %in% 'Ultrahypermutatant'] = 'Hypermutant'
        if (length(unique(hMutData$burden)) > 2) {
          stop('Problem in converting the hypermutation status to a binary variable.')
        }
      }
      
      # need to trim off end of TCGA barcode and keep only the patient and sample type portion, because the barcodes from 
      # the hypermutation data correspond to DNA samples instead of RNA samples, and will therefore not match
      hMutData$barcode <- substr(hMutData$barcode,1,16)
      
      # initialize data with NAs (samples without hypermutation status will remain NA)
      sampleHypMutStatus <- data.frame(x = matrix(NA,nrow = nrow(sample.data),ncol = 1))
      colnames(sampleHypMutStatus) <- s_var
      
      # iterate through TCGA barcodes associated with hypermutation data, and map to barcodes in sample.data
      for (i in seq(nrow(hMutData))) {
        # this is done in a loop because some barcodes have multiple matches
        rowinds <- grep(hMutData$barcode[i],sample.data$Barcode)
        sampleHypMutStatus[rowinds,1] <- hMutData$burden[i]
      }

      # merge hypermutation status with sample.data
      sample.data <- cbind(sample.data,sampleHypMutStatus)
      
    } else {
      message(paste('Sample variable "',s_var,'" not recognized.', sep = ''))
    }
  }
  
  # # prefix sample.data column names with "ClassVar_"
  # colnames(sample.data) <- paste('ClassVar_',colnames(sample.data),sep = '')
  message('Done.')
  
  # extract count data from SummarizedExperiment object
  message('Extracting data for selected genes... ', appendLF = FALSE)
  if ( any(grep('ENSG000',gene_list[1])) ) {
    rowdataind = 1
  } else {
    rowdataind = 2
  }
  keepRows <- rowData(allcancer.data)[,rowdataind] %in% gene_list
  allcancer.data.subset <- allcancer.data[keepRows,]
  count_data <- assay(allcancer.data.subset)
  
  # rename rows and columns
  rownames(count_data) <- rowData(allcancer.data.subset)[,rowdataind]
  colnames(count_data) <- seq(1,ncol(count_data))
  
  # transpose count matrix
  count_data <- t(count_data)
  
  # append sample data to count matrix as last columns
  count_data <- cbind(count_data,sample.data)
  message('Done.')
  
  # write to csv and/or rds
  if (!is.null(file_name)) {
    write.csv(count_data, file = paste(file_name,'.csv',sep=''), quote = TRUE)
    saveRDS(count_data, file = paste(file_name,'.rds',sep=''))
  } else {
    return(count_data)
  }
  
  
}