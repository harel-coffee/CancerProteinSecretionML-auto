generate_allcancerdata_csv <- function(gene_list=NULL, sample_vars=NULL, file_name='allcancerdata', merge_coadread=FALSE, main_dir=NULL) {
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
  # merge_coadread    Merge the TCGA-COAD and TCGA-READ projects into one (TCGA-COADREAD)
  #                   DEFAULT = FALSE
  #
  # main_dir          Path to the CancerProteinSecretionML directory
  #                   e.g., 'usr/yourname/Documents/CancerProteinSecretionML'
  
  if (is.null(main_dir)) {
    stop('You must specify the path to the main CancerProteinSecretionML directory!')
  }
  
  # load required packages
  library(SummarizedExperiment)
  
  if (is.null(gene_list)) {
    # get list of PSP genes from txt file
    psp.data <- read.delim(file.path(main_dir, 'data', 'PSPgenes.txt'), sep = '\t')
    gene_list <- as.character(psp.data$name)
  }
  
  if (is.null(sample_vars)) {
    sample_vars <- c('CancerStatus','Project','TumorStage','TumorStageMerged','TumorStageBinary','OverallSurvival',
                     'Race','Gender','Barcode','Mutations')
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
  allcancer.data <- readRDS(file.path(main_dir, 'data', 'TCGA_data_tpm.rds'))
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
      
    } else if (s_var == 'Mutations') {  # add mutation data for all genes listed in "common_mut_genes.txt"
      
      # load mutation data (loads dataframe 'all_mut_map')
      all_mut_map <- readRDS(file.path(main_dir, 'data', 'TCGA_data_mutations.rds'))
      
      # load list of genes that are in the top-5 most mutated for at least one cancer type
      genes <- read.delim(file.path(main_dir, 'data', 'common_mut_genes.txt'), stringsAsFactors=F, header=F)
      keepRows <- rownames(all_mut_map) %in% genes$V1
      all_mut_map <- t(all_mut_map[keepRows, ])
      
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