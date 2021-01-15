
# CancerProteinSecretionML: retrieve and process the UniProt PTM data for all human proteins
# *** CHANGE THIS PATH TO REPOSITORY PATH IN YOUR LOCAL COMPUTER ***
# path2secAI <- '/path/to/local/CancerProteinSecretionML/'
path2secAI <- '/Users/rasools/drive/projects/secAI/CancerProteinSecretionML/'
# packages used for this  
# sessionInfo()
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Big Sur 10.16

# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
# [1] stringr_1.4.0  UniprotR_2.0.2 biomaRt_2.45.6  
# 2021-01-14
####################################################################

library(biomaRt)
library(UniprotR)
library(stringr)

#first get all human genes
mart = useMart('ENSEMBL_MART_ENSEMBL')

#find corresponding protein(s) for each gene
mart2 = useDataset('hsapiens_gene_ensembl', mart)
allprot = getBM(mart = mart2,  attributes = c(
  'ensembl_gene_id',
  'uniprot_gn_id'),
  uniqueRows=TRUE)
allprot = allprot[which(allprot$uniprot_gn_id != ""),]

#find the status of proteins if you want to only keep reviewed one. It's a very time consuming step
#almost 1000ms for each entry.
st <- Sys.time()
st
allprot$status = NA
for (i in 1:dim(allprot)[1]) {
    remove(sta)
    sta = GetMiscellaneous(allprot$uniprot_gn_id[i])
    allprot$status[i] = sta$Status
}
en <- Sys.time()
en-st

write.table(allprot, file=paste(path2secAI, 'data/all_uniprot_proteins_status.txt', sep=''), row.names=F, col.names=T, sep='\t', quote=F)

#extract PTM information for proteins. Using loop helps to save retrieved data in case of connection instability and continue the loop from stopped entry
st <- Sys.time()
psim = allprot[which(allprot$status == 'reviewed'),]
for (i in 1:dim(psim)[1]) {
  remove(ptm)
  #retrieve and parse Post translation modification information
  ptm = GetPTM_Processing(psim$uniprot_gn_id[i])
  psim$N.gly[i] = str_count(ptm$Glycosylation, pattern = "N-linked")
  psim$O.gly[i] = str_count(ptm$Glycosylation, pattern = "O-linked")
  psim$ds[i] = str_count(ptm$Disulfide.bond, pattern = "DISULFID")
  psim$signal.peptide[i] = !is.na(ptm$Signal.peptide)
  
  #retrieve and parse general proteins information
  #sub cellular location
  #scl = GetSubcellular_location(psim$uniprot_gn_id[i])
  #psim$sub_cell_loc[i] = str_replace(scl$Subcellular.location..CC., 'SUBCELLULAR LOCATION: ', '')
  
  #miscellaneous info
  #mis = GetMiscellaneous(psim$uniprot_gn_id[i])
  #psim$evidence = mis$Protein.existence
  
  #protein name
  #taxa = GetNamesTaxa(psim$uniprot_gn_id[i])
  #psim$names = taxa$Protein.names
  #retrieve and parse proteins sub cellular location
}
psim[is.na(psim)] = 0
write.table(psim, file=paste(path2secAI, 'data/uniprot/uniprot_reviewed_info.txt', sep=''), row.names=F, col.names=T, sep='\t', quote=F)
en <- Sys.time()
en-st

