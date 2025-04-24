library('sys')
library('readr')
library('WriteXLS')
library('EnhancedVolcano')

source('/data/vrc_his/douek_lab/snakemakes/Utility_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/DE_functions.R')

args = commandArgs(trailingOnly=TRUE)
if(interactive()){
  project <- '2021600_kristin'
  qc_name <- 'Run2023-05-14'
  rds <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Cluster_DE.RDS')
  outfile <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Cluster_DE.xls')
  pthresh <- 0.1
}else{
  args = commandArgs(trailingOnly=TRUE)
  rds <- args[1]
  outfile <- args[2]
  pthresh <- args[3]  
}

### Read the FindAllMarkers output RDS object
dat <- readRDS(rds)

### Subset each by input p-value threshold
p_val <- c('padj', 'p_val_adj')
p_val <- p_val[which(p_val %in% colnames(dat_list[[1]]))]
### Threshold based on adjusted p value
dat <- dat[which(dat[, p_val] < as.numeric(pthresh)),]
### Divide results by reference cluster
clusters <- unique(dat$cluster)
dats <- lapply(clusters, function(cluster) dat[which(dat$cluster == cluster),])

WriteXLS(x = dats, ExcelFileName = outfile, SheetNames = paste0('Cluster_', clusters), row.names = T)
