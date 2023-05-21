library('sys')
library('readr')
library('WriteXLS')
library('EnhancedVolcano')

source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/DE_functions.R')

args = commandArgs(trailingOnly=TRUE)
if(interactive()){
  project <- '2021600_kristin'
  qc_name <- 'Run2022-11-14'
  rds <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Cluster_DE.RDS')
  outfile <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Cluster_DE.xls')
  pthresh <- 0.1
}else{
  args = commandArgs(trailingOnly=TRUE)
  rds <- args[1]
  outfile <- args[2]
  pthresh <- args[3]  
}

dat <- readRDS(rds)

dat <- dat[which(dat$p_val_adj < as.numeric(pthresh)),]

clusters <- unique(dat$cluster)
dats <- lapply(clusters, function(cluster) dat[which(dat$cluster == cluster),])

WriteXLS(x = dats, ExcelFileName = outfile, SheetNames = paste0('Cluster_', clusters), row.names = T)
