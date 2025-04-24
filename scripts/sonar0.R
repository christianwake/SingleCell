library('sys')
library('Seurat')
library('ggplot2')
library('dplyr')
library('viridis')
library('data.table')
library('patchwork')
#library('PKI')
library('tinytex')
library('reticulate')
library('umap')
library('gridExtra')
library('cowplot')
library('logspline')
library('readxl')
library('WriteXLS')
library('Matrix')
library('pastecs')

source('/data/vrc_his/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/sample_sheet_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/sc_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/dehash_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/Utility_functions.R')

if(interactive()){
  project <- '2021614_21-002'
  qc_name <- '2024-01-20'
  runs_dir <- '/data/vrc_his/douek_lab/Runs/'
  covs_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, 
                      '/Sample_sheet.csv')
  method <- 'Custom'
  calls_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, 
                        '/data/Dehash_CRint_calls_', method, '.tsv')
  out_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, 
                     '/results/', qc_name, '/For_SONAR.csv')
}else{
  args <- commandArgs(trailingOnly=TRUE)
  covs_file <- args[1]
  calls_file <- args[2]
  project <- args[3]
  method <- args[4]
  out_file <- args[5]
}

### Read covariates
covs <- read.csv(covs_file, stringsAsFactors = F, header = T,
                check.names = F, 
                colClasses = 'character',
                sep = ',')

if(!file.exists(calls_file)){
  calls_file <- paste0('data/Dehash_calls_', method, 'tsv')
}
calls <- read.table(calls_file, sep = '\t', header = T, quote = '')

### Remove Negatives and Doublets
calls <- calls[which(!is.na(calls$Assignment)), ]
calls <- calls[which(!(calls$Assignment %in% c('Negative', 'Ambiguous', 'Doublet'))), ]

### Remove unneeded columns
colnames(calls)
calls <- calls[, c('CR_ID', 'Assignment', 'Hashtag', 'bc')]
#calls <- calls[, c('CR_ID', 'Assignment')]

for(id in unique(calls$Assignment)){
  calls[which(calls$Assignment == id), c("Subject", 'Arm', 'Visit', 'Cell_subset', 'VDJ_index_name')] <- covs[which(covs$Sample_ID == id),c("Subject", 'Arm', 'Visit', 'Cell_subset', 'VDJ_index_name')]
}
### Remove those without VDJ index
calls <- calls[which(!(calls$VDJ_index_name == 'no VDJ')), ]

#calls$VDJ_dir <- paste0(getwd(), '/multi_output/', calls$CR_ID, '/outs/per_sample_outs/', calls$CR_ID, '/vdj_b/')
calls$VDJ_dir <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, 
                        '/multi_output/', calls$CR_ID, '/outs/per_sample_outs/', calls$CR_ID, 
                        '/vdj_b/')

all(dir.exists(calls$VDJ_dir))

write.table(calls, out_file, sep = ',', quote = F, row.names = F, col.names = T)
