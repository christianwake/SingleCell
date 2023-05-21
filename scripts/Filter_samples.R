library('sys')
library('Seurat')
library('ggplot2')
library('dplyr')
library('viridis')
library('data.table')
library('patchwork')
library('PKI')
library('tinytex')
library('dsb')
#library('tidyverse')
#library('pastecs')
library('umap')
library('gridExtra')
library('cowplot')
library('logspline')
library('readxl')
library('WriteXLS')
library('Matrix')

if(interactive()){
  project <- '2021614_21-002'
  #qc_name <- 'Both_celltypes'
  qc_name <- 'DSB_by_sample'
  sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/All_data.RDS')
  filter_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/QC_steps/Sample_filters.csv')
  out_rds <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/PostQC1.RDS')
}else{
  args = commandArgs(trailingOnly=TRUE)
  
  sdat_file <- args[1]
  #batch_txt <- args[2] ### Currently not even used, but need to be snakefile input
  filter_file <- args[2]
  out_rds <- args[3]
  
}

filters <- read.table(filter_file, header = T, sep = ',')
filters <- filters[which(filters$type == 'batch_remove'), ]

sdat <- readRDS(sdat_file)

#print(table(sdat$Date_sort, sdat$Sample_Name))
#print('Su7_03_B_cells' %in% sdat$Sample_Name)

for(filter in row.names(filters)){
  feature <- filters[filter, 'feature']
  value <- filters[filter, 'value']
  #print(paste0('Filter: ', value, ' of ', feature))
  #print('Su7_03_B_cells' %in% sdat$Sample_Name)
  
  if(feature %in% colnames(sdat@meta.data)){
    ### Boolean
    exclude <- substr(value, 1,1) == '!'
    value <- gsub('!', '', value)
    ###
    values <- strsplit(value, ',')[[1]]
    if(all(values %in% names(table(sdat[[feature]])))){
      if(exclude){
        cell_ids <- colnames(sdat)[which(!sdat@meta.data[, feature] %in% values)]
      } else{
        cell_ids <- colnames(sdat)[which(sdat@meta.data[, feature] %in% values)]
      }
      sdat <- subset(sdat, cells = cell_ids)
    }
    else{
      warning(paste0("Input value(s) ", value, " does not exist. Cannot filter."))
    }
  } else{
    warning(paste0('Input feature ', feature, ' does not exist. Cannot filter.'))
  }
}

saveRDS(sdat, file = out_rds)
