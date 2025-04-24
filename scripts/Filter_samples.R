library('sys')
library('Seurat')
library('ggplot2')
library('dplyr')
library('viridis')
library('data.table')
library('patchwork')
#library('PKI')
#library('tinytex')
#library('dsb')
#library('tidyverse')
#library('pastecs')
#library('umap')
library('gridExtra')
library('cowplot')
#library('logspline')
library('readxl')
library('WriteXLS')
library('Matrix')

if(interactive()){
  #project <- '2021614_21-002'
  #qc_name <- 'DSB_by_sample'
  project <- '2021600_kristin'
  qc_name <- 'Run2025-03-03'
  sdat_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/All_data.RDS')
  filter_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/QC_steps/Sample_filters.csv')
  out_rds <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/PostQC2.RDS')
  out_csv <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Sample_sheet.csv')
  }else{
  args = commandArgs(trailingOnly=TRUE)
  
  sdat_file <- args[1]
  #batch_txt <- args[2] ### Currently not even used, but need to be snakefile input
  filter_file <- args[2]
  out_rds <- args[3]
  out_csv <- args[4]
  
}

### File should have columns named type, feature, value
filters <- read.table(filter_file, header = T, sep = ',')
filters <- filters[which(filters$type == 'batch_remove'), ]

sdat <- readRDS(sdat_file)
print(row.names(sdat@assays$prot))


for(filter in row.names(filters)){
  feature <- filters[filter, 'feature']
  value <- filters[filter, 'value']
  print(paste0('Filter: ', value, ' of ', feature))
  
  mdat <- sdat@meta.data
  if(feature %in% colnames(mdat)){
    ### Set NAs to strings to simplify (but only in the mdat, not in sdat)
    mdat[which(is.na(mdat[,feature])), feature] <- 'NA'
    ### Boolean
    exclude <- substr(value, 1,1) == '!'
    value <- gsub('!', '', value)
    ###
    values <- strsplit(value, ',')[[1]]
  
    if(all(values %in% names(table(mdat[,feature])))){
      if(exclude){
        cell_ids <- row.names(mdat[which(!mdat[,feature] %in% values),])
      } else{
        cell_ids <- row.names(mdat[which(mdat[,feature] %in% values),])
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
#write.table(sdat@meta.data, file = out_csv, quote = F, sep = ',')
