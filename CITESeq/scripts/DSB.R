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
library('reticulate')
library('umap')
library('gridExtra')
library('cowplot')
library('logspline')

source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')

# source('/Volumes/VRC1_DATA/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
# source('/Volumes/VRC1_DATA/douek_lab/snakemakes/sc_functions.R')

if(interactive()){
  base_dir <- '/Volumes/VRC1_DATA/'
  base_dir <- '/hpcdata/vrc/vrc1_data/'
  
  project <- '2021614_21-002'
  qc_name <- 'DSB_by_sample'
  qc_name <- '2023-05-30'
  # batch_value <- '2021-12-02'
  # batch_name <- 'Date_sort'
  # batch_value <- 'Su7_03_Innate'
  batch_value <- 'Su1_01_B_cells'
  batch_name <- 'Sample_Name'
  hto_string <- 'C0251,C0252,C0253,C0254,C0255,C0256,C0257,C0258,C0259,C0260'
  isotype_controls <- 'C0090,C0091,C0092,C0095'
  
  # project <- '2022619_857.3b'
  # qc_name <- 'QC_first_pass'
  # batch_name <- 'Date_Sort'
  # batch_value <- '2022-02-22'
  # hto_string <- 'C0251,C0252,C0253,C0254,C0255,C0256,C0257,C0258,C0259,C0260'

  sdat_file <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/results/', qc_name, 
                      '/batches/', batch_value, '/Cell_filtered.RDS')
  br_rds <- paste0(base_dir,'/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Background.RDS')
  covs_file <- paste0(base_dir,'/douek_lab/projects/RNASeq/', project, '/Sample_sheet.csv')
  labels_file <- paste0(base_dir,'douek_lab/projects/RNASeq/', project, '/data/Cell_data.csv')
  out_file <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/', batch_value, '/DSB_normalized_data.RDS')
  
  }else{
  args <- commandArgs(trailingOnly=TRUE)
  
  sdat_file <- args[1]
  br_rds <- args[2]
  covs_file <- args[3]
  labels_file <- args[4]
  batch_value <- args[5]
  batch_name <- args[6]
  hto_string <- args[7]
  isotype_controls <- args[8]
  out_file <- args[9]
}

isotype_controls <- strsplit(isotype_controls, ',')[[1]]
### This file uses the "count" slot of the prot assay for DSB, and prints DSB results to a separate file (not Seurat)

### Read covariates
covs <- read.csv(covs_file, stringsAsFactors = F, header = T,
                check.names = F, 
                colClasses = 'character',
                sep = ',')
row.names(covs) <- covs$Sample_ID

### Seurat object does not have negative droplets
sdat <- readRDS(sdat_file)
### Background has only negative droplets
background <- readRDS(br_rds)
#row.names(sdat@assays$prot@data) <- gsub('-totalseq', '', row.names(sdat@assays$prot))
#row.names(background) <- gsub('_totalseq', '', row.names(background))
DefaultAssay(sdat) <- 'RNA'

print(all(row.names(background) == row.names(sdat@assays$prot@data)))

### Metadata is from dehashing, and has Negative droplets
# metadata <- read.table(labels_file, header = T, sep = ',')
# table(metadata$Assignment_CR, metadata$Assignment_simple)
# metadata <- metadata[which(metadata$Assignment != 'Doublet'),]
# table(metadata$Assignment_CR, metadata$Assignment_simple)

### Subset to the input batch
keep_cells <- colnames(sdat)[which(sdat@meta.data[, batch_name] == batch_value)]
sdat <- subset(sdat, cells = keep_cells)

### Subset background to this batch. This assumes that the first _ section of the colnames is the cellranger ID
covs <- covs[which(covs[, batch_name] == batch_value),]
background <- background[, colnames(background)[sapply(colnames(background), function(x) strsplit(x, '_')[[1]][1]) %in% unique(covs$CR_ID)]]

md <- sdat@meta.data

prot <- sdat@assays$prot@data
print(mean(colSums(prot)))

is.na(background)[1]
### Save as data frame. When read should be added to Seurat assay "data" slot
print(isotype_controls)
print(row.names(prot))
dsb <- DSB_once(prot, md, negative_mtx_rawprot = background, isotype_controls = isotype_controls)

print(mean(colSums(dsb)))
saveRDS(dsb, file = out_file)

