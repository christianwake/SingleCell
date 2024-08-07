library('sys')
library('Seurat')
library('stringr')
library('pheatmap')
library('ggplot2')
library('textshape')
library('dplyr')
library('biomaRt')
library('grid')
library('scales')
library('data.table')
library('VennDiagram')
library('scuttle')

source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R')

if(interactive()){
  # project <- '2021614_21-002'
  # qc_name <- '2023-March'
  # qc_name <- '2023-Jan'
  # qc_name <- 'Both_celltypes'
  # sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/All_data.RDS')
  # covs_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/Sample_sheet.csv')
  # filter_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/QC_steps/RNA_filters.csv')
  # batch_name <- 'Date_sort'
  # batch_value <- '2021-11-09'
  # batch_value <- '2022-08-12'
  # #batch_value <- '2021-12-02'
  # base_dir <- '/hpcdata/vrc/vrc1_data/'
  
  base_dir <- '/Volumes/VRC1_DATA/'
  #base_dir <- '/hpcdata/vrc/vrc1_data/'
  project <- '2022619_857.3b'
  qc_name <- 'QC_SampleName_pass'
  sdat_file <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/data/All_data.RDS')
  covs_file <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/Sample_sheet.csv')
  filter_file <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/QC_steps/step3_cell_and_feature_filters.csv')
  batch_name <- 'Sample_Name'
  batch_value <- '15C225_w31_Probe+'
  gtf_file <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/data/gtf.RDS')
  out_rds <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/', batch_value, '/RNA_cell_filtered.RDS')
  out_pdf <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/', batch_value, '/MT_nFeatures.pdf')
  out_txt1 <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/', batch_value, '/Excluded_genes.txt')
  out_txt2 <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/', batch_value, '/Filtered_Ncells.txt')
  
  setwd(paste0(base_dir, 'douek_lab/projects/RNASeq/2022619_857.3b'))
} else {
  args = commandArgs(trailingOnly=TRUE)
  sdat_file <- args[1]
  filter_file <- args[2]
  out_rds <- args[3]
}

filters <- read.table(filter_file, header = T, sep = ',')
downsample_RNA <- filters[which(grepl('downsample', tolower(filters$type)) & filters$feature == 'RNA'), 'value']
downsample_prot <- filters[which(grepl('downsample', tolower(filters$type)) & filters$feature == 'prot'), 'value']

sdat <- readRDS(sdat_file)
DefaultAssay(sdat) <- 'RNA'

### DOWNSAMPLE
if(length(downsample_RNA) > 0){
  print('Downsample RNA')
  ### Default is by batch, unless filter file has 'cell' in the 'value' line e.g. "downsample_by_cell"
  ds_by <- ifelse(grepl('cell', filters[which(grepl('downsample', tolower(filters$type)) & filters$feature == 'RNA'), 'type']), 'cell', 'batch')
  bycol <- ifelse(ds_by == 'cell', T, F)
  print(paste0('Down sampling each ', ds_by , ' by ', downsample_RNA, ' to RNA UMI total.'))
  downsampled <- downsampleMatrix(GetAssayData(object = sdat, assay = "RNA", slot = "counts"), prop = as.numeric(downsample_RNA), bycol = bycol)
  sdat@assays$RNA@counts <- downsampled
  sdat$nFeature_RNA <- colSums(sdat@assays$RNA@counts > 0)
  sdat$nCount_RNA <- colSums(sdat@assays$RNA@counts)
}
if(length(downsample_prot) > 0){
  print('Downsample protein')
  ### Default is by batch, unless filter file has 'cell' in the 'value' line e.g. "downsmaple_by_cell"
  ds_by <- ifelse(grepl('cell', filters[which(grepl('downsample', tolower(filters$type)) & filters$feature == 'prot'), 'type']), 'cell', 'batch')
  bycol <- ifelse(ds_by == 'cell', T, F)
  print(paste0('Down sampling each ', ds_by, ' by ', downsample_prot, ' to prot UMI total.'))
  ### Print mean nCount_prot before
  print(mean(colSums(sdat@assays$prot@data)))
  downsampled <- downsampleMatrix(GetAssayData(object = sdat, assay = "prot", slot = "data"), prop = as.numeric(downsample_prot), bycol = bycol)
  sdat@assays$prot@data <- downsampled
  sdat$nFeature_prot <- colSums(sdat@assays$prot@data > 0)
  sdat$nCount_prot <- colSums(sdat@assays$prot@data)
  ### Print mean nCount_prot after
  print(mean(colSums(sdat@assays$prot@data)))
}

print('Saving new seurat object')
saveRDS(sdat, file = out_rds)
