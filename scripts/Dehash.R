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
library('readxl')
library('WriteXLS')
library('Matrix')
library('pastecs')

source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/sample_sheet_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R')

#source('/Volumes/VRC1_DATA/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
#source('/Volumes/VRC1_DATA/douek_lab/wakecg/sample_sheet_functions.R')
#source('/Volumes/VRC1_DATA/douek_lab/snakemakes/sc_functions.R')
#source('/Volumes/VRC1_DATA/douek_lab/snakemakes/Utility_functions.R')

if(interactive()){
  #project <- '2021617_mis-c'
  # runs_dir <- '/hpcdata/vrc/vrc1_data/douek_lab/Runs/210511_A00243_0102_AHC3YJDSX2/wakecg/'
  # covs_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021617_mis-c/All_covariates.csv'
  # out_rds <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021617_mis-c/seurat_object.RDS'
  # out_labels <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021617_mis-c/labels.csv'
  # out_csv <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021617_mis-c/RNA_counts.csv'

  project <- '2021614_21-002'
  runs_dir <- '/hpcdata/vrc/vrc1_data/douek_lab/Runs/'
  covs_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/Sample_sheet.csv')
  out_labels <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/Cell_data.csv')
  out_metrics <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/', project, '/cellranger_count_metrics.tsv')
  quant_func <- 'multi'
  batch <- 'Date_sort'

  # project <- '2022619_857.3b'
  # quant_func <- 'multi'
  # batch <- 'Date_Sort'
  # runs_dir <- '/hpcdata/vrc/vrc1_data/douek_lab/Runs/'
  # covs_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/Sample_sheet.csv')
  # out_labels <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/Cell_data.csv')
  # out_metrics <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/BHLN7YDSX3/cellranger_count_metrics.tsv')
 # runs_dir <- '/Volumes/VRC1_DATA/douek_lab/Runs/'
 # covs_file <- paste0('/Volumes/VRC1_DATA/douek_lab/projects/RNASeq/', project, '/Sample_sheet.csv')
 # out_labels <- paste0('/Volumes/VRC1_DATA/douek_lab/projects/RNASeq/', project, '/data/Cell_data.csv')
 # out_metrics <- paste0('/Volumes/VRC1_DATA/douek_lab/projects/RNASeq/', project, '/BHLN7YDSX3/cellranger_count_metrics.tsv')
  quant_func <- 'count'
  batch <- 'Date_Sort'
  
  gtf_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/gtf.RDS')
}else{
  args <- commandArgs(trailingOnly=TRUE)
  
  covs_file <- args[1]
  gtf_file <- args[2]
  runs_dir <- args[3]
  project <- args[4]
  out_labels <- args[5]
  #out_metrics <- args[4]
}

demux_method <- 'seuratDemux'
project_dir <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project)
# project_dir <- paste0('/Volumes/VRC1_DATA/douek_lab/projects/RNASeq/', project)
### Read covariates
dat <- read.csv(covs_file, stringsAsFactors = F, header = T,
                 check.names = F, 
                 colClasses = 'character',
                 sep = ',')
### If the sheet is a replicate sheet rather than sample sheet
# if('Rep_ID' %in% colnames(dat)){
#   dat$Sample_ID <- dat$Rep_ID
# }
row.names(dat) <- dat$Sample_ID

### Do direct matching from gtf
if(grepl('\\.gtf', gtf_file)){
  gtf <- read_gtf(gtf_file, feature_type = 'gene', atts_of_interest = c('gene_id', 'gene_name',  'gene_biotype'))
} else if(grepl('\\.RDS', gtf_file)){
  gtf <- readRDS(gtf_file)
} else{
  warning('No gtf file')
}
row.names(gtf) <- gtf$gene_id

### Combine cellranger count output metrics
# if(quant_func == 'multi'){
#   runs_subdir <- '/outs/multi/count/raw_feature_bc_matrix/'
#   metrics <- lapply(unique(dat$CR_ID), function(x) as.data.frame(read.table(paste0(runs_dir, x, '/outs/per_sample_outs/', x, '/metrics_summary.csv'), sep = ',')[2,]))
#   metrics <- rbindlist(metrics)
#   colnames(metrics) <- gsub(' ', '_', unlist(read.table(paste0(runs_dir, dat$CR_ID[1], '/outs/per_sample_outs/', dat$CR_ID[1], '/metrics_summary.csv'), sep = ',')[1,]))
# } else if(quant_func == 'count'){
#   runs_subdir <- '/outs/raw_feature_bc_matrix/'
#   metrics <- lapply(unique(dat$CR_ID), function(x) as.data.frame(read.table(paste0(runs_dir, x, '/outs/metrics_summary.csv'), sep = ',')[2,]))
#   metrics <- rbindlist(metrics)
#   colnames(metrics) <- gsub(' ', '_', unlist(read.table(paste0(runs_dir, dat$CR_ID[1], '/outs/metrics_summary.csv'), sep = ',')[1,]))
# }
# metrics <- as.data.frame(t(metrics))
# colnames(metrics) <- unique(dat$CR_ID)
# write.table(metrics, out_metrics, sep = '\t', quote = F,  row.names = T, col.names = T)

demux_method <- 'seuratDemux'
label_list <- as.list(rep(NA, length(unique(dat$CR_ID))))
names(label_list) <- unique(dat$CR_ID)
### For each cellranger count output, read and demultiplex
for(id in unique(dat$CR_ID)){
  print(id)
  pdf_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/', id, '_hash.pdf')
  label_list[[id]] <- dehash_CR_ID(dat, id, runs_dir, project_dir, demux_method = demux_method, autoThresh = F, gtf = gtf, pdf_file)
}

#downdir <- "/hpcdata/vrc/vrc1_data/douek_lab/Runs/211123_A00243_0139_BHLKCFDSX2/count_output/CR1/outs/downsampled_prot/"
labels <- rbindlist(label_list)
labels <- as.data.frame(labels)

print(table(labels$Assignment_CR, labels$Assignment_simple))


row.names(labels) <- unlist(lapply(label_list, function(x) row.names(x)))
table(labels$Assignment)
labels[, 'cell_id'] <- row.names(labels)
nums <- as.data.frame(table(labels$Assignment))
colnames(nums) <- c('Assignment', 'N_demultiplexed')
table(nums$Assignment)
write.table(labels, out_labels, sep = ',', quote = F, row.names = F, col.names = T)
