### This script assumes you've already completed bcl2fastq, Trimmomatic, STAR and maybe BALDR
library('sys')
library('Seurat')
library('ggplot2')
library('dplyr')
library('viridis')
library('data.table')
library('PKI')
#library('tinytex')
#library('biomaRt')
#library('harmony')
library('readxl')
library('stringr')
library('stringi')
#library('limma')
library('vegan')

source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/Utility_functions.R')

sheet <- NA
if(interactive()){
  mtx_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021600_kristin/data/raw_mtx0.rds'
  covs_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021600_kristin/data/Covariates_QC_metrics.csv'
  gtf_file <- '/hpcdata/vrc/vrc1_data/douek_lab/wakecg/data/GRCh38_reference.RDS'
  out_seurat <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021600_kristin/data/seurat1.rds'
  project <- '2021600_kristin'
} else{
  args = commandArgs(trailingOnly=TRUE)
  mtx_file <- args[1]
  covs_file <- args[2]
  project <- args[3]
  gtf_file <- args[4]
  out_seurat <- args[5]
}

### Read in covariates
sample_info <- read.csv(covs_file, stringsAsFactors = F, header = T,
                        check.names = F, 
                        colClasses = 'character',
                        sep = ',')
row.names(sample_info) <- sample_info$Cell_ID


mtx <- readRDS(mtx_file)
print(paste0('Reading mtx file with ', dim(mtx)[1], ' features and ', dim(mtx)[2], ' cells'))
mtx <- mtx[, sample_info$Cell_ID]

gtf <- readRDS(gtf_r_file)

## Complete list of possible covariates to bring from sample_info to seurat object
#covs <- c('sample', 'animal', 'timepoint', 'plate_full', 'plate', 'flowcell', 'flowcell_plot', 'well', 'cell_id', 'Stim', 'CD4_CD8', 'Confidence',
#         'specificity', 'counts', 'after trimming', 'trim loss', 'mapped reads', 'mapped freq', 'counter', 'diversity')
#covs <- c(covs, colnames(trims), colnames(stars))
covs <- colnames(sample_info)
sdat <- mtx2seurat(mtx, sample_info, stars = NA, gtf, project = project, covs = covs)

saveRDS(sdat, out_seurat)


