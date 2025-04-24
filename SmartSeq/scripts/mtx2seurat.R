### This script assumes you've already completed bcl2fastq, Trimmomatic, STAR and maybe BALDR
library('sys')
library('Seurat')
library('ggplot2')
library('dplyr')
library('viridis')
library('data.table')
#library('PKI')
#library('tinytex')
#library('biomaRt')
#library('harmony')
library('readxl')
library('stringr')
library('stringi')
#library('limma')
library('vegan')

source('/data/vrc_his/douek_lab/snakemakes/sc_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/Utility_functions.R')

if(interactive()){
  project <- '2021600_kristin'
  project <- '2024615_Boswell'
  #project <- '2022620_857.1'
  
  mtx_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/data/raw_mtx.RDS')
  covs_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/data/Covariates_QC_metrics.csv')
  gtf_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/data/gtf.RDS')
  out_seurat <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/data/All_data.RDS')
} else{
  args = commandArgs(trailingOnly=TRUE)
  project <- args[1]
  mtx_file <- args[2]
  covs_file <- args[3]
  gtf_file <- args[4]
  out_seurat <- args[5]
}
print(project)
print(mtx_file)
print(covs_file)
print(gtf_file)
print(out_seurat)
### Read in covariates
sample_info <- read.csv(covs_file, stringsAsFactors = F, header = T,
                        check.names = F, 
                        colClasses = 'character',
                        sep = ',')
row.names(sample_info) <- sample_info$Cell_ID


mtx <- readRDS(mtx_file)
print(paste0('Reading mtx file with ', dim(mtx)[1], ' features and ', dim(mtx)[2], ' cells'))
### Subset to those in the QC filtered sample file
mtx <- mtx[, sample_info$Cell_ID]

### gtf is only needed for "sum_chromosome(mtx, chr ='MT', gtf)"
gtf <- readRDS(gtf_file)

## Complete list of possible covariates to bring from sample_info to seurat object
#covs <- c('sample', 'animal', 'timepoint', 'plate_full', 'plate', 'flowcell', 'flowcell_plot', 'well', 'cell_id', 'Stim', 'CD4_CD8', 'Confidence',
#         'specificity', 'counts', 'after trimming', 'trim loss', 'mapped reads', 'mapped freq', 'counter', 'diversity')
#covs <- c(covs, colnames(trims), colnames(stars))
covs <- colnames(sample_info)
sdat <- mtx2seurat(mtx, sample_info, stars = NA, gtf, project = project, covs = covs)
### Because the smartseq quants are TPMs
sdat[['nCount_RNA']] <- NULL

saveRDS(sdat, out_seurat)


