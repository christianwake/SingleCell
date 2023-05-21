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

source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R')

if(interactive()){
  project <- '2021600_kristin'
  covs_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021600_kristin/Cell_sheet.csv'
  out_mtx <-  '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021600_kristin/data/raw_mtx.rds'
  
  project <- '2022620_857.1'
  covs_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/Cell_sheet.csv'
  out_mtx <-  '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/data/raw_mtx.rds'
} else{
  args = commandArgs(trailingOnly=TRUE)
  project <- args[1]
  covs_file <- args[2]
  out_mtx <- args[3]
}

print(project)
print(covs_file)
print(out_mtx)
base_dir  <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/')
### Read in covariates
sample_info <- read.csv(covs_file, stringsAsFactors = F, header = T,
                 check.names = F, 
                 colClasses = 'character',
                 sep = ',')
### Remove those without FCID
FCID <- c('FCID', 'Flowcell_ID', 'flowcell_full')
FCID <- FCID[FCID %in% colnames(sample_info)][1]
sample_info <- sample_info[which(!is.na(sample_info[,FCID])),]
sample_info <- sample_info[which(sample_info[,FCID] != ''),]
### Define more useful columns from parts of existing columns
sample_info$cell_dir <- file.path(base_dir, 'data', 'count', sample_info$Cell_ID)
### Recreate tab_files from sample_sheet
sample_info$tab_files <- file.path(sample_info$cell_dir, 'gene_abundances.tab')
any(file.exists(sample_info$tab_files))
all(file.exists(sample_info$tab_files))

row.names(sample_info) <- sample_info$Cell_ID

print(paste0(toString(sample_info[which(!file.exists(sample_info$tab_files)), 'Cell_ID']), ' does not have tab output.'))
### Can't do anything without the count files, so remove cells (from the sample sheet) if their tab file is missing.
sample_info <- sample_info[which(file.exists(sample_info$tab_files)),]

mtx <- tab2mtx(sample_info$tab_files, sample_info$Cell_ID)
saveRDS(mtx, out_mtx)

