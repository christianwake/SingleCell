library('sys')
library('Seurat')
library('stringr')
library('ggplot2')
#library('umap')
library('plyr')
library('dplyr')
library('grid')
library('scales')
library('data.table')

print(version)
print(.libPaths())
#.libPaths(.libPaths()[c(2,3,1)])
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R')

if(interactive()){
  # project <- '2021617_mis-c'
  
  # project <- '2021600_kristin'
  # qc_name <- 'Run2022-11-14'
  # sdat_step <- paste0('data/All_data.RDS')
  # norm_method <- 'TPM'
  
  project <- '2022620_857.1'
  qc_name <- '2022-11-16'
  sdat_step <- paste0('data/All_data.RDS')
  
  # project <- '2022619_857.3b'
  # qc_name <- 'QC_first_pass'
  # sdat_step <- paste0('data/All_data.RDS')
  
  sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/PostQC3.RDS')
  exclude_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/Excluded_genes.txt')
  gtf_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/gtf.RDS')
  out_rds <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/PostQC3.RDS')
} else{
  args = commandArgs(trailingOnly=TRUE)
  
  sdat_file <- args[1]
  exclude_file <- args[2]
  norm_method <- args[3]
  #gtf_file <- args[4]
  out_rds <- args[4]
}
if(norm_method == ''){
  norm_method <- 'Seurat'
}
print(norm_method)
sdat <- readRDS(sdat_file)
DefaultAssay(sdat) <- 'RNA'

### If the file exists and is not empty, read it
if(file.exists(exclude_file) & file.info(exclude_file)$size != 0){
  exclude_gene_names <- read.table(exclude_file)[,1] 
} else{
  exclude_gene_names <- c()
}
### Remove features

### Standard Seurat log normalization, should be used on cellranger-processed data like from TenX
if(norm_method == 'Seurat'){
  ### NormalizeData uses 'count' slot and overrides the 'data' slot
  print('Using Seurats NormalizeData function')
  sdat <- NormalizeData(sdat) 
} else if(norm_method == 'TPM'){
  print('Assuming input is TPM, taking log + 1')
  sdat@assays$RNA@data <- as.matrix(log(sdat@assays$RNA@counts + 1))
} #else if(norm_method == 'quminorm'){
  
#}

### 

print(out_rds)
saveRDS(sdat, file = out_rds)
