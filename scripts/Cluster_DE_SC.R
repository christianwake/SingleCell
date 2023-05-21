library('sys')
library('Seurat')
library('stringr')
library('pheatmap')
library('ggplot2')
library('umap')
library('textshape')
library('dplyr')
library('biomaRt')
library('grid')
library('scales')
library('rsconnect')

source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

if(interactive()){
  project <- '2021618_galt'
  qc_name <- 'Go1'
  
  project <- '2022620_857.1'
  qc_name <- '2022-11-01'
  sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Filtered_clustered.RDS')
  exclude_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Excluded_genes.txt')
  out_rdata <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/App/Data.RData')
} else{
  args = commandArgs(trailingOnly=TRUE)
  sdat_file <- args[1]
  exclude_file <- args[2]
  out_rdata <- args[3] 
}

sdat <- readRDS(sdat_file)
### If the file exists and is not empty, read it
if(file.exists(exclude_file) & file.info(exclude_file)$size != 0){
  exclude_gene_names <- read.table(exclude_file)[,1] 
} else{
  exclude_gene_names <- c()
}

sdat <- ScaleData(sdat, assay = "RNA", features = row.names(sdat))

### Exclude genes, e.g. IG and ribosomal RNAs
rna_features <- row.names(sdat@assays$RNA@counts)[which(!row.names(sdat@assays$RNA@counts) %in% exclude_gene_names)]
rna_by_rna <- run_fam(sdat, 'RNA', 'RNA_clusters', 5, features = rna_features)
### Order
rna_by_rna <- rna_by_rna[order(rna_by_rna$p_val),]

saveRDS(rna_by_rna, out_rdata)

