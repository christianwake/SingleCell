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
library('lisi')

source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

if(interactive()){
  project <- '2021600_kristin'
  qc_name <- 'Run2023-05-14'
  sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Mapped.RDS')
  exclude_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Excluded_genes.txt')
  cluster_info <- 'RNA_clusters-0'
  out_rds <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/SubClusters/', cluster_info, '.RDS')
  
}else{
  args = commandArgs(trailingOnly=TRUE)
  sdat_file <- args[1]
  exclude_file <- args[2]
  cluster_info <- args[3]
  out_rds <- args[4]
}
  
### Genes to exclude
if(file.info(exclude_file)[, 'size'] > 0){
  exclude_gene_names <- read.table(exclude_file)[,1]
} else{
  exclude_gene_names <- c()
}

### Read Seurat object
sdat <- readRDS(sdat_file)
#DefaultAssay(sdat) <- assay

## Parse subcluster information
cluster_info <- strsplit(cluster_info, '-')[[1]]
subset_type<- cluster_info[1]
subset_clusters <- cluster_info[2:length(cluster_info)]

### Will ouput a Seurat object but with basically only the outputs of the function included within
sdat <- umap_cluster_subset(sdat, subset_type, subset_clusters, exclude_gene_names, 0.1, output_full = FALSE)

saveRDS(sdat, file = out_rds)
