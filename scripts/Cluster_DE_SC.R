library('sys')
library('Seurat')
#library('Seurat', lib.loc = .libPaths()[2])
library('stringr')
library('pheatmap')
library('ggplot2')
library('umap')
#library('textshape')
library('dplyr')
library('biomaRt')
library('grid')
library('scales')
library('rsconnect')

source('/data/vrc_his/douek_lab/snakemakes/sc_functions.R')
source('/data/vrc_his/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

if(interactive()){
  project <- '2021618_galt'
  qc_name <- 'Go1'
  
  project <- '2022620_857.1'
  qc_name <- '2022-11-01'
  sdat_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Filtered_clustered.RDS')
  exclude_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Excluded_genes.txt')
  out_rds <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Cluster_DE.RDS')
  out_tsv <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Cluster_DE.tsv')
} else{
  args = commandArgs(trailingOnly=TRUE)
  sdat_file <- args[1]
  exclude_file <- args[2]
  exclude_file2 <- args[3]
  gtf_file <- args[4]
  out_rds <- args[5]
  out_tsv <- args[6]
}

sdat <- readRDS(sdat_file)
### If the file exists and is not empty, read it
if(file.exists(exclude_file) & file.info(exclude_file)$size != 0){
  exclude_gene_names <- read.table(exclude_file)[,1] 
} else{
  exclude_gene_names <- c()
}
### If the gene-exclusion file exists and is not empty, read it
if(file.exists(exclude_file2) & file.info(exclude_file2)$size != 0){
  exclude_gene_names <- c(read.table(exclude_file)[,1])
} 

### Do direct matching from gtf
if(grepl('\\.gtf', gtf_file)){
  gtf <- read_gtf(gtf_file, feature_type = 'gene', atts_of_interest = c('gene_id', 'gene_name',  'gene_biotype'))
} else if(grepl('\\.RDS', gtf_file)){
  gtf <- readRDS(gtf_file)
} else{
  warning('No gtf file')
}
row.names(gtf) <- gtf$gene_id

sdat <- ScaleData(sdat, assay = "RNA", features = row.names(sdat))

### Exclude genes, e.g. IG and ribosomal RNAs
rna_features <- row.names(sdat@assays$RNA@layers$counts)[which(!row.names(sdat@assays$RNA@layers$counts) %in% exclude_gene_names)]
print('Beginning run_fam()')
de <- run_fam(sdat, 'RNA', 'RNA_clusters', 5, features = rna_features)
print('Done with run_fam()')
### Order
de <- de[order(de$p_val),]
### Add gene name
### Determine whether 'gene_name' or 'gene_id' better matches those in the seurat object
gtf_cols <- c('gene_name', 'gene_id')
de$gene_name <- gtf[row.names(de), 'gene_name']
de <- de[, c('gene_name', colnames(de)[which(colnames(de) != 'gene_name')])]
saveRDS(de, out_rds)
write.table(de, file = out_tsv, sep = '\t', quote = F, row.names = T)
