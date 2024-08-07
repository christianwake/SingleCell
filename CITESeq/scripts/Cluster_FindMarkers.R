library('sys')
library('Seurat')
library('stringr')
library('pheatmap')
library('ggplot2')
#library('umap')
library('textshape')
library('dplyr')
library('biomaRt')
library('grid')
library('scales')

source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

args = commandArgs(trailingOnly=TRUE)
cseq_file <- args[1]
exclude_file <- args[2]
out_rds <- args[3]
assay <- args[4]
cluster <- args[5]

# cseq_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/Go2/Mapped.RDS'
# exclude_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/Go2/Excluded_genes.txt'
# out_rds <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/Go2/Cluster_DE/One-one/test.RDS'
# assay = 'RNA'
# cluster = 'predicted.celltype.l1-B-T'

cseq <- readRDS(cseq_file)
DefaultAssay(cseq) <- assay
exclude_gene_names <- read.table(exclude_file)[,1]
### Exclude the IG and ribosomal RNAs
if(assay == 'RNA'){
  features <- row.names(cseq@assays$RNA@counts)[which(!row.names(cseq@assays$RNA@counts) %in% exclude_gene_names)]
} else{
  features <- NULL
}

cluster <- strsplit(cluster, '-')[[1]]
cluster_name <- cluster[1]
cluster1 <- c(cluster[2])
cluster2 <- c(cluster[3])
### Replace "T" with its constituent parts
if('T' %in% cluster1){
  cluster1 <- c(cluster1[!which(cluster1 != 'T')], 'CD4 T', 'CD8 T', 'other T')
}
if('T' %in% cluster2){
  cluster2 <- c(cluster2[!which(cluster2 != 'T')], 'CD4 T', 'CD8 T', 'other T')
}
  
### Check that the input data is present
if(cluster_name %in% colnames(cseq@meta.data)){
  print(levels(cseq[[cluster_name]][,1]))
  clusters <- unique(cseq[[cluster_name]][,1])
  if(all(cluster1 %in% clusters) & all(cluster2 %in% clusters)){
    Idents(cseq) <- cseq[[cluster_name]]
    
    cluster_de <- FindMarkers(cseq, assay = assay, slot = 'data', features = features, logfc.threshold = 0, ident.1 = cluster1, ident.2 = cluster2)
    saveRDS(cluster_de, out_rds)
  } else{
    print(paste0('Warning: Input cluster not present:', c(cluster1, cluster2)[which(!c(cluster1, cluster2) %in% clusters)]))
  }
} else {
  print('Warning: Input cluster name is not present.')
}
