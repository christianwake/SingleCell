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

source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

if(interactive()){
  project <- '2021600_kristin'
  qc_name <- 'Run2023-05-14'
  assay <- 'RNA'
  cluster <- 'RNA_clusters-0-1'
  # project <- '2021614_21-002'
  # qc_name <- 'Go2'
  # assay <- 'RNA'
  # cluster <- 'predicted.celltype.l1-B-T'
  sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Mapped.RDS')
  exclude_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Excluded_genes.txt')
  out_rds <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Cluster_DE/One-one/test.RDS')
  out_tsv <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Cluster_DE/One-one/test.tsv')
} else{
  args = commandArgs(trailingOnly=TRUE)
  sdat_file <- args[1]
  exclude_file <- args[2]
  gtf_file <- args[3]
  out_rds <- args[4]
  out_tsv <- args[5]
  assay <- args[6]
  cluster <- args[7]
}

sdat <- readRDS(sdat_file)
DefaultAssay(sdat) <- assay

exclude_gene_names <- read.table(exclude_file)[,1]
### Exclude the IG and ribosomal RNAs
if(assay == 'RNA'){
  features <- row.names(sdat@assays$RNA@counts)[which(!row.names(sdat@assays$RNA@counts) %in% exclude_gene_names)]
} else{
  features <- NULL
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

### Parse cluster type and numbers from the input
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
if(cluster_name %in% colnames(sdat@meta.data)){
  print(levels(sdat[[cluster_name]][,1]))
  clusters <- unique(sdat[[cluster_name]][,1])
  if(all(cluster1 %in% clusters) & all(cluster2 %in% clusters)){
    Idents(sdat) <- sdat[[cluster_name]]
    
    de <- FindMarkers(sdat, assay = assay, slot = 'data', features = features, logfc.threshold = 0, ident.1 = cluster1, ident.2 = cluster2)
    ### Order
    de <- de[order(de$p_val),]
    ### Add gene name
    ### Determine whether 'gene_name' or 'gene_id' better matches those in the seurat object
    gtf_cols <- c('gene_name', 'gene_id')
    de$gene_name <- gtf[row.names(de), 'gene_name']
    de <- de[, c('gene_name', colnames(de)[which(colnames(de) != 'gene_name')])]
    ### Save
    saveRDS(de, out_rds)
    write.table(de, file = out_tsv, sep = '\t', quote = F, row.names = T)
    
  } else{
    print(paste0('Warning: Input cluster not present:', c(cluster1, cluster2)[which(!c(cluster1, cluster2) %in% clusters)]))
  }
} else {
  print('Warning: Input cluster name is not present.')
}
