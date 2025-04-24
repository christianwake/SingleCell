library('sys')
library('Seurat')
library('stringr')
library('pheatmap')
library('ggplot2')
#library('reticulate'); use_virtualenv("r-reticulate")
library('umap')
library('textshape')
library('dplyr')
library('biomaRt')
library('grid')
library('scales')
library('readxl')
library('cowplot')

source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

args = commandArgs(trailingOnly=TRUE)

sdat_file <- args[1]
exclude_file <- args[2]
negative_markers <- args[3]
gene_file <- args[4]
clust_resolution <- as.numeric(args[5])
out_rds <- args[6]
out_pdf <- args[7]

# sdat_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021618_galt/results/RNA_cell_filtered.RDS'
# exclude_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021618_galt/results/Excluded_genes.txt'
# negative_markers <- ''
# gene_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021618_galt/genes.xlsx'
# clust_resolution <- as.numeric('0.2')
# out_rds <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021618_galt/Clustered.RDS'
# out_pdf <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021618_galt/RNA_clusters.pdf'

### Defaults
if(is.na(clust_resolution) | clust_resolution == ''){
  clust_resolution <- 0.8
}
sdat <- readRDS(sdat_file)
DefaultAssay(sdat) <- 'RNA'

### If the file exists and is not empty, read it
if(file.exists(exclude_file) & file.info(exclude_file)$size != 0){
  exclude_gene_names <- read.table(exclude_file)[,1] 
} else{
  exclude_gene_names <- c()
}

negative_markers <- trimws(strsplit(negative_markers, ',')[[1]])

sdat <- NormalizeData(sdat)
sdat <- FindVariableFeatures(sdat, assay = 'RNA', selection.method = "vst")
#VariableFeaturePlot(sdat)

variable_genes <- sdat@assays$RNA@var.features
variable_genes <- variable_genes[which(!variable_genes %in% exclude_gene_names)]

sdat <- ScaleData(sdat, assay = "RNA", features = variable_genes)
sdat <- RunPCA(sdat, verbose = F, assay = "RNA", features = variable_genes)

sdat <- FindNeighbors(sdat, dims = 1:25)
sdat <- FindClusters(sdat, resolution = clust_resolution)
sdat <- RunUMAP(sdat, dims = 1:26, reduction = "pca", assay = "RNA")

### Save the RNA PCA, clusters and UMAP so that they can be overridden by the protein process
#sdat@reductions$RNA_pca <- sdat@reductions$pca

saveRDS(sdat, file = out_rds)

pdf(out_pdf)
DimPlot(sdat, group.by = 'seurat_clusters', label = T, reduction = 'umap') + ggtitle('RNA UMAP, RNA clusters')
FeaturePlot(sdat, feature = 'nFeature_RNA')
FeaturePlot(sdat, feature = 'nCount_RNA')
FeaturePlot(sdat, feature = 'MT_prop')
dev.off()

