library('sys')
library('Seurat')
library('SeuratWrappers')
library('harmony')
library('stringr')
library('pheatmap')
library('ggplot2')
library('umap')
library('textshape')
library('dplyr')
library('biomaRt')
library('grid')
library('scales')
library('vegan')
library('data.table')

source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R')

args = commandArgs(trailingOnly=TRUE)

sdat_file <- args[1]
exclude_file <- args[2]
out_rds <- args[3]
# sdat_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021617_mis-c/results/Imputed.RDS'
# exclude_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021617_mis-c/results/Excluded_genes.txt'
# out_rds <- "/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021617_mis-c/results/Imputed_scaled.RDS"
# out_pdf <-  "/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021617_mis-c/results/Dropout.pdf"

sdat <- readRDS(sdat_file)
DefaultAssay(sdat) <- 'RNA'

### If the file exists and is not empty, read it
if(file.exists(exclude_file) & file.info(exclude_file)$size != 0){
  exclude_gene_names <- read.table(exclude_file)[,1] 
} else{
  exclude_gene_names <- c()
}

### NormalizeData uses 'count' slot and overrides the 'data' slot
sdat <- NormalizeData(sdat, assay = 'RNA')

sdat <- FindVariableFeatures(sdat, assay = 'RNA', selection.method = "vst")
#VariableFeaturePlot(sdat)
variable_genes <- sdat@assays$RNA@var.features
variable_genes <- variable_genes[which(!variable_genes %in% exclude_gene_names)]

### data slot is unchanged, scaled.data slot is added
### By default, it seems that ScaleData wants the raw 'counts', so if it is absent it will use 'data' but throw a warning
sdat <- ScaleData(sdat, assay = "RNA", vars.to.regress = c('MT_prop', 'nFeature_RNA'))
sdat <- RunPCA(sdat, verbose = F, assay = "RNA", features = variable_genes)
#sdat <- FindNeighbors(sdat, dims = 1:25)
#sdat <- FindClusters(sdat, resolution = 0.3)
#sdat <- RunUMAP(sdat, dims = 1:26, reduction = "pca", assay = "RNA")

saveRDS(sdat, file = out_rds)

