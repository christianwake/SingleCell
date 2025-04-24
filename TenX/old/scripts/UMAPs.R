library('sys')
library('Seurat')
library('stringr')
library('pheatmap')
library('ggplot2')
library('reticulate'); use_virtualenv("r-reticulate")
library('umap')
library('textshape')
library('dplyr')
library('biomaRt')
library('grid')
library('scales')

source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

args = commandArgs(trailingOnly=TRUE)

sdat_file <- args[1]
out_pdf <- args[2]

#sdat_file <- '/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/snakemake/results/Biolegend/Filtered_clustered.RDS'
#out_pdf <- '/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/snakemake/results/Biolegend/UMAP.pdf'

sdat <- readRDS(sdat_file)
DefaultAssay(sdat) <- 'prot'
config = umap.defaults
config$n_neighbors = 40
config$min_dist = 0.4

# run umap directly on dsb normalized values
dsb <- as.matrix(sdat@assays$prot@data)
ump <- umap(t(dsb), config = config)
umap_res <- as.data.frame(ump$layout)
colnames(umap_res) <- c("UMAP_1", "UMAP_2")

# save results dataframe 
df_dsb <- cbind(sdat@meta.data, umap_res, as.data.frame(t(dsb)))
##### Clusters - RNA_clusters, prot_clusters, p_dist_res.0.5, seurat_clusters
pdf(out_pdf)
### RNA UMAP
DimPlot(sdat, group.by = 'RNA_clusters', label = T, reduction = 'RNA_umap') + ggtitle('RNA UMAP, RNA clusters')
DimPlot(sdat, group.by = 'prot_clusters', label = T, reduction = 'RNA_umap') + ggtitle('RNA UMAP, protein clusters')

### Protein UMAP
DimPlot(sdat, group.by = 'prot_clusters', label = T, reduction = 'prot_umap') + ggtitle('Protein UMAP, protein clusters')
DimPlot(sdat, group.by = 'RNA_clusters', label = T, reduction = 'prot_umap') + ggtitle('Protein UMAP, RNA clusters')
dev.off()
