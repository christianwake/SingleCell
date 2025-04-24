
library('sys')
library('Seurat')
library('stringr')
library('pheatmap')
library('ggplot2')
#library('umap')
#library('textshape')
library('dplyr')
library('biomaRt')
library('grid')
library('scales')

source('/data/vrc_his/douek_lab/snakemakes/sc_functions.R')
source('/data/vrc_his/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

args = commandArgs(trailingOnly=TRUE)
cseq_file <- args[1]
out_rds <- args[2]
out_pdf <- args[3]
# cseq_file <- '/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/Go2/Filtered.RDS'
# out_rds <- '/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/Go2/Processed_data.RDS'
# out_pdf <- '/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/Go2/WNN.pdf'

# cseq <- readRDS('/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/Downsample_0.33/Processed_data.RDS')
# saveRDS(cseq$RNA.weight, '/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/Downsample_0.33/weights.RDS')
# rm(cseq)
# cseq <- readRDS('/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/Go2/Processed_data.RDS')
# saveRDS(cseq$RNA.weight, '/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/Go2/weights.RDS')
# 
# ds <- readRDS('/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/Downsample_0.33/weights.RDS')
# wall <-  readRDS('/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/Go2/weights.RDS')
# 
# ds <- as.data.frame(as.numeric(ds))
# ds$downsample <- 'True'
# colnames(ds) <- c('RNA_weight', 'downsample')
# wall <- as.data.frame(as.numeric(wall))
# wall$downsample <- 'False'
# colnames(wall) <- c('RNA_weight', 'downsample')
# dat <- rbind(ds, wall, make.row.names = T)
# 
# p <- ggplot(dat, aes(x=downsample, y=RNA_weight)) + 
#   geom_boxplot()

resolution <- 0.1
cseq <- readRDS(cseq_file)

print(colnames(cseq@meta.data))
print(names(cseq@graphs))
print(names(cseq@reductions))
print(cseq@neighbors)
# ### Will follow the tutorial which uses PCA of variable proteins, although we didn't previoulsy do PCA on prot, rather used UMAP directly on the prot features.
# DefaultAssay(cseq) <- 'prot'
# # we will use all ADT features for dimensional reduction
# # we set a dimensional reduction name to avoid overwriting the
# VariableFeatures(cseq) <- rownames(cseq[["prot"]])
# cseq <- NormalizeData(cseq, normalization.method = 'CLR', margin = 2) %>% 
#   ScaleData() %>% RunPCA(reduction.name = 'apca')
# # Identify multimodal neighbors. These will be stored in the neighbors slot, 
# # and can be accessed using cseq[['weighted.nn']]
# # The WNN graph can be accessed at cseq[["wknn"]], 
# # and the SNN graph used for clustering at cseq[["wsnn"]]
# # Cell-specific modality weights can be accessed at cseq$RNA.weight
# cseq <- FindMultiModalNeighbors(
#   cseq, reduction.list = list("pca", "apca"), 
#   dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
# )
# cseq <- RunUMAP(cseq, nn.name = "weighted.nn", reduction.name = "wnn_umap", reduction.key = "wnnUMAP_")
# ### Rename the last clusters before they are overridden
# cseq <- FindClusters(cseq, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)
# cseq$wnn_clusters <- cseq$seurat_clusters


##### DSB WNN tutorial

# set up dsb values to use in WNN analysis 
DefaultAssay(cseq) = "prot"
# hack seurat to use normalized protein values as a dimensionality reduction object.
prots <- rownames(cseq@assays$prot@data)
VariableFeatures(cseq) <- prots

print(summary(cseq$nCount_prot))
pdf('/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/ADT.pdf')
VlnPlot(cseq, features = 'nCount_dsb')
dev.off()
# run true pca to initialize dr pca slot for WNN 
cseq <- ScaleData(cseq, assay = 'prot', verbose = T)
cseq <- RunPCA(cseq, reduction.name = 'pdsb', reduction.key = 'pdsb_', features = VariableFeatures(cseq), approx = F, verbose = T)

# make matrix of norm values to add as dr embeddings
pseudo = t(cseq@assays$prot@data)
pseudo_colnames = paste('pseudo', 1:length(prots), sep = "_")
colnames(pseudo) = pseudo_colnames
# add to object 
cseq@reductions$pdsb@cell.embeddings = pseudo

print(names(cseq@reductions))

# run WNN 
cseq = FindMultiModalNeighbors(
  object = cseq,
  reduction.list = list("pca", "pdsb"),
  weighted.nn.name = "dsb_wnn", 
  knn.graph.name = "dsb_knn",
  #modality.weight.name = "dsb_weight",
  modality.weight.name = list('RNA.weight', 'prot.weight'),
  snn.graph.name = "dsb_snn",
  dims.list = list(1:50, 1:length(prots)), 
  verbose = T
)

print(table(Distances(object = cseq[["dsb_wnn"]]) %in% "NaN"))
cseq = FindClusters(cseq, graph.name = "dsb_knn", algorithm = 3, resolution = resolution,
                 random.seed = 1990,  verbose = T)
cseq$wnn_clusters <- cseq$seurat_clusters

cseq = RunUMAP(cseq, nn.name = "dsb_wnn", reduction.name = "dsb_wnn_umap", 
            reduction.key = "dsb_wnnUMAP_", seed.use = 1990, verbose = T)


print(colnames(cseq@meta.data))
print(names(cseq@graphs))
print(names(cseq@reductions))
print(cseq@neighbors)

# plot 
# p1 = Seurat::DimPlot(cseq, reduction = 'dsb_wnn_umap', group.by = 'dsb_knn',
#                      label = TRUE, repel = TRUE, label.size = 2.5, pt.size = 0.1) + 
#   theme_bw() + 
#   xlab('dsb protein RNA multimodal UMAP 1') + 
#   ylab('dsb protein RNA multimodal UMAP 2') + 
#   ggtitle('WNN using dsb normalized protein values')
# 
# p1


saveRDS(cseq, file = out_rds)


p1 <- DimPlot(cseq, reduction = 'dsb_wnn_umap', group.by = 'wnn_clusters', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p2 <- VlnPlot(cseq, features = "RNA.weight", group.by = 'wnn_clusters', sort = TRUE, pt.size = 0.1) +
  NoLegend()

pdf(out_pdf)
p1
p2
plot(cseq$RNA.weight, cseq$nCount_RNA)
boxplot(cseq$RNA.weight)
dev.off()

