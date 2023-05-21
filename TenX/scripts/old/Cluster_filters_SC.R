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

### Filters
### 1. Mostly empty clusters
### 2. Clusters with > 3 times the expression of a negative marker than the next highest cluster
### 3. Cells with 0 counts of any positive marker

args = commandArgs(trailingOnly=TRUE)
sdat_file <- args[1]
negative_markers <- args[2]
positive_markers <- args[3]
out_rds <- args[4]
out_txt <- args[5]
out_pdf <- args[6]
#sdat_file <- '/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/snakemake/results/Biolegend/Clustered.RDS'
#negative_markers <- 'AB-CD3,B-CD3,CD4, CD14'
#positive_markers <- 'CD19,CD20'
#out_rds <- '/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/snakemake/results/Biolegend/Filtered.RDS'
#out_txt <- '/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/snakemake/results/Biolegend/Filters.txt'
#out_pdf <- '/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/snakemake/results/Biolegend/Cluster_filters.pdf'

sdat <- readRDS(sdat_file)

fullN <- dim(sdat)[2]
txt <- as.data.frame(matrix(ncol = 2, nrow = 5))
colnames(txt) <- c('filter', 'n')
txt[, 'filter'] <- c('Low expression cluster', 'negative marker cluster', 'positive marker cell', 'remaining', 'total_filtered')

negative_markers <- trimws(strsplit(negative_markers, ',')[[1]])
positive_markers <- trimws(strsplit(positive_markers, ',')[[1]])

### Get the cell_ids  for each cluster
cluster_names <- levels(sdat$seurat_clusters)
### Cell membership of each cluster as a list
cluster_list <- lapply(cluster_names, function(cluster) names(sdat$seurat_clusters[which(sdat$seurat_clusters == cluster)]))
names(cluster_list) <- cluster_names

### Filter any largely empty clusters
cluster_means <- sapply(cluster_list, function(cluster) mean(as.numeric(sdat@assays$RNA@data[, cluster])))
small_clust <- names(cluster_means)[which(cluster_means == min(cluster_means))]
### Ratio of smallest cluster to second smallest
small_clust_frac <- cluster_means[order(cluster_means, decreasing = F)][1]/ cluster_means[order(cluster_means, decreasing = F)][2]
print(paste0('Cluster ', small_clust, ' has ', small_clust_frac, ' the average expression than the next highest cluster'))
### Add this cluster to list to filter, if it meets some relative and mean thresholds
if(small_clust_frac < 0.4 & cluster_means[small_clust] < 10){
  print(paste0('Removing cluster ', small_clust))
  to_remove <- colnames(sdat[,which(sdat$seurat_clusters %in% small_clust)])
  to_keep <- colnames(sdat)[which(!colnames(sdat) %in% to_remove)]
  
  ncells1 <- dim(sdat)[2]
  sdat <- subset(sdat, cells = to_keep)
  #sdat <- subset(sdat, subset = seurat_clusters != small_clust)
  print(paste0('N cells removed: ', (ncells1 - dim(sdat)[2]), ' of ', ncells1))
  txt[1, 'n'] <- (ncells1 - dim(sdat)[2])
  ### (Redo) Get the cell_ids  for each cluster
  cluster_names <- levels(sdat$seurat_clusters)[table(sdat$seurat_clusters) !=0]
  ### Cell membership of each cluster as a list
  cluster_list <- lapply(cluster_names, function(cluster) names(sdat$seurat_clusters[which(sdat$seurat_clusters == cluster)]))
  names(cluster_list) <- cluster_names
}

### Only keep those actually in the data
negative_markers <- negative_markers[which(negative_markers %in% row.names(sdat@assays$RNA@counts))]
### For each negative marker, find which cluster has the highest mean expression and by how much
negative_clusters <- c()
for(negative in negative_markers){
  cluster_means <- sapply(cluster_list, function(cluster) mean(as.numeric(sdat@assays$RNA[negative, cluster])))
  names(cluster_means) <- cluster_names
  ### cluster name with the largest negative marker value
  nclust <- names(cluster_means)[which(cluster_means == max(cluster_means))]
  ### Ratio of highest negative cluster to second highest
  nratio <- cluster_means[order(cluster_means, decreasing = T)][1]/ cluster_means[order(cluster_means, decreasing = T)][2]
  print(paste0('Cluster ', nclust, ' has ', nratio, ' times the ', negative, ' expression than the next highest cluster'))
  ### Add this cluster to list to filter, if it meets some relative and mean thresholds
  if(nratio > 3 & cluster_means[nclust] > 10){
    negative_clusters <- c(negative_clusters, nclust) 
  }
}
negative_clusters <- unique(negative_clusters)
if(length(negative_clusters) > 1){
  print(paste0('Removing cluster ', negative_clusters))
  to_remove <- colnames(sdat[,which(sdat$seurat_clusters %in% negative_clusters)])
  to_keep <- colnames(sdat)[which(!colnames(sdat) %in% to_remove)]
  
  ncells1 <- dim(sdat)[2]
  sdat <- subset(sdat, cells = to_keep)
  print(paste0('N cells removed by cluster marker in total: ', (ncells1 - dim(sdat)[2]), ' of ', ncells1))
  txt[2, 'n'] <- (ncells1 - dim(sdat)[2])
}

### Checking N cells to be removed for flat presence/absence of the markers
if(length(positive_markers > 0)){
  ncells1 <- dim(sdat)[2]
  sdat <- marker_cell_filter(sdat, output = 'Seurat', negatives = c(), positive_markers)
  above0 <- marker_cell_filter(sdat, output = 'cells', negative_markers, positives = c())
  txt[3, 'n'] <- (ncells1 - dim(sdat)[2])
  length(above0)  
}

DefaultAssay(sdat) <- 'RNA'
# VlnPlot(sdat, group.by = "RNA_clusters", features = "nFeature_RNA")
# VlnPlot(sdat, group.by = "RNA_clusters", features = 'MT_sum')

saveRDS(sdat, file = out_rds)
txt[which(txt$filter == 'remaining'), 'n'] <- dim(sdat)[2]
txt[which(txt$filter == 'total_filtered'), 'n'] <- fullN - dim(sdat)[2]
write.table(txt, out_txt, sep = ',', quote = F, row.names =F, col.names = T)

pdf(out_pdf)
for(m in negative_markers){
  print(VlnPlot(sdat, features = m) )
}
dev.off()

