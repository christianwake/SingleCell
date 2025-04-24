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
### 1. Clusters with very few reads
### 2. Clusters with > 3 times the expression of a negative marker than the next highest cluster
### 3. Cells with 0 counts of any positive marker

args = commandArgs(trailingOnly=TRUE)
cseq_file <- args[1]
filter_file <- args[2]
out_rds <- args[3]
out_txt <- args[4]
out_pdf <- args[5]

# cseq_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/Test_subset.RDS'
# cseq_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/Go1/Clustered.RDS'
# filter_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/QC_steps/cluster_filters.csv'
# out_rds <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/Go1/Filtered.RDS'
# out_txt <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/Go1/Filters.txt'
# out_pdf <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/Go1/Clustered.pdf'

filters <- read.table(filter_file, header = T, sep = ',')
min_cells <- filters[which(filters$type == 'RNA' & filters$feature == 'min_cells'), 'value']
min_expr_flat <- filters[which(filters$type == 'RNA' & filters$feature == 'mean_expression_flat'), 'value']
min_expr_ratio <- filters[which(filters$type == 'RNA' & filters$feature == 'mean_expression_ratio'), 'value']
negative_markers <- filters[which(filters$type == 'RNA' & filters$feature == 'negative'), 'value']
positive_markers <- filters[which(filters$type == 'RNA' & filters$feature == 'positive'), 'value']
negative_markers <- trimws(strsplit(negative_markers, ',')[[1]])
positive_markers <- trimws(strsplit(positive_markers, ',')[[1]])

cseq <- readRDS(cseq_file)


datatype = 'RNA'
cluster_name = 'RNA_clusters'
### To create a test set for interactive use
# csub <- subset(x = cseq, downsample = 1000)
# saveRDS(csub, file = '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/Test_subset.RDS')

fullN <- dim(cseq)[2]
txt <- as.data.frame(matrix(ncol = 2, nrow = 5))
colnames(txt) <- c('filter', 'n')
txt[, 'filter'] <- c('Low expression cluster', 'negative marker cluster', 'positive marker cell', 'remaining', 'total_filtered')

### Get the cell_ids  for each cluster
cluster_names <- levels(cseq[[cluster_name]][,1])
### Cell membership of each cluster as a list
cluster_list <- lapply(cluster_names, function(cluster) row.names(cseq[[cluster_name]][which(cseq[[cluster_name]][,1] == cluster),, drop = F]))
names(cluster_list) <- cluster_names

### Filter any largely empty clusters
cluster_means <- sapply(cluster_list, function(cluster) mean(as.numeric(cseq@assays[[datatype]]@data[, cluster])))
cluster_sizes <- table(cseq[[cluster_name]])
### Ratio of smallest cluster to second smallest
small_clust_frac <- cluster_means[order(cluster_means, decreasing = F)][1]/ cluster_means[order(cluster_means, decreasing = F)][2]
small_clust <- names(small_clust_frac)
print(paste0('Cluster ', small_clust, ' has mean ', cluster_means[small_clust],', fractionally ', small_clust_frac, ' the average expression than the next highest cluster'))
### Add this cluster to list to filter, if it meets some relative and mean thresholds
if(small_clust_frac < 0.6 | cluster_means[small_clust] < 10){
  print(paste0('Removing cluster ', small_clust))
  to_remove <- colnames(cseq)[which(cseq[[cluster_name]][,1] %in% small_clust)]
  to_keep <- colnames(cseq)[which(!(colnames(cseq) %in% to_remove))]
  ncells1 <- dim(cseq)[2]
  cseq <- subset(cseq, cells = to_keep)
  #cseq <- subset(cseq, subset = prot_clusters != small_clust)
  print(paste0('N cells removed: ', (ncells1 - dim(cseq)[2]), ' of ', ncells1))
  txt[1, 'n'] <- (ncells1 - dim(cseq)[2])
  ### (Redo) Get the cell_ids  for each cluster
  cluster_names <- levels(cseq[[cluster_name]][,1])[table(cseq[[cluster_name]]) !=0]
  ### Cell membership of each cluster as a list
  cluster_list <- lapply(cluster_names, function(cluster) row.names(cseq[[cluster_name]][which(cseq[[cluster_name]][,1] == cluster),, drop = F]))
  names(cluster_list) <- cluster_names
}

### Only keep those actually in the data
negative_markers <- negative_markers[which(negative_markers %in% row.names(cseq@assays[[datatype]]@counts))]
### For each negative marker, find which cluster has the highest mean expression and by how much
negative_clusters <- c()
for(negative in negative_markers){
  cluster_means <- sapply(cluster_list, function(cluster) mean(as.numeric(cseq@assays[[datatype]][negative, cluster])))
  names(cluster_means) <- cluster_names
  ### cluster name with the largest negative marker value
  tclust <- names(cluster_means)[which(cluster_means == max(cluster_means))]
  ### Ratio of highest negative cluster to second highest
  nratio <- cluster_means[order(cluster_means, decreasing = T)][1]/ cluster_means[order(cluster_means, decreasing = T)][2]
  print(paste0('Cluster ', tclust, ' has ', nratio, ' times the ', negative, ' expression than the next highest cluster'))
  ### Add this cluster to list to filter, if it meets some relative and mean thresholds
  if(nratio > 3 & cluster_means[tclust] > 10){
    negative_clusters <- c(negative_clusters, tclust) 
  }
}

negative_clusters <- unique(negative_clusters)
if(length(negative_clusters) > 1){
  print(paste0('Removing cluster ', negative_clusters))
  to_remove <- colnames(cseq[,which(cseq[[cluster_name]][,1] %in% negative_clusters)])
  to_keep <- colnames(cseq)[which(!colnames(cseq) %in% to_remove)]
  
  ncells1 <- dim(cseq)[2]
  cseq <- subset(cseq, cells = to_keep)
  print(paste0('N cells removed by cluster marker in total: ', (ncells1 - dim(cseq)[2]), ' of ', ncells1))
  txt[2, 'n'] <- (ncells1 - dim(cseq)[2])
}

### Checking N cells to be removed for flat presence/absence of the markers
if(length(positive_markers > 0)){
  ncells1 <- dim(cseq)[2]
  cseq <- marker_cell_filter(cseq, output = 'Seurat', negatives = c(), positive_markers)
  above0 <- marker_cell_filter(cseq, output = 'cells', negative_markers, positives = c())
  txt[3, 'n'] <- (ncells1 - dim(cseq)[2])
  length(above0)  
}

# VlnPlot(cseq, group.by = "RNA_clusters", features = "nFeature_RNA")
# VlnPlot(cseq, group.by = "RNA_clusters", features = 'MT_sum')

saveRDS(cseq, file = out_rds)
txt[which(txt$filter == 'remaining'), 'n'] <- dim(cseq)[2]
txt[which(txt$filter == 'total_filtered'), 'n'] <- fullN - dim(cseq)[2]
write.table(txt, out_txt, sep = ',', quote = F, row.names =F, col.names = T)

pdf(out_pdf)
for(m in negative_markers){
  print(VlnPlot(cseq, features = m) )
}
dev.off()

