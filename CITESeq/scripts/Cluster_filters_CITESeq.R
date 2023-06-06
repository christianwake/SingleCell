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

source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

# source('/Volumes/VRC1_DATA/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
# source('/Volumes/VRC1_DATA/douek_lab/snakemakes/sc_functions.R')

### Filters
### 1. Clusters with very few reads
### 2. Clusters with > 3 times the expression of a negative marker than the next highest cluster
### 3. Cells with 0 counts of any positive marker

if(interactive()){
  # project <- '2021614_21-002'
  # qc_name <- '2023-Jan'
  # cseq_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/Test_subset.RDS')
  # cseq_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Clustered.RDS')
  # filter_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/QC_steps/cluster_filters.csv')
  # out_rds <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Filtered.RDS')
  # out_txt <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Filters.txt')
  # out_pdf <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Clustered.pdf')
 
  base_dir <- '/Volumes/VRC1_DATA/'
  #base_dir <- '/hpcdata/vrc/vrc1_data/' 
  project <- '2022619_857.3b'
  qc_name <- 'QC_first_pass'
  cseq_file <- paste0(base_dir,'/douek_lab/projects/RNASeq/', project, '/results/Test_subset.RDS')
  cseq_file <- paste0(base_dir,'/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Clustered.RDS')
  filter_file <- paste0(base_dir,'douek_lab/projects/RNASeq/', project, '/QC_steps/step5_cluster_filters.csv')
  out_rds <- paste0(base_dir,'douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Filtered.RDS')
  out_txt <- paste0(base_dir,'/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Filters.txt')
  out_pdf <- paste0(base_dir,'/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Clustered.pdf')
  
  
}else{
  args = commandArgs(trailingOnly=TRUE)
  cseq_file <- args[1]
  filter_file <- args[2]
  out_rds <- args[3]
  out_txt <- args[4]
  out_pdf <- args[5]  
}

filters <- read.table(filter_file, header = T, sep = ',')

cseq <- readRDS(cseq_file)
fullN <- dim(cseq)[2]

### To create a test set for interactive use
# csub <- subset(x = cseq, downsample = 1000)
# saveRDS(csub, file = '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/Test_subset.RDS')

### For each input filter (row in csv file)
to_remove <- c()
for(r in row.names(filters)){
  res <- filter_clusters(cseq, filters[r, 'data'], filters[r, 'threshold'], filters[r, 'features'], filters[r, 'value'])
  filters[r, 'cluster removed'] <- toString(res[['clusters']])
  filters[r, 'N cells'] <- res[['N']]
  to_remove <- c(to_remove, res[['cell_ids']])
}
to_remove <- unique(to_remove)
to_keep <- colnames(cseq)[which(!(colnames(cseq) %in% to_remove))]

pdf(out_pdf)

dat <- as.data.frame(table(cseq$prot_clusters))
colnames(dat) <- c('Cluster', 'N_cells')
bar_prot <-ggplot(data=dat, aes(x = Cluster, y = log(N_cells))) +
  geom_bar(stat="identity") + ggtitle('protein clusters')
bar_prot

dat <- as.data.frame(table(cseq$RNA_clusters))
colnames(dat) <- c('Cluster', 'N_cells')
bar_RNA <-ggplot(data=dat, aes(x = Cluster, y = log(N_cells))) +
  geom_bar(stat="identity") + ggtitle('RNA clusters')
bar_RNA

# filters[which(filters$threshold %in% c('negative', 'positive')),]
# ### Split on commas
# negative_markers <- trimws(strsplit(negative_markers, ',')[[1]])
# ### Only keep those actually in the data
# negative_markers <- negative_markers[which(negative_markers %in% row.names(cseq@assays[[datatype]]@counts))]
# for(m in negative_markers){
#   print(VlnPlot(cseq, features = m) )
# }

VlnPlot(cseq, group.by = "RNA_clusters", features = "nCount_RNA")
VlnPlot(cseq, group.by = "RNA_clusters", features = "nFeature_RNA")
VlnPlot(cseq, group.by = "RNA_clusters", features = 'MT_sum')
VlnPlot(cseq, group.by = "prot_clusters", features = "nCount_prot")
VlnPlot(cseq, group.by = "RNA_clusters", features = "nFeature_prot")
VlnPlot(cseq, group.by = "prot_clusters", features = 'MT_sum')
dev.off()

## Subset and save
cseq <- subset(cseq, cells = to_keep)
saveRDS(cseq, file = out_rds)
### Add remaining info to csv file and save
filters[dim(filters)[1] + 1, c('threshold', 'N cells')] <- c('Remaining', dim(cseq)[2])
filters[dim(filters)[1] + 1, c('threshold', 'N cells')] <- c('Total removed', fullN - dim(cseq)[2])
write.table(filters, out_txt, sep = ',', quote = T, row.names =F, col.names = T)



