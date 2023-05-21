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
library('data.table')

source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

### Filters
### 1. Clusters with very few reads
### 2. Clusters with > x the expression of a negative marker than the next highest cluster
### 3. Cells with x counts of any positive marker

if(interactive()){
  # project <- '2021614_21-002'
  # qc_name <- 'Go1'
  # filter_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/QC_steps/cluster_filters.csv')

  project <- '2021600_kristin'
  qc_name <- 'Run2023-05-14'
  filter_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/QC_steps/step5_cluster_filters.csv')
  
  # project <- '2022620_857.1'
  # qc_name <- '2023-01-09'
  # qc_name <- '2022-12-09'
  # qc_name <- '2023-05-02' ### A few more cells salvaged
  
  filter_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/QC_steps/step5_cluster_filters.csv')
   
  gtf_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/gtf.RDS')
  sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/Test_subset.RDS')
  sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Clustered.RDS')
  out_rds <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Filtered.RDS')
  out_txt <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Filters.txt')
  out_pdf <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Clustered.pdf')
}else{
  args = commandArgs(trailingOnly=TRUE)
  sdat_file <- args[1]
  filter_file <- args[2]
  gtf_file <- args[3]
  out_rds <- args[4]
  out_txt <- args[5]
  out_pdf <- args[6] 
}

filters <- read.table(filter_file, header = T, sep = ',')

sdat <- readRDS(sdat_file)
fullN <- dim(sdat)[2]
### Find all cluster options from the seurat object and choose from our list, the first in our order of preference
cloptions <- c('RNA_clusters', 'clusters_seurat_integrated', 'clusters_harmony_integrated', 'cluster_unintegrated', 'seurat_clusters')
cloptions <- cloptions[cloptions %in% colnames(sdat@meta.data)[grepl('clusters', colnames(sdat@meta.data))]]
print(cloptions)
cluster_name <- cloptions[1]

### For each input filter (row in csv file)
filters$N <- NA
to_remove <- c()
for(r in row.names(filters)){
  res <- filter_clusters(sdat, filters[r, 'data'], filters[r, 'threshold'], filters[r, 'features'], filters[r, 'value'], cluster_name)
  to_remove <- c(to_remove, res[['cell_ids']])
  filters[r, 'N'] <- length(to_remove)
}

negative_markers <- filters[which(filters$threshold %in% c('negative')), 'features']
if('negative' %in% tolower(filters$threshold)){
  negative_markers <- trimws(strsplit(negative_markers, ',')[[1]])
} else{
  negative_markers <- c()
}

positive_markers <- filters[which(filters$threshold %in% c('positive')), 'features']
if('positive' %in% tolower(filters$threshold)){
  positive_markers <- trimws(strsplit(positive_markers, ',')[[1]])
} else{
  positive_markers <- c()
}

pdf(out_pdf)
if(length(c(negative_markers, positive_markers)) > 0){
  ### Gene_name -> gene_id conversion
  if(!all(c(negative_markers, positive_markers) %in% row.names(sdat))){
    ### Do directy matching from gtf
    if(grepl('\\.gtf', gtf_file)){
      gtf <- read_gtf(gtf_file, feature_type = 'gene', atts_of_interest = c('gene_id', 'gene_name',  'gene_biotype'))
    } else if(grepl('.RDS', gtf_file)){
      gtf <- readRDS(gtf_file)
    }
    
    ### For converting from gene_name to gene_id
    res <- gene_matching_v1(negative_markers, row.names(sdat), gtf, one2multi = 'exclude')
    negative_markers <- res[which(!is.na(res$gene_id)), 'gene_id']
    names(negative_markers) <- res[which(!is.na(res$gene_id)), 'gene_name']
    
    res <- gene_matching_v1(positive_markers, row.names(sdat), gtf, one2multi = 'exclude')
    positive_markers <- res[which(!is.na(res$gene_id)), 'gene_id']
    names(positive_markers) <- res[which(!is.na(res$gene_id)), 'gene_name']
  }
  
  ### List of cells per cluster
  cluster_names <- levels(sdat[[cluster_name]][,1])
  cluster_list <- lapply(cluster_names, function(cluster) row.names(sdat@meta.data[which(sdat[[cluster_name]] == cluster),]))

  names(cluster_list) <- cluster_names
  ### Only keep those actually in the data
  negative_value <- filters[which(filters$threshold == 'negative'), 'value']
  #negative_value <- '5(SD)'
  ### Determine what method to use for negative thresholding
  if(grepl('\\(', negative_value)){
    negative_method <- gsub('\\)', '', strsplit(negative_value, '\\(')[[1]][2])
    negative_value <- as.numeric(strsplit(negative_value, '\\(')[[1]][1])
  } else{ ### If nothing specied in parentheses, default is flat threshold
    negative_method <- 'relative'
  }
  negative_markers <- negative_markers[which(negative_markers %in% row.names(sdat))]
  negs <- c()
  for(i in 1:length(negative_markers)){
    print(VlnPlot(sdat, features = negative_markers[i]) + ggtitle(names(negative_markers[i])))
    
    means <- sapply(cluster_list, function(cluster) mean(as.numeric(sdat@assays$RNA[negative_markers[i], cluster])))
    names(means) <- cluster_names
    means <- means[order(means, decreasing = T)]
    if(negative_method == 'relative'){
      ##### 
      ### 1st most, 1st + 2nd most, 1st + 2nd + 3rd most, etc.
      clusters1 <- lapply(2:length(means), function(x) names(means)[1:(x-1)])
      ### Next most to compare to each of the above groups of clusters
      clusters2 <- names(means)[2:length(means)]
      comps <- sapply(1:(length(means)-1), function(ii) compare_cluster_mean(sdat, cluster_list, negative_markers[i], clusters1[[ii]], clusters2[ii], ratio = negative_value))
      if(any(comps)){
        negative_clusters <- clusters1[[which(comps == T)[1]]]
        print(paste0('Negative marker clusters: ', toString(negative_clusters)))
        ### Add cell IDs to negatives array
        negs <- c(negs, unlist(lapply(negative_clusters, function(c) cluster_list[[c]])))
      }
    } else if(negative_method == 'flat'){
       clusters <- names(means)[which(means > negative_value)]
       negs <- unlist(lapply(clusters, function(c) cluster_list[[c]]))
    }
  }
  
  filters[which(filters$threshold == 'negative'), 'N'] <- length(negs)
  
  to_remove <- c(to_remove, negs)
  to_remove <- unique(to_remove)
}
to_keep <- colnames(sdat)[which(!(colnames(sdat) %in% to_remove))]

quants <- c('nCount_RNA', 'nFeature_RNA', 'fraction_MT')
quants <- quants[quants %in% colnames(sdat@meta.data)]
for(quant in quants){
   print(VlnPlot(sdat, group.by = cluster_name, features = quant))
}

### Barplot size of clusters
dat <- as.data.frame(table(sdat[[cluster_name]]))
colnames(dat) <- c('Cluster', 'N_cells')
bar_RNA <- ggplot(data=dat, aes(x = Cluster, y = N_cells)) +
  geom_bar(stat="identity") + ggtitle('RNA clusters')
bar_RNA

dev.off()
filters <- rbind(filters, c('RNA', 'total filtered', '', '', length(unique(to_remove))))
filters <- rbind(filters, c('RNA', 'total remaining', '', '', length(colnames(sdat)) - length(unique(to_remove))))
sdat <- subset(sdat, cells = to_keep)
sdat[[cluster_name]] <- droplevels(sdat[[cluster_name]])


## Subset and save
saveRDS(sdat, file = out_rds)
### Add remaining info to csv file and save
write.table(filters, out_txt, sep = ',', quote = F, row.names =F, col.names = T)



