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

# source('/Volumes/VRC1_DATA/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
# source('/Volumes/VRC1_DATA/douek_lab/snakemakes/sc_functions.R')

### Filters
### 1. Clusters with very few reads
### 2. Clusters with > x times the expression of a negative marker than the next highest cluster
### 3. Cells with y counts of any positive marker

if(interactive()){
  project <- '2021614_21-002'
  qc_name <- '2024-01-20'
  base_dir <- '/data/vrc_his/' 
  sdat_file <- paste0(base_dir,'/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Clustered.RDS')
  filter_file <- paste0(base_dir,'douek_lab/projects/RNASeq/', project, '/QC_steps/cluster_filters.csv')
  out_rds <- paste0(base_dir,'douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Filtered.RDS')
  out_txt <- paste0(base_dir,'/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Filters.txt')
  out_pdf <- paste0(base_dir,'/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Clustered.pdf')
  gtf_file <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/data/gtf.RDS')
  
  
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
### Remove NA values
filters <- filters[which(!is.na(filters$value)),]

#sdat <- readRDS(sdat_file)
downsample_var <- 0.2
if(!interactive()){
  ### Read Seurat data
  sdat <- readRDS(sdat_file)
  print('Done reading Seurat object')
  ### Downsample so I can work interactiveley for testing
  assays <- c('RNA', 'prot')
  assays <- assays[which(assays %in% names(sdat@assays))]
  layers <- c('counts', 'data')
  layers <- layers[which(layers %in% names(sdat@assays$RNA@layers))]
  ssub <- DietSeurat(sdat, layers = layers, assays = assays)
  cells <- sample(x = colnames(ssub), size = (length(colnames(ssub)) * downsample_var), replace = F)
  
  ssub <- subset(ssub, cells = cells)
  saveRDS(ssub, file = gsub('.RDS', paste0('_DownSampledTo', downsample_var, '.RDS'), sdat_file))
  print('saved downsampled version for testing')
} else{
  sdat <- readRDS(gsub('.RDS', paste0('_DownSampledTo', downsample_var, '.RDS'), sdat_file))
}
fullN <- dim(sdat)[2]

### Do directy matching from gtf
if(grepl('\\.gtf', gtf_file)){
  gtf <- read_gtf(gtf_file, feature_type = 'gene', atts_of_interest = c('gene_id', 'gene_name',  'gene_biotype'))
} else if(grepl('.RDS', gtf_file)){
  gtf <- readRDS(gtf_file)
}

pdf(out_pdf)
### For each input filter (row in csv file)
to_remove <- c()
for(r in row.names(filters)){
  res <- filter_clusters(sdat, filters[r, 'data'], filters[r, 'threshold'], 
                         filters[r, 'features'], filters[r, 'value'], gtf)
  filters[r, 'cluster removed'] <- toString(res[['clusters']])
  filters[r, 'N cells'] <- res[['N']]
  to_remove <- c(to_remove, res[['cell_ids']])
}
to_remove <- unique(to_remove)
to_keep <- colnames(sdat)[which(!(colnames(sdat) %in% to_remove))]

if('RNA' %in% names(sdat@assays)){
  VlnPlot(sdat, group.by = "RNA_clusters", features = "nFeature_RNA") + xlab('RNA clusters')
  VlnPlot(sdat, group.by = "RNA_clusters", features = 'MT_sum') + xlab('RNA clusters')
}

if('prot' %in% names(sdat@assays)){
  VlnPlot(sdat, group.by = "prot_clusters", features = "nFeature_prot") + xlab('protein clusters')
  VlnPlot(sdat, group.by = "prot_clusters", features = 'MT_sum') + xlab('protein clusters')
}
dev.off()

## Subset and save
sdat <- subset(sdat, cells = to_keep)
saveRDS(sdat, file = out_rds)
### Add remaining info to csv file and save
filters[dim(filters)[1] + 1, c('threshold', 'N cells')] <- c('Remaining', dim(sdat)[2])
filters[dim(filters)[1] + 1, c('threshold', 'N cells')] <- c('Total removed', fullN - dim(sdat)[2])
write.table(filters, out_txt, sep = '\t', quote = F, row.names = F, col.names = T)



