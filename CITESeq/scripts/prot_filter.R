library('sys')
library('Seurat')
library('stringr')
library('pheatmap')
library('ggplot2')
#library('textshape')
library('dplyr')
library('biomaRt')
library('grid')
library('scales')
library('data.table')
library('VennDiagram')

source('/data/vrc_his/douek_lab/snakemakes/sc_functions.R')
source('/data/vrc_his/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/Utility_functions.R')

args = commandArgs(trailingOnly=TRUE)
batch_name <- args[1]
out_rds <- args[2]
out_pdf <- args[3]
out_txt <- args[4]
sdat_file <- args[5]
covs_file <- args[6]
filter_file <- args[7]
batch_files <- args[8:length(args)]

# batch_name <- 'Date_Sort'
# out_rds <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/RNA_cell_filtered.RDS')
# out_pdf <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/MT_nFeatures.pdf')
# out_txt <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/Filtered_Ncells.txt')
# sdat_file <- '/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/DSB_normalized_data.RDS'
# covs_file <- '/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/All_covariates.csv'
# filter_file <- '/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/filters.csv'
batch_files <- c('/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/2021-11-09/Filters_Ncells.txt',
                 '/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/2021-11-10/Filters_Ncells.txt',
                 '/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/2021-12-02/Filters_Ncells.txt')

filters <- read.table(filter_file, header = T, sep = ',')
nFeature_RNA <- filters[which(filters$type == 'cell' & filters$feature == 'nFeature_RNA'), 'value']
fraction_MT <- filters[which(filters$type == 'cell' & filters$feature == 'fraction_MT'), 'value']
key_words <- filters[which(filters$type == 'gene' & filters$feature == 'pattern'), 'value']
if(length(nFeature_RNA) == 0){
  nFeature_RNA <- '20,3(SD)'
}
#if(is.na(fraction_MT) | fraction_MT == ''){
if(length(fraction_MT) == 0){
  fraction_MT <- '0,0.4'
}
if(length(key_words) == 0){
  key_words <- ''
}

read_ncells <- function(batch_files, batch){
  batch_file <- batch_files[batch]
  tab <- read.table(batch_file, sep = ',')
  tab$batch <- batch
  tab$filter <- row.names(tab)
  tab <- tab[, c('batch', 'step', 'filter', 'N', 'N_by_step')]
}
### Get the batch name from the file path
names(batch_files) <- sapply(batch_files, function(batch_file) tail(strsplit(dirname(batch_file), '/')[[1]], n = 1))
### Read and combine the csv outputs from per-batch RNA_cell_filters.R
rna_filter <- rbindlist(lapply(names(batch_files), function(batch) read_ncells(batch_files, batch)))

### Read Seurat object
sdat <- readRDS(sdat_file)
DefaultAssay(sdat) <- 'RNA'

### Get array of lower,upper bound thresholds from the comma-delimited string. Options are flat values, or if in values of standard deviation, with a trailing (SD)
MT_thresh <- threshold_string_seurat(sdat, fraction_MT, 'fraction_MT', T)
NFeat_thresh <- threshold_string_seurat(sdat, nFeature_RNA, 'nFeature_RNA', T)


#### RNA plots again but by batch_name
p1 <- VlnPlot(sdat, features = "nFeature_RNA", group.by = 'Date_Sort', combine = T, log = T, pt.size = 0) +
p2 <- VlnPlot(sdat, features = 'fraction_MT', group.by = 'Date_Sort', combine = T, pt.size = 0) 
p3 <- VlnPlot(sdat, features = 'nFeature_dsb', group.by = 'Date_Sort', combine = T, pt.size = 0) 
p4 <- VlnPlot(sdat, features = 'nCount_dsb', group.by = 'Date_Sort', combine = T, pt.size = 0) 
  
pdf(out_pdf)
p1
p2
p3
p4

### Set up data frame for text summary
filter_names <- c('nFeature_RNA', 'fraction_MT')
filter_names <- paste(rep(filter_names, each = 2), c('low', 'high'), sep = '_')
txt <- as.data.frame(matrix(ncol = 3, nrow = length(filter_names) + 2))
row.names(txt) <- c('original', filter_names, 'remaining')
colnames(txt) <- c('step', 'N', 'N_by_step')
txt[, 'step'] <- 0:(length(filter_names)+1)

sdat_full <- sdat
NF_low <- colnames(sdat@assays$RNA@counts)[which(sdat$nFeature_RNA < NFeat_thresh[1])]
NF_high <- colnames(sdat@assays$RNA@counts)[which(sdat$nFeature_RNA > NFeat_thresh[2])]
MT_low <- colnames(sdat@assays$RNA@counts)[which(sdat$fraction_MT < MT_thresh[1])]
MT_high <- colnames(sdat@assays$RNA@counts)[which(sdat$fraction_MT > MT_thresh[2])]
txt['nFeature_RNA_low', 'N'] <- length(NF_low)
txt['nFeature_RNA_high', 'N'] <- length(NF_high)
txt['fraction_MT_low', 'N'] <- length(MT_low)
txt['fraction_MT_high', 'N'] <- length(MT_high)

plot.new()
my_venn(names = c('NFeature too low', 'NFeature too high', 'MT too high'), list(NF_low, NF_high, MT_high))
dev.off()

### Recording the number of cells removed at each step
sdat <- subset(sdat, subset = nFeature_RNA >=  NFeat_thresh[1])
print((dim(sdat_full)[2] - dim(sdat)[2]))
txt['nFeature_RNA_low', 'N_by_step'] <- (dim(sdat_full)[2] - dim(sdat)[2])
sdat$nFeature_RNA <- as.numeric(sdat$nFeature_RNA)
sdat <- subset(sdat, subset = nFeature_RNA <=  NFeat_thresh[2])
txt['nFeature_RNA_high', 'N_by_step'] <- (dim(sdat_full)[2] - txt[2, 'N_by_step'] - dim(sdat)[2])
sdat <- subset(sdat, subset = fraction_MT >=  MT_thresh[1])
txt['fraction_MT_low', 'N_by_step'] <- (dim(sdat_full)[2] - txt[2, 'N_by_step'] - txt[3, 'N_by_step'] - dim(sdat)[2])
sdat <- subset(sdat, subset = fraction_MT <= MT_thresh[2])
txt['fraction_MT_high', 'N_by_step'] <- (dim(sdat_full)[2] - txt[2, 'N_by_step'] - txt[3, 'N_by_step'] - txt[4, 'N_by_step']  - dim(sdat)[2])
plot(sdat$nFeature_RNA, sdat$fraction_MT, xlab = 'N genes', ylab = 'Sum MT TPM')
dev.off()

txt['original', 'N'] <- dim(sdat_full)[2]
txt['original', 'N_by_step'] <- NA
txt['remaining', 'N'] <- dim(sdat)[2]
txt['remaining', 'N_by_step'] <- NA
write.table(txt, out_txt, sep = ',', quote = F, row.names = T, col.names = T)
#sdat_removed <- subset(sdat_full, cells = colnames(sdat_full)[!colnames(sdat_full) %in% colnames(sdat)])

### Fraction kept by each sample. The first two lines are necessary in the case that an entire sample is filtered
# temp <- table(sdat$sample)
# temp[unique(sdat_full$sample)[which(!unique(sdat_full$sample) %in% unique(sdat$sample))]] <- 0
# temp/table(sdat_full$sample)[names(temp)]

saveRDS(sdat, file = out_rds)
