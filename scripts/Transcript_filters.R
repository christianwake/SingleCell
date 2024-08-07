library('sys')
library('Seurat')
library('stringr')
library('pheatmap')
library('ggplot2')
#library('umap')
library('textshape')
library('dplyr')
library('biomaRt')
library('grid')
library('scales')
library('data.table')
library('VennDiagram')
#library('scuttle')

source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R')

if(interactive()){
  project <- '2021614_21-002'
  qc_name <- '2024-01-20'
  # qc_name <- '2023-Jan'
  # qc_name <- 'Both_celltypes'
  # sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/All_data.RDS')
  # covs_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/Sample_sheet.csv')
  # filter_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/QC_steps/RNA_filters.csv')
  batch_name <- 'Sample_Name'
  batch_value <- 'Su1_03_B_cells'
  # batch_name <- 'Date_sort'
  # batch_value <- '2021-11-09'
  # batch_value <- '2022-08-12'
  # #batch_value <- '2021-12-02'
  # base_dir <- '/hpcdata/vrc/vrc1_data/'
  
  #base_dir <- '/Volumes/VRC1_DATA/'
  base_dir <- '/hpcdata/vrc/vrc1_data/'
  # project <- '2022619_857.3b'
  # qc_name <- 'QC_SampleName_pass'
  sdat_file <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/data/All_data.RDS')
  #covs_file <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/Sample_sheet.csv')
  #filter_file <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/QC_steps/step3_cell_and_feature_filters.csv')
  filter_file <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/QC_steps/Transcript_filters.csv')
  
  #batch_name <- 'Sample_Name'
  #batch_value <- '15C225_w31_Probe+'
  gtf_file <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/data/gtf.RDS')
  out_rds <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/', batch_value, '/RNA_cell_filtered.RDS')
  out_pdf <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/', batch_value, '/MT_nFeatures.pdf')
  out_txt1 <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/', batch_value, '/Excluded_genes.txt')
  out_txt2 <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/', batch_value, '/Filtered_Ncells.txt')
} else {
  args = commandArgs(trailingOnly = TRUE)
  sdat_file <- args[1]
  filter_file <- args[2]
  out_txt1 <- args[3]
  out_txt2 <- args[4]
  out_rds <- args[5]
  out_tsv <- args[6]
}

filters <- read.table(filter_file, header = T, sep = ',')
analysis_only <- filters[which(filters$type == 'Analysis only'), ]
filters <- filters[which(filters$type != 'Analysis only'), ]
### Remove NA values
filters <- filters[which(!is.na(filters$value)),]
fraction_cells <- filters[which(filters$type %in% c('feature', 'gene', 'trancript') & filters$feature == 'fraction_cells'), 'value']
cell_min <- gsub('<', '', gsub('=', '', filters[which(filters$type == 'gene' & filters$feature == 'cells'), 'value']))
count_min <- gsub('<', '', gsub('=', '', filters[which(filters$type == 'gene' & filters$feature == 'counts'), 'value']))
key_words <- filters[which(filters$type %in% c('feature', 'gene', 'trancript') & filters$feature == 'pattern'), 'value']
exceptions <- filters[which(filters$type %in% c('feature', 'gene', 'trancript') & filters$feature == 'pattern_exceptions'), 'value']

if(length(key_words) == 0){
  key_words <- ''
}
if(length(exceptions) == 0){
  exceptions <- ''
}

### Record N filtered
filters$N <- NA
#### Read Seurat object
sdat <- readRDS(sdat_file)
print(row.names(sdat@assays$prot))
DefaultAssay(sdat) <- 'RNA'

exclude_gene_names <- c()
### Feature exclusion based on fraction of cells with at least one count
if(length(fraction_cells) > 0){
  fraction_cells <- as.numeric(gsub('<', '', fraction_cells))
  print(paste0('Excluding features present in fewer than ', fraction_cells, ' of cells.'))
  
  min_cells <- length(colnames(sdat)) * fraction_cells
  a1 <- Sys.time()
  a <- CreateSeuratObject(sdat@assays$RNA@data, assay = 'RNA', min.cells = min_cells)
  genes <- setdiff(row.names(sdat), row.names(a))
  rm(a)
  a2 <- Sys.time()
  a2-a1
  
  exclude_gene_names <- c(exclude_gene_names, genes)
  hist(rowSums(sdat@assays$RNA@data > 0)/length(row.names(sdat)), xlab = 'Fraction of cells with > 0 counts', main = 'N features by fraction of cells')
  ### Update filters file
  filters[which(filters$feature == 'fraction_cells'), 'N'] <- length(genes)
  filters[which(filters$feature == 'fraction_cells'), 'value'] <- 
    paste0(filters[which(filters$feature == 'fraction_cells'), 'value'], '(', min_cells, ')')
}

### Feature exclusion based on minimum cell counts
if(length(cell_min) > 0){
  print(paste0('Excluding features with fewer than ', cell_min, ' total cells'))
  cell_min <- as.numeric(cell_min)
  
  a <- CreateSeuratObject(sdat@assays$RNA@data, assay = 'RNA', min.cells = cell_min)
  genes <- setdiff(row.names(sdat), row.names(a))
  rm(a)

  exclude_gene_names <- c(exclude_gene_names, genes)
  hist(rowSums(sdat@assays$RNA@data > 0), xlab = 'N cells with > 0 counts', main = 'N features by fraction of cells')
  filters[which(filters$feature == 'cells'), 'N'] <- length(genes)
}
### Feature exclusion based on minimum read counts 
if(length(count_min) > 0){
  print(paste0('Excluding features with fewer than ', count_min, ' total counts.'))
  count_min <- as.numeric(count_min) 
  genes <- row.names(sdat)[which(rowSums(sdat@assays$RNA@data) <= count_min)]
  exclude_gene_names <- c(exclude_gene_names, genes)
  filters[which(filters$feature == 'counts'), 'N'] <- length(genes)
}

### Feature exclusion based on gene name pattern
all_names <- row.names(sdat@assays$RNA@counts)
key_words <- strsplit(key_words, ',')[[1]]
key_word_matches <- all_names[sapply(all_names, function(x) any(sapply(key_words, function(y) grepl(y, x))))]
key_word_matches <- sort(key_word_matches)
### Exclude those specified as exceptions
key_words <- strsplit(exceptions, ',')[[1]]
if (length(key_words) > 0) {
  print(paste0('Excluding features whose names match key word(s):  ', toString(key_words), ', with exceptions: ', exceptions))
  exception_gene_names <- key_word_matches[sapply(key_word_matches, function(x) any(sapply(key_words, function(y) grepl(y, x))))]
  key_word_matches <- key_word_matches[!key_word_matches %in% exception_gene_names]
  filters[which(filters$feature == 'pattern'), 'N'] <- length(key_word_matches)
  filters[which(filters$feature == 'pattern_exceptions'), 'N'] <- length(exception_gene_names)
}
exclude_gene_names <- unique(c(exclude_gene_names, key_word_matches))

print(filters)
filters[(length(row.names(filters)) +1),] <- c('gene', 'total_filtered', '', length(exclude_gene_names))
### Write
print(paste0('Excluding ', length(exclude_gene_names), ' genes.'))
write.table(exclude_gene_names, out_txt1, quote = F, sep = ',', row.names = F, col.names = F)
write.table(filters, out_tsv, quote = F, sep = '\t', row.names = F, col.names = T)

saveRDS(sdat, file = out_rds)

filters <- analysis_only
### Feature exclusion based on gene name pattern
key_words <- filters[which(filters$type %in% c('Analysis only') & filters$feature == 'pattern'), 'value']
exceptions <- filters[which(filters$type %in% c('Analysis only') & filters$feature == 'pattern_exceptions'), 'value']
if(length(key_words) == 0){
  key_words <- ''
}
if(length(exceptions) == 0){
  exceptions <- ''
}
all_names <- row.names(sdat@assays$RNA@counts)
key_words <- strsplit(key_words, ',')[[1]]
key_word_matches <- all_names[sapply(all_names, function(x) any(sapply(key_words, function(y) grepl(y, x))))]
key_word_matches <- sort(key_word_matches)
### Exclude those specified as exceptions
key_words <- strsplit(exceptions, ',')[[1]]
if (length(key_words) > 0) {
  print(paste0('Excluding features whose names match key word(s):  ', toString(key_words), ', with exceptions: ', exceptions))
  exception_gene_names <- key_word_matches[sapply(key_word_matches, function(x) any(sapply(key_words, function(y) grepl(y, x))))]
  key_word_matches <- key_word_matches[!key_word_matches %in% exception_gene_names]
  filters[which(filters$feature == 'pattern'), 'N'] <- length(key_word_matches)
  filters[which(filters$feature == 'pattern_exceptions'), 'N'] <- length(exception_gene_names)
}
egn <- c()
egn <- unique(c(egn, key_word_matches))
print(length(egn))
egn <- egn[which(!(egn %in% exclude_gene_names))]
write.table(egn, out_txt2, quote = F, sep = ',', row.names = F, col.names = F)
