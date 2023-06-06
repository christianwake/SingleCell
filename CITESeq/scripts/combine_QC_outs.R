print(.libPaths())
print(version)
library('sys')
library('Seurat')
library('ggplot2')
library('dplyr')
library('viridis')
library('data.table')
library('patchwork')
library('PKI')
library('tinytex')
#library('tidyverse')
#library('pastecs')
library('reticulate')
library('gridExtra')
library('cowplot')
library('logspline')
library('WriteXLS')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

if(interactive()){
  project <- '2021614_21-002'
  qc_name <- '2023-Jan'
  batches <- c('2021-11-09', '2021-11-10', '2021-12-02', '2022-08-11', '2022-08-12')
  
  project <- '2022619_857.3b'
  qc_name <- 'QC_first_pass'
  batches <- c('2022-02-14', '2022-02-15', '2022-02-16', '2022-02-17', '2022-02-18', '2022-02-22', '2022-02-24', '2022-02-25')
  
  txt_out <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Filters_Ncells.txt')
  xls_out <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Filters_Ncells.xls')
  
  ex_out <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Excluded_genes.txt')
  filter_files <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/', batches, '/Filters_Ncells.txt')
  ex_files <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/', batches, '/Excluded_genes.txt')
  
  
  
}else{
  args <- commandArgs(trailingOnly=TRUE)
  
  txt_out <- args[1]
  xls_out <- args[2]
  #ex_out <- args[3]
  #all_files <- args[4:length(args)]
  #filter_files <- all_files[1:(length(all_files)/2)]
  #ex_files <- all_files[((length(all_files)/2)+1):length(all_files)]
  filter_files <- args[4:length(args)]
}
# ### Read gene exclusion files into list
# ex <- lapply(ex_files, function(ex_file) read.table(ex_file)[,1])
# ### Double check that each batch has the same exclusions
# if(length(unique(ex)) != 1){
#   warning('Batches have not excluded the same genes.')
# } else{
#   write.table(ex[[1]], ex_out, quote = F, sep = ',', row.names = F, col.names = F)
# }

read_ncells <- function(batch_files, batch){
  batch_file <- batch_files[batch]
  tab <- read.table(batch_file, sep = ',')
  tab$batch <- batch
  tab$filter <- row.names(tab)
  tab <- tab[, c('batch', 'step', 'filter', 'N_fail', 'N_by_step', 'freq_fail', 'freq_by_step')]
}
### Get the batch name from the file path
names(filter_files) <- sapply(filter_files, function(batch_file) tail(strsplit(dirname(batch_file), '/')[[1]], n = 1))
### Read and combine the csv outputs from per-batch RNA_cell_filters.R
dat_list <- lapply(names(filter_files), function(batch) read_ncells(filter_files, batch))
names(dat_list) <- names(filter_files)

### Sum
dats <- lapply(dat_list, function(dat) dat[, c('N_fail', 'N_by_step')])
summed <- dats[[1]]
for(i in 2:length(dats)){
  summed <- summed + dats[[i]]
}
### Add row.names
summed$filter <- dat_list[[1]][, 'filter']
summed <- summed[, c('filter', 'N_fail', 'N_by_step')]
summed$fraq_fail <- summed$N_fail /  summed[1,2]
summed$fraq_by_step <- summed$N_by_step /  summed[1,2]

out_list <- append(list(summed), dat_list)
names(out_list) <- c('Total', names(filter_files))
WriteXLS(out_list, ExcelFileName = xls_out, row.names = F)

### Simple csv
rna_filter <- as.data.frame(rbindlist(dat_list))
write.table(rna_filter, txt_out, sep = ',', quote = F, row.names = F, col.names = T)

