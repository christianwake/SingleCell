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
#library('logspline')
library('WriteXLS')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

if(interactive()){
  project <- '2021614_21-002'
  # qc_name <- '2023-Jan'
  # batches <- c('2021-11-09', '2021-11-10', '2021-12-02', '2022-08-11', '2022-08-12')
  qc_name <- '2024-01-20'
  batches <- c('2021-11-09', '2021-11-10', '2021-12-02', '2022-08-11', '2022-08-12')
  batches <- c('Su13_03_Innate', 'Su19_3_Innate', 'Su9_03_B_cells')
  
  # project <- '2022619_857.3b'
  # qc_name <- 'QC_first_pass'
  # batches <- c('2022-02-14', '2022-02-15', '2022-02-16', '2022-02-17', '2022-02-18', '2022-02-22', '2022-02-24', '2022-02-25')
  
  txt_out <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Cell_filtered.txt')
  txt_out2 <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Batches_cells_remaining.txt')
  xls_out <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Cell_filtered.xls')
  
  ex_out <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Excluded_genes.txt')
  filter_files <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, 
                         '/results/', qc_name, '/batches/', batches, '/Cell_filtered.txt')
  ex_files <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, 
                     '/results/', qc_name, '/', batches, '/Excluded_genes.txt')
  
  
  
}else{
  args <- commandArgs(trailingOnly=TRUE)
  txt_out <- args[1]
  txt_out2 <-args[2]
  xls_out <- args[3]
  pdf_out <- args[4]
  #ex_out <- args[3]
  #all_files <- args[4:length(args)]
  #filter_files <- all_files[1:(length(all_files)/2)]
  #ex_files <- all_files[((length(all_files)/2)+1):length(all_files)]
  filter_files <- args[5:length(args)]
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
### Read and combine the csv outputs from per-batch Cell_filters.R
dat_list <- lapply(names(filter_files), function(batch) read_ncells(filter_files, batch))
names(dat_list) <- names(filter_files)

### Sum
dats <- lapply(dat_list, function(dat) dat[, c('N_fail', 'N_by_step')])
summed <- dats[[1]]
if(length(dats) > 1){
  for(i in 2:length(dats)){
    summed <- summed + dats[[i]]
  }
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


sample_dat <- t(sapply(names(dats), function(su) 
  dats[su][[1]][c('remaining', 'original'), 'N_fail']))
colnames(sample_dat) <- c('remaining', 'original')
sample_dat <- sample_dat[order(sample_dat[, 'remaining']), ]

write.table(sample_dat, txt_out2, sep = ',', quote = F, row.names = T, col.names = T)

pdf(pdf_out)
plot_dat <- as.data.frame(t(as.data.frame(lapply(names(dats), 
                     function(su) dats[[su]][c('original', 'remaining'), 'N_fail']))))
colnames(plot_dat) <- c('original', 'remaining')
row.names(plot_dat) <- names(dats)
plot_dat$batch <- names(dats)
plot_dat$N_filtered <- plot_dat$original - plot_dat$remaining
sample_order <- plot_dat[order(plot_dat$N_filtered, decreasing = T), 'batch']
plot_dat <- plot_dat[, c('batch', 'original', 'remaining')]

plot_dat <- melt(plot_dat)
colnames(plot_dat) <- c('batch', 'pre_or_post', 'cells')
ggplot(plot_dat, aes(x = cells, color = pre_or_post)) + geom_histogram() +
  ylab('N samples') + xlab('N cells')


#sample_order <- c("Su9_03_B_cells", "Su13_03_Innate", "Su19_3_Innate")
plot_dat <- mutate(plot_dat, batch=factor(plot_dat$batch, levels = sample_order))

p <- ggplot(plot_dat, aes(x = cells, y = batch)) + 
  geom_line() + theme(axis.text=element_text(size = 6)) + ggtitle('Cell filter N before/after') + 
  geom_point(size = 2, alpha = 0.3, aes(color = pre_or_post))
print(p)
dev.off()
