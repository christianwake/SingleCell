
library('sys')
library('Seurat')
library('ggplot2')
library('dplyr')
library('viridis')
library('data.table')
library('patchwork')
library('PKI')
library('tinytex')
library('dsb')
#library('tidyverse')
#library('pastecs')
library('reticulate')
library('umap')
library('gridExtra')
library('cowplot')
library('logspline')
library('readxl')
library('WriteXLS')
library('Matrix')
library('pastecs')

# source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
# source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/sample_sheet_functions.R')
# source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
# source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R')

source('/Volumes/VRC1_DATA/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
source('/Volumes/VRC1_DATA/douek_lab/wakecg/sample_sheet_functions.R')
source('/Volumes/VRC1_DATA/douek_lab/snakemakes/sc_functions.R')
source('/Volumes/VRC1_DATA/douek_lab/snakemakes/Utility_functions.R')

project <- '2022619_857.3b'
# runs_dir <- '/hpcdata/vrc/vrc1_data/douek_lab/Runs/'
# covs_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022619_857.3b/Sample_sheet.csv'
# ### Auto thresh off
# multiseq_off_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022619_857.3b/data/Cell_data_old.csv'
# ### Auto thresh on
# multiseq_on_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022619_857.3b/data/Cell_data_new.csv'
# ### Custom method
# custom_labels_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022619_857.3b/hash_calls.tsv'


runs_dir <- '/Volumes/VRC1_DATA/douek_lab/Runs/'
covs_file <- '/Volumes/VRC1_DATA/douek_lab/projects/RNASeq/2022619_857.3b/Sample_sheet.csv'
### Auto thresh off
multiseq_off_file <- '/Volumes/VRC1_DATA/douek_lab/projects/RNASeq/2022619_857.3b/data/Cell_data_old.csv'
### Auto thresh on
multiseq_on_file <- '/Volumes/VRC1_DATA/douek_lab/projects/RNASeq/2022619_857.3b/data/Cell_data_new.csv'
### Custom method
custom_labels_file <- '/Volumes/VRC1_DATA/douek_lab/projects/RNASeq/2022619_857.3b/hash_calls.tsv'

# out_metrics <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022619_857.3b/BHLN7YDSX3/cellranger_count_metrics.tsv'
# quant_func <- 'count'
# batch <- 'Date_Sort'

project_dir <- paste0('/Volumes/VRC1_DATA/douek_lab/projects/RNASeq/', project)
# project_dir <- paste0('/Volumes/VRC1_DATA/douek_lab/projects/RNASeq/', project)
### Read covariates
dat <- read.csv(covs_file, stringsAsFactors = F, header = T,
                check.names = F, 
                colClasses = 'character',
                sep = ',')
row.names(dat) <- dat$Sample_ID

labels_mson <- read.csv(multiseq_on_file)
labels_msoff <- read.csv(multiseq_off_file)
row.names(labels_mson) <- labels_mson$cell_id
row.names(labels_msoff) <- labels_msoff$cell_id
colnames(labels_mson) <- c('Assignment_mson', colnames(labels_mson)[2:10])
labels_msoff <- labels_msoff[, c('Assignment'), drop = F]
colnames(labels_msoff) <- c("Assignment_msoff")
labels_ms <- cbind(labels_msoff, labels_mson)
labels_ms$Assignment_simple_msoff <- ifelse(labels_ms$Assignment_msoff %in% c('Doublet', 'Negative'), labels_ms$Assignment_msoff, 'Positive')
labels_ms$Assignment_simple_mson <- ifelse(labels_ms$Assignment_mson %in% c('Doublet', 'Negative'), labels_ms$Assignment_mson, 'Positive')
labels_cust$Assignment_simple_cust <- ifelse(labels_cust$Assignment %in% c('Doublet', 'Negative'),labels_cust$Assignment, 'Positive')

cust_list <- lapply(unique(dat$CR_ID), function(cr) {
       # custom_labels_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022619_857.3b/count_output/', cr, '/outs/all_hash_calls.tsv')
       custom_labels_file <- paste0('/Volumes/VRC1_DATA/douek_lab/projects/RNASeq/2022619_857.3b/count_output/', cr, '/outs/all_hash_calls.tsv')
       read.table(custom_labels_file, sep = '\t', header = T)
       })
names(cust_list) <- unique(dat$CR_ID)
labels_cust <- rbindlist(cust_list, idcol = T)
colnames(labels_cust) <- c('CR_ID', 'cell_id', 'assignment')


#### Stacked bar plots of C02XY per CR ID
a <-lapply(unique(dat$CR_ID), function(id) dehash_plot_dat(dat, id, runs_dir, project_dir))
names(a) <- unique(dat$CR_ID)
pdat <- rbindlist(lapply(a, function(x) { 
  pdat <- as.data.frame(x)
  pdat$hash <- row.names(pdat)
  pdat
}), idcol = T)
colnames(pdat) <- c('CR_ID', 'N_UMI', 'Hash')

pdat <- mutate(pdat, CR_ID = factor(CR_ID, levels=c(paste0('CR', 1:16))))
p0 <- ggplot(as.data.frame(pdat), aes(x = CR_ID, y = N_UMI, fill = Hash)) + geom_col(width = 0.3)
p1 <- ggplot(as.data.frame(pdat), aes(x = CR_ID, y = N_UMI, fill = Hash)) + geom_col(width = 0.3, position = 'fill')


### Fraction of prot-called Positives that are CR-called Cells, per CR_ID
CRID_disp <- sapply(unique(labels$CR_ID), function(cr) {
  sub <- labels[which(labels$CR_ID == cr & labels$Assignment_simple == 'Positive'), ]
  length(which(sub$Assignment_CR == 'Cells')) / dim(sub)[1]
})
names(CRID_disp) <- gsub('.Cells', '', names(CRID_disp))

SID_disp <- lapply(unique(labels$CR_ID), function(cr) {
  sub <- labels[which(labels$CR_ID == cr & labels$Assignment_simple == 'Positive'), ]
  sapply(unique(sub$Assignment), function(id) {
    sub2 <- sub[which(sub$Assignment == id), ]
    out <- length(which(sub2$Assignment_CR == 'Cells')) / dim(sub2)[1]
    #names(out) <- gsub('.Cells', '', names(out))
    out
  })
})
names(SID_disp) <- unique(labels$CR_ID)


p2 <- ggplot(labels, aes(x = rna_size, y = prot_size)) + geom_bin2d(bins = 150)+ labs(title = id, y='log10(protein count)', x = 'log10(RNA count)')

p3 <- ggplot(labels, aes(x = rna_size, y = prot_size)) + labs(title = paste0(id, ' (by cellranger call)'), y='log10(protein count)', x = 'log10(RNA count)') +
  geom_bin2d(bins = 150) +
  scale_fill_viridis_c(option = "C") +
  facet_wrap(~Assignment_CR)

p4 <- ggplot(labels, aes(x = rna_size, y = prot_size, col = Assignment_simple)) + 
  labs(title = paste0(id, ' (by cellranger call)'), y='log10(protein count)', x = 'log10(RNA count)', color = 'Dehash call') +
  geom_point(size =0.5 ) +
  #scale_fill_viridis_c(option = "C") +
  facet_wrap(~Assignment_CR)


print(p2)
print(p3)
print(p4)
