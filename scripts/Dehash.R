library('sys')
library('Seurat')
library('ggplot2')
library('dplyr')
library('viridis')
library('data.table')
library('patchwork')
library('PKI')
library('tinytex')
#library('dsb')
#library('tidyverse')
#library('pastecs')
library('reticulate')
#library('umap')
library('gridExtra')
library('cowplot')
#library('logspline')
library('readxl')
#library('WriteXLS')
library('Matrix')
library('pastecs')

source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sample_sheet_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/dehash_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R')

if(interactive()){
  # project <- '2021614_21-002'
  # runs_dir <- '/hpcdata/vrc/vrc1_data/douek_lab/Runs/'
  # #demux_method <- 'Custom'
  # demux_method <- 'MULTIseqDemux'
  # covs_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/Sample_sheet.csv')
  # out_labels <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/Dehash_calls_', demux_method, '.tsv')
  # thresh_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/Dehash_threshs_', demux_method, '.csv')
  # out_metrics <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/', project, '/cellranger_count_metrics.tsv')
  # batch <- 'Date_sort'
  
  # project <- '2023600_21-0012'
  # runs_dir <- '/hpcdata/vrc/vrc1_data/douek_lab/Runs/'
  # demux_method <- 'Custom'
  # covs_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/Sample_sheet.csv')
  # out_labels <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/Dehash_calls_', demux_method, '.tsv')
  # thresh_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/Dehash_threshs_', demux_method, '.csv')
  # out_metrics <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/', project, '/cellranger_count_metrics.tsv')
  # batch <- 'Date_sort'
  
  
  project <- '2024605_Hillary_test'
  runs_dir <- '/hpcdata/vrc/vrc1_data/douek_lab/Runs/'
  #demux_method <- 'Custom'
  demux_method <- 'Trough'
  covs_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/Sample_sheet.csv')
  out_labels <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/Dehash_calls_', demux_method, '.tsv')
  thresh_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/Dehash_threshs_', demux_method, '.csv')
  out_metrics <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/', project, '/cellranger_count_metrics.tsv')
  batch <- 'FC_ID'
  
  out_pdf1 <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/Dehash_feature_counts_', demux_method, '.pdf')
  out_pdf2 <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/Dehash_cell_counts_', demux_method, '.pdf')
  gtf_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/gtf.RDS')
}else{
  args <- commandArgs(trailingOnly=TRUE)
  
  covs_file <- args[1]
  gtf_file <- args[2]
  runs_dir <- args[3]
  project <- args[4]
  demux_method <- as.character(args[5])
  thresh_file <- args[6]
  out_labels <- args[7]
  out_pdf1 <- args[8]
  out_pdf2 <- args[9]
}

### Determine whether method is simply reading thresholds from a file and applying them.
demux_method0 <- demux_method
### If thresh_file exists, read it
if(file.exists(thresh_file)){
  print('Threshold file already exists, so will be read and used.')
  thresh_dat <- read.table(thresh_file, sep = ',', header = T, quote = '')
  ### If it seems that there was no header, re read it and assume header names
  if(!any(c('CR_ID', 'hash_id', 'Trough', "MULTIseqDemux") %in% colnames(thresh_dat))){
    thresh_dat <- read.table(thresh_file, sep = ',', header = F, quote = '')
    colnames(thresh_dat) <- c('CR_ID', 'hash_id', demux_method)
  } else{
    ### Allow for the input method to refer to values from the other methods in the file
    demux_methods <- c('Trough', 'MULTIseqDemux', 'Custom')
    demux_methods <- demux_methods[which(demux_methods %in% colnames(thresh_dat))]
    if(!demux_method %in% colnames(thresh_dat)){
      warning('Input demux method is not within the input threshold file columns')
    }
    ### Replace references to other columns with their value in that column
    thresh_dat[, demux_method] <- sapply(1:length(row.names(thresh_dat)), function(i)
      fetch_threshold(demux_methods, demux_method, thresh_dat, 
                      thresh_dat[i, 'CR_ID'], thresh_dat[i, 'hash_id']))
    
    thresh_dat <- thresh_dat[c('CR_ID', 'hash_id', demux_method)]
  }
  demux_method <- 'InputThresholdFile'
} else{
  thresh_dat <- NA
  print('No threshold csv file exists.')
  ### Should not exist if thresh_dat == NA but just in case it does, remove it
  if(!is.na(thresh_file)){
    if(file.exists(thresh_file)){
      file.remove(thresh_file, showWarnings = F)
    }
  }
}

### Laurens' first run of 21-0012 used an old version of CR which output a different path name
if(project == '2023105_21-0012'){
  CR_version = '1.1.1'
} else{
  CR_version = '7.1.1'
}
method_options <- c('MULTIseqDemux', 'Trough', 'Custom', 'InputThresholdFile')
if(!demux_method %in% method_options){
  warning(paste0('Input demultiplex method(s) are not among the options (', toString(method_options), ')'))
}
### Convert inputs to the exact case of the method_options
demux_method <- method_options[which(toupper(method_options) == toupper(demux_method))]

print(paste0('Using method ', demux_method))
project_dir <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project)
### Read covariates
dat <- read.csv(covs_file, stringsAsFactors = F, header = T,
                 check.names = F, 
                 colClasses = 'character',
                 sep = ',')
### If the sheet is a replicate sheet rather than sample sheet
# if('Rep_ID' %in% colnames(dat)){
#   dat$Sample_ID <- dat$Rep_ID
# }
row.names(dat) <- dat$Sample_ID
##### Standardize protein names
### space then parentheses <- parentheses
dat$HTO_index_name <- gsub(' \\(', '\\(', dat$HTO_index_name)
### spaces, underscores to '.'
dat$HTO_index_name<- make.names(dat$HTO_index_name, allow_ = 'F')

### Do direct matching from gtf
if(grepl('\\.gtf', gtf_file)){
  gtf <- read_gtf(gtf_file, feature_type = 'gene', atts_of_interest = c('gene_id', 'gene_name',  'gene_biotype'))
} else if(grepl('\\.RDS', gtf_file)){
  gtf <- readRDS(gtf_file)
} else{
  warning('No gtf file')
}
row.names(gtf) <- gtf$gene_id

label_list <- as.list(rep(NA, length(unique(dat$CR_ID))))
names(label_list) <- unique(dat$CR_ID)
### For each cellranger count output, read and demultiplex
for(id in unique(dat$CR_ID)){
  print(id)
  #pdf_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/', id, '_hash.pdf')
  pdf_file <- NA
  ### Dehash Wrapper
  label_list[[id]] <- dehash_CR_ID(dat, id, runs_dir, project_dir, demux_method = demux_method, 
                                   autoThresh = F, gtf = gtf, thresh_file, thresh_dat, pdf_file, 
                                   CR_version)
}

#downdir <- "/hpcdata/vrc/vrc1_data/douek_lab/Runs/211123_A00243_0139_BHLKCFDSX2/count_output/CR1/outs/downsampled_prot/"
labels <- rbindlist(label_list)
labels <- as.data.frame(labels)
print(table(labels$Assignment_CR, labels$Assignment_simple))

row.names(labels) <- unlist(lapply(label_list, function(x) row.names(x)))
table(labels$Assignment)
labels[, 'cell_id'] <- row.names(labels)

# if(project == '2021614_21-002'){
#   disconcordance <- which(labels$Assignment != 'Negative' & labels$Assignment_CR == 'Negative')
#   labels[disconcordance, 'Assignment'] <- 'Negative'
#   labels[disconcordance, 'Assignment_simple'] <- 'Negative'
# }
nums <- as.data.frame(table(labels$Assignment))
colnames(nums) <- c('Assignment', 'N_demultiplexed')
table(nums$Assignment)
write.table(labels, out_labels, sep = '\t', quote = F, row.names = F, col.names = T)

# ### log10(RNA count) x log10(protein count) density plot, split by CellRanger call
# p2 <- ggplot(labels, aes(x = rna_size, y = prot_size)) + labs(title = paste0(id, ' (by cellranger call)'), y='log10(protein count)', x = 'log10(RNA count)') +
#   geom_bin2d(bins = 150) +
#   scale_fill_viridis_c(option = "C") +
#   facet_wrap(~Assignment_CR)
pdf(out_pdf1)
freqtable <- table(labels$Assignment_simple, labels$Assignment_CR)
# ann_text <- data.frame(rna_size = 3, prot_size = 3,
#                        Assignment_simple = factor("Cells",levels = c("Background", "Cells")))
### log10(RNA count) x log10(protein count) by Dehash call, split by CellRanger call
p3 <- ggplot(labels, aes(x = rna_size, y = prot_size, col = Assignment_simple)) + 
  labs(title = ' All data (by cellranger call)', y='log10(protein count)', x = 'log10(RNA count)', color = 'Dehash call') +
  geom_point(size =0.5 ) +
  facet_wrap(~Assignment_CR) +
  annotation_custom(tableGrob(freqtable), xmin=1, xmax=4, ymin=0, ymax=2)
  #geom_table() ### In ggpmisc?
print(p3)

### Pngs are required if data is large
png_path <- paste0(dirname(out_pdf1), '/Dehash_feature_counts/')
dir.create(png_path)
png(paste0(png_path, 'All_CRs_', demux_method0, '.png'))
print(p3)
dev.off()
for(id in unique(dat$CR_ID)){
  sublab <- labels[which(labels$CR_ID == id), ]
  freqtable <- table(sublab$Assignment_simple, sublab$Assignment_CR)
  p3 <- ggplot(sublab, aes(x = rna_size, y = prot_size, col = Assignment_simple)) + 
    labs(title = paste0(id, ' (by cellranger call)'), y='log10(protein count)', x = 'log10(RNA count)', 
         color = paste0(demux_method0, '\ndehash call')) +
    geom_point(size =0.5 ) +
    facet_wrap(~Assignment_CR) +
    annotation_custom(tableGrob(freqtable), xmin=1, xmax=4, ymin=0, ymax=2)
  
  png(paste0(png_path, id, '_', demux_method0, '.png'))
  print(p3)
  dev.off()
  print(p3)
}

dev.off()
  
### For each CR_ID, get the number of cells according to Cell Ranger
CR_counts <- sapply(unique(dat$CR_ID), function(cr_id) table(labels[which(labels$CR_ID == cr_id), 'Assignment_CR'])[['Cells']])
### For each CR_ID, for each Sample_ID in that CR_ID, get the number of cells according to the dehashing
Dehash_counts <- lapply(unique(dat$CR_ID), function(cr_id) 
  table(labels[which(labels$CR_ID == cr_id), 'Assignment'])[dat[which(dat$CR_ID == cr_id), 'Sample_ID']]
)
### Sum to cells per CR_ID (according to dehashing method)
Dehash_counts <- sapply(Dehash_counts, function(x) sum(as.numeric(x)))
names(Dehash_counts) <- unique(dat$CR_ID)

### Plot 
plotdat <- data.frame(Dehash_counts = Dehash_counts, CR_counts= CR_counts)

pdf(out_pdf2)
p3 <- ggplot(plotdat, aes_string(x = 'CR_counts', y = 'Dehash_counts', label = 'row.names(plotdat)')) + 
  geom_text(size = 3) + geom_abline(slope = 1)
print(p3) 

if('Cell_Count' %in% colnames(dat)){
  Sort_counts <- lapply(unique(dat$CR_ID), function(cr_id) dat[which(dat$CR_ID == cr_id), 'Cell_Count'])
  Sort_counts <- sapply(Sort_counts, function(x) sum(as.numeric(x)))
  names(Sort_counts) <- unique(dat$CR_ID)
  plotdat$Sort_counts <- Sort_counts
  ###
  p4 <- ggplot(plotdat, aes_string(x = 'CR_counts', y = 'Sort_counts', label = 'row.names(plotdat)')) + 
    geom_text(size = 3) + geom_abline(slope = 1)
  print(p4) 
  p5 <- ggplot(plotdat, aes_string(x = 'Dehash_counts', y = 'Sort_counts', label = 'row.names(plotdat)')) + 
    geom_text(size = 3) + geom_abline(slope = 1)
  print(p5) 
  
  plotdat$CR_ID <- row.names(plotdat)
  pdat <- melt(plotdat)
  colnames(pdat) <- c('CR_ID', 'type', 'count')
  pdat$type_abbrev <- sapply(pdat$type, function(x) substr(x, 1, 1))
  pdat$log2_count <- log2(as.numeric(pdat$count))
    
  p6 <- ggplot(pdat, aes_string(x = 'CR_ID', y = 'count', label = 'type_abbrev', color = 'type')) +
    geom_text(size = 3)
  print(p6)    
  p7 <- ggplot(pdat, aes_string(x = 'CR_ID', y = 'log2_count', label = 'type_abbrev', color = 'type')) +
    geom_text(size = 3)
  print(p7)   
}

dev.off()
