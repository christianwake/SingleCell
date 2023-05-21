library('sys')
library('Seurat')
library('SeuratWrappers')
library('harmony')
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
library('umap')
library('gridExtra')
library('cowplot')
library('logspline')
library('readxl')
library('WriteXLS')
library('Matrix')

source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/sample_sheet_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R')

dehash_CR_ID <- function(dat, id, data_dir){
  hashtags <- as.data.frame(dat[which(dat$CR_ID == id),])[,'HTO_index_name']
  specimens <- as.data.frame(dat[which(dat$CR_ID == id),])[,'Sample_ID']
  names(specimens) <- hashtags
  specimens <- c(specimens, 'Negative', 'Doublet')
  names(specimens) <- c(hashtags, 'Negative', 'Doublet')
  ### Read filtered data
  CR_filtered <- Read10X(data.dir = gsub('raw', 'filtered', data_dir))
  stained_cells <- colnames(CR_filtered$`Gene Expression`)
  
  rm(CR_filtered)
  CR_raw <- Read10X(data.dir = data_dir)
  
  ### Get some RNA metadata
  rna <- CR_raw$`Gene Expression`
  # create metadata of droplet QC stats used in standard scRNAseq processing
  rna_size <- log10(Matrix::colSums(rna))
  ngene <- Matrix::colSums(rna > 0)
  mtgene <- grep(pattern = "^MT-", rownames(rna), value = TRUE)
  propmt <- Matrix::colSums(rna[mtgene, ]) / Matrix::colSums(rna)
  rm(rna)
  
  prot <- CR_raw[[names(CR_raw)[which(names(CR_raw) != 'Gene Expression')]]]
  rm(CR_raw)
  #colnames(prot) <- paste0(colnames(prot), '_', id)
  row.names(prot) <- gsub('_totalseq', '', row.names(prot))
  prot_size <- log10(Matrix::colSums(prot))
  
  ### hash
  labels <- solve_hashes_main(prot, hashtags, out_dir, demux_method)
  rm(prot)
  ### If it is not within the standards, make it Negative
  labels[which(!(labels$Assignment %in% c('Doublet', 'Negative', hashtags))), 'Assignment'] <- 'Negative'
  ### Convert from index ID to specimen ID
  labels[, 'Assignment'] <- specimens[labels[, 'Assignment']]
  table(labels$Assignment)
  
  labels[, c('propmt', 'rna_size', 'ngene', 'prot_size')] <- as.data.frame(cbind(propmt, rna_size, ngene, prot_size))
  labels$bc <- rownames(labels)
  labels$Assignment_CR = ifelse(test = row.names(labels) %in% stained_cells, yes = 'Cells', no = 'Background')
  table(labels$Assignment_CR)
  table(labels$Assignment_CR, labels$Assignment)
  # filter barcodes to only include those with data for both assays 
  labels <- labels %>% dplyr::filter(rna_size > 0 & prot_size > 0 )
  table(labels$Assignment)
  ### Add batch info
  # b <- names(table(dat[which(dat$CR_ID == id), batch]))
  # if(length(b) > 1){ 
  #   warning('More than one batch for a cellranger sample.')
  # }
  # labels[, batch] <- b
  labels[, 'CR_ID'] <- id
  labels$Assignment_simple <- ifelse(!labels$Assignment %in% c("Doublet", "Negative"), 'Positive', labels$Assignment)
  row.names(labels) <- paste0(row.names(labels), '_', id)
  
  p <- ggplot(labels, aes(x = log10(ngene), y = prot_size)) +
    theme_bw() + 
    geom_bin2d(bins = 300) + 
    scale_fill_viridis_c(option = "C") + 
    facet_wrap(~Assignment_simple)
  print(p)
  return(labels)
}

args <- commandArgs(trailingOnly=TRUE)

covs_file <- args[1]
runs_dir <- args[2]
out_labels <- args[3]

# runs_dir <- '/hpcdata/vrc/vrc1_data/douek_lab/Runs/210511_A00243_0102_AHC3YJDSX2/wakecg/'
# covs_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021617_mis-c/All_covariates.csv'
# out_rds <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021617_mis-c/seurat_object.RDS'
# out_labels <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021617_mis-c/labels.csv'
# out_csv <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021617_mis-c/RNA_counts.csv'

# runs_dir <- '/hpcdata/vrc/vrc1_data/douek_lab/Runs/'
# covs_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/All_covariates.csv'
# out_labels <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/data/labels_allbackground.csv'
# out_metrics <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/cellranger_count_metrics.tsv'
# quant_func <- 'count'
# batch <- 'Date_Sort'

### Read covariates
dat <- read.csv(covs_file, stringsAsFactors = F, header = T,
                check.names = F, 
                colClasses = 'character',
                sep = ',')
row.names(dat) <- dat$Sample_ID

### Combine cellranger count output metrics
# if(quant_func == 'multi'){
#   runs_subdir <- '/outs/multi/count/raw_feature_bc_matrix/'
#   metrics <- lapply(unique(dat$CR_ID), function(x) as.data.frame(read.table(paste0(runs_dir, x, '/outs/per_sample_outs/', x, '/metrics_summary.csv'), sep = ',')[2,]))
#   metrics <- rbindlist(metrics)
#   colnames(metrics) <- gsub(' ', '_', unlist(read.table(paste0(runs_dir, dat$CR_ID[1], '/outs/per_sample_outs/', dat$CR_ID[1], '/metrics_summary.csv'), sep = ',')[1,]))
# } else if(quant_func == 'count'){
#   runs_subdir <- '/outs/raw_feature_bc_matrix/'
#   metrics <- lapply(unique(dat$CR_ID), function(x) as.data.frame(read.table(paste0(runs_dir, x, '/outs/metrics_summary.csv'), sep = ',')[2,]))
#   metrics <- rbindlist(metrics)
#   colnames(metrics) <- gsub(' ', '_', unlist(read.table(paste0(runs_dir, dat$CR_ID[1], '/outs/metrics_summary.csv'), sep = ',')[1,]))
# }
# metrics <- as.data.frame(t(metrics))
# colnames(metrics) <- unique(dat$CR_ID)
# write.table(metrics, out_metrics, sep = '\t', quote = F,  row.names = T, col.names = T)

demux_method <- 'seuratDemux'
label_list <- as.list(rep(NA, length(unique(dat$CR_ID))))
names(label_list) <- unique(dat$CR_ID)
### For each cellranger count output, read and demultiplex
for(id in unique(dat$CR_ID)){
  ### Find the raw feature barcode matrix file for this cell ranger sample
  flowcell <- unique(dat[which(dat$CR_ID == id), 'flowcell_full'])
  ### Use multi if it exists, count if not. If neither exist, print warning and skip
  multi_dir <- c(file.path(runs_dir, flowcell, 'multi_output', id, 'outs', 'multi', 'count','raw_feature_bc_matrix'),
                 file.path(runs_dir, flowcell, 'multi_output', id, 'outs', 'multi', 'count','filtered_feature_bc_matrix'))
  count_dir <- c(file.path(runs_dir, flowcell, 'count_output', id, 'outs', 'raw_feature_bc_matrix'),
                 file.path(runs_dir, flowcell, 'count_output', id, 'outs', 'filtered_feature_bc_matrix'))
  if(all(file.exists(multi_dir))){
    data_dir <- multi_dir[1]
  } else if(all(file.exists(count_dir))){
    data_dir <- count_dir[1]
  } else{
    warning(paste0(id, ' Does Not Exist! Skipping.'))
    next
  }
  
  label_list[[id]] <- dehash_CR_ID(dat, id, data_dir)
}

labels <- rbindlist(label_list)
labels <- as.data.frame(labels)

row.names(labels) <- unlist(lapply(label_list, function(x) row.names(x)))
table(labels$Assignment)
labels[, 'cell_id'] <- row.names(labels)
nums <- as.data.frame(table(labels$Assignment))
colnames(nums) <- c('Assignment', 'N_demultiplexed')

write.table(labels, out_labels, sep = ',', quote = F, row.names = F, col.names = T)

