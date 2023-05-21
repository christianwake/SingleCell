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

source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/sample_sheet_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/Utility_functions.R')

args = commandArgs(trailingOnly=TRUE)

runs_dir <- arg[1]
covs_file <- args[2]
labels_file <- args[3]
hto_string <- args[4]
out_rds <- args[5]

runs_dir <- '/hpcdata/vrc/vrc1_data/douek_lab/Runs/211123_A00243_0139_BHLKCFDSX2/multi_output/'
covs_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/All_covariates.csv'
labels_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/data/labels.csv'

out_rds <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/Biolegend/Cells.RDS'
quant_func <- 'count'
### Read covariates
dat <- read.csv(covs_file, stringsAsFactors = F, header = T,
                 check.names = F, 
                 colClasses = 'character',
                 sep = ',')
row.names(dat) <- dat$Sample_ID

labels <- read.table(labels_file, header = T, sep = ',')

### 
if(quant_func == 'multi'){
  runs_subdir <- '/outs/multi/count/raw_feature_bc_matrix/'
} else if(quant_func == 'count'){
  runs_subdir <- '/outs/raw_feature_bc_matrix/'
}

rlist <- as.list(rep(NA, length(unique(dat$CR_ID))))
names(rlist) <- unique(dat$CR_ID)
### For each cellranger count output, read and demultiplex
for(id in unique(dat$CR_ID)){
  #data_dir <- paste0(runs_dir, id, '/outs/filtered_feature_bc_matrix/')
  data_dir <- paste0(runs_dir, id, runs_subdir)
  CR_dat <- Read10X(data.dir = data_dir)
  rna <- CR_dat[['Gene Expression']]
  colnames(rna) <- paste0(colnames(rna), '_', id)
  rm(CR_dat)
  sub_labels <- labels[which(labels$cell_id %in% colnames(rna)),]
  ### Remove Negatives and Doublets (don't do this if you will use DSB for CITESeq)
  specs <- unlist(dat[which(dat$Sample_ID != "20A-782-32"), 'Sample_ID'])
  positive_cells <- sub_labels[which(sub_labels$Assignment %in% c('Positive', specs)), , drop = F]$cell_id
  rna <- rna[, which(colnames(rna) %in% positive_cells)]
  
  rlist[[id]] <- rna
  rm(rna)
}

### Combine RNA data. Conversion to dataframe crashes, so can't use rblindlist
rna <- rlist[[1]]
for(i in 2:length(rlist)){
  rna <- cbind(rna, rlist[[i]])
}
colnames(rna) <- unlist(lapply(rlist, function(x) colnames(x)))

### create Seurat object with cell-containing drops (min.cells is a gene filter, not a cell filter)
# cseq <- Seurat::CreateSeuratObject(counts = rna, meta.data = labels, assay = "RNA", min.cells = 20)
# all(colnames(cseq) == row.names(labels))

saveRDS(rna, file = out_rds)

