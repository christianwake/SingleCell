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
labels_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/data/Demultiplex_labels.csv'
out_rds <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/Biolegend/DSB_normalized_data.RDS'

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

plist <- as.list(rep(NA, length(unique(dat$CR_ID))))
names(plist) <- unique(dat$CR_ID)
### For each cellranger count output, read and demultiplex
for(id in unique(dat$CR_ID)){
  data_dir <- paste0(runs_dir,  id, runs_subdir)
  CR_dat <- Read10X(data.dir = data_dir)
  
  prot <- CR_dat[[names(CR_dat)[which(names(CR_dat) != 'Gene Expression')]]]
  colnames(prot) <- paste0(colnames(prot), '_', id)
  row.names(prot) <- gsub('_totalseq', '', row.names(prot))
  rm(CR_dat)
  
  plist[[id]] <- prot
  rm(prot)
}

### Combine protein data
prot <- plist[[1]]
for(i in 2:length(plist)){
  prot <- cbind(prot, plist[[i]])
}
colnames(prot) <- unlist(lapply(plist, function(x) colnames(x)))

