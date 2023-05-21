
library('umap')
library('gridExtra')
library('cowplot')
library('logspline')
library('readxl')
library('WriteXLS')
library('Matrix')
library('Seurat')
library('data.table')
#source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/sample_sheet_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R')

args = commandArgs(trailingOnly=TRUE)
runs_dir <- args[1]
covs_file <- args[2]
gtf_file <- args[3]
out_rds <- args[4]

# barcodes <- args[1]
# features <- args[2]
# matrix <- args[3]
# labels_file <- args[4]
# RDS_file <- args[5]
# hto_string <- args[6]

### The Read10X R function uses input directory path, but for snakemake it is better to have files directly as the input/output.
# if(dirname(barcodes) == dirname(features) & dirname(features) == dirname(matrix)){
#   data_dir <- dirname(barcodes)
# } else{
#   print('Input data files are not in the same directory. This will cause an error.')
# }
# runs_dir <- '/hpcdata/vrc/vrc1_data/douek_lab/Runs/CIENI-MX/211015_NB501187_0098_AHFCK7BGXG/count_output/'
# runs_dir <- '/hpcdata/vrc/vrc1_data/douek_lab/Runs/'
# covs_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021618_galt/C7_10x_15oct21.csv'
# covs_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021618_galt/All_covariates.csv'
# out_rds <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021618_galt/data/Seurat_raw.RDS'

dat <- read.table(covs_file, sep = ',', header = T)
row.names(dat) <- dat$Sample_ID

### instatiate lists to hold dataframes for RNA data and metadata
rlist <- as.list(rep(NA, length(unique(dat$Sample_ID))))
names(rlist) <- unique(dat$Sample_ID)
mlist <- as.list(rep(NA, length(unique(dat$Sample_ID))))
names(mlist) <- unique(dat$Sample_ID)
### For each cellranger count output, read
for(id in unique(dat$Sample_ID)){
  print(id)
  ### Find the raw feature barcode matrix file for this cell ranger sample
  flowcell <- unique(dat[which(dat$Sample_ID == id), 'flowcell_full'])
  ### Use multi if it exists, count if not. If neither exist, print warning and skip
  multi_dir <- file.path(runs_dir, flowcell, 'multi_output', id, 'outs', 'multi', 'count','filtered_feature_bc_matrix')
  count_dir <- file.path(runs_dir, flowcell, 'count_output', id, 'outs', 'filtered_feature_bc_matrix')
  if(all(file.exists(multi_dir))){
    data_dir <- multi_dir
  } else if(all(file.exists(count_dir))){
    data_dir <- count_dir
  } else{
    warning(paste0(id, ' Does Not Exist! Skipping.'))
    next
  }
  
  
  rna <- Read10X(data.dir = data_dir)
  colnames(rna) <- paste0(colnames(rna), '_', id)
  rlist[[id]] <- rna
  
  ### Create data frame for meta data in 'dat', per cell
  meta <- as.data.frame(matrix(ncol = length(colnames(dat)), nrow = length(colnames(rna))))
  row.names(meta) <- colnames(rna)
  colnames(meta) <- colnames(dat)
  for(cov in colnames(dat)){
    meta[, cov] <- dat[id, cov]
  }
  mlist[[id]] <- meta
}

### Combine RNA data. Conversion to dataframe crashes, so can't use rblindlist
rna <- rlist[[1]]
for(i in 2:length(rlist)){
  rna <- cbind(rna, rlist[[i]])
}
#colnames(rna) <- unlist(lapply(rlist, function(x) colnames(x)))
### Combine meta info
meta <- as.data.frame(rbindlist(lapply(mlist, function(x) as.data.frame(x))))

### create Seurat object with cell-containing drops (min.cells is a gene filter, not a cell filter)
sdat <- Seurat::CreateSeuratObject(counts = rna, meta.data = meta, assay = "RNA", min.cells = 20)

saveRDS(sdat, file = out_rds)
