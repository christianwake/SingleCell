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

if(interactive()){
  # project <- '2021614_21-002'
  # covs_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/All_covariates.csv'
  
  project <- '2022619_857.3b'
  covs_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022619_857.3b/Sample_sheet.csv'
  
  runs_dir <- '/hpcdata/vrc/vrc1_data/douek_lab/Runs/'
  labels_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/Cell_data.csv')
  gtf_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/gtf.RDS')
  out_rds <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/All_data.RDS')
} else{
  args = commandArgs(trailingOnly=TRUE)
  runs_dir <- args[1]
  project <- args[2]
  covs_file <- args[3]
  labels_file <- args[4]
  gtf_file <- args[5]
  out_rds <- args[6]  
}

project_dir <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project)
### Read covariates
dat <- read.csv(covs_file, stringsAsFactors = F, header = T,
                 check.names = F, 
                 colClasses = 'character',
                 sep = ',')
row.names(dat) <- dat$Sample_ID

### Hash info
labels <- read.table(labels_file, header = T, sep = ',')
row.names(labels) <- labels$cell_id

### Read gtf
if(grepl('\\.gtf', gtf_file)){
  gtf <- read_gtf(gtf_file, feature_type = 'gene', atts_of_interest = c('gene_id', 'gene_name',  'gene_biotype'))
} else if(grepl('\\.RDS', gtf_file)){
  gtf <- readRDS(gtf_file)
}
row.names(gtf) <- gtf$gene_id

rlist <- as.list(rep(NA, length(unique(dat$CR_ID))))
names(rlist) <- unique(dat$CR_ID)
### For each cellranger count output, read and filter based on info in 'labels'
for(id in unique(dat$CR_ID)){
  print(id)
  ### Get desired directory of data
  data_dir <- my_read_10x(id, dat, runs_dir, project_dir, cellranger = 'count', filtered = 'raw')
  
  CR_dat <- Read10X(data.dir = data_dir)
  ### RNA data
  rna <- CR_dat[['Gene Expression']]
  colnames(rna) <- paste0(colnames(rna), '_', id)
  
  rm(CR_dat)
  ### Keep cells by filtering on the dehashing summary column, Assignment_simple
  positive_cells <- labels[which(labels$Assignment_simple == 'Positive' & labels$cell_id %in% colnames(rna)), , drop = F]$cell_id
  rna <- rna[, which(colnames(rna) %in% positive_cells)]
  ### Check ID type. 
  #gtf_cols <- c('gene_name', 'gene_id')
  #id_type <- get_id_name(row.names(rna), gtf)
  #a <- sapply(row.names(rna), function(gene) ifelse(gene %in% gtf[, 'gene_id'], gene, gtf[which(gtf[, 'gene_name'] == gene), 'gene_id']))

  rlist[[id]] <- rna
  rm(rna)
}
print('combining RNA')
### Combine RNA data. Conversion to dataframe crashes, so can't use rblindlist
rna <- rlist[[1]]
for(i in 2:length(rlist)){
  rna <- cbind(rna, rlist[[i]])
}
colnames(rna) <- unlist(lapply(rlist, function(x) colnames(x)))

print('editing labels')
### Subset labels
labels <- labels[which(labels$cell_id %in% colnames(rna)),]
colnames(labels)[which(colnames(labels) == 'Assignment')] <- 'Sample_ID'
##### Add useful info from dat (sample file) to labels (which will become Seurat labelsdata)
### Remove non-unique columns
dat <- dat[, sapply(colnames(dat), function(x) length(unique(dat[, x])) > 1)]
### Remove columns that are already in labels
dat <- dat[, !colnames(dat) %in% colnames(labels)]
### Insantiate new columns in labels
labels[, colnames(dat)] <- NA
### Loop over Sample_IDs
for(s in unique(labels$Sample_ID)){
  labels[which(labels$Sample_ID == s), colnames(dat)] <- dat[s, ]
}

### create Seurat object with cell-containing drops (min.cells is a gene filter, not a cell filter)
sdat <- Seurat::CreateSeuratObject(counts = rna, meta.data = labels, assay = "RNA", min.cells = 20)
saveRDS(sdat, file = out_rds)

VlnPlot(sdat, 'nCount_RNA', split.by = 'Assignment_CR' )
