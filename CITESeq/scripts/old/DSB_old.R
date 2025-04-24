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

source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

if(interactive()){
  qc_name <- '2023-Jan'
  batch_value <- '2021-11-10'
  batch_name <- 'Date_Sort'
  sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/', qc_name, '/', batch_value, '/RNA_cell_filtered.RDS')
  br_rds <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/', qc_name, '/Background.RDS')
  covs_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/All_covariates.csv'
  labels_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/data/Cell_data.csv'
  hto_string <- 'C0251,C0252,C0253,C0254,C0255,C0256,C0257,C0258,C0259,C0260'
  out_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/', qc_name, '/', batch_value, '/DSB_normalized_data.RDS')
}else{
  args <- commandArgs(trailingOnly=TRUE)
  
  sdat_file <- args[1]
  br_rds <- args[2]
  covs_file <- args[3]
  labels_file <- args[4]
  batch_value <- args[5]
  batch_name <- args[6]
  hto_string <- args[7]
  ### These next two will be the same if there was only 1 panel in the config file (string, rather than file name input)
  # panel_file <- args[5]
  # panel <- args[6]
  out_file <- args[8]
  
}

### This file uses the "count" slot of the prot assay for DSB, and prints DSB results to a separate file (not Seurat)
panel_file <- ''
panel <- ''
### Read covariates
covs <- read.csv(covs_file, stringsAsFactors = F, header = T,
                check.names = F, 
                colClasses = 'character',
                sep = ',')
row.names(covs) <- covs$Sample_ID

sdat <- readRDS(sdat_file)
background <- readRDS(br_rds)
#row.names(sdat@assays$prot@data) <- gsub('-totalseq', '', row.names(sdat@assays$prot))
#row.names(background) <- gsub('_totalseq', '', row.names(background))
DefaultAssay(sdat) <- 'RNA'

print(all(row.names(background) == row.names(sdat@assays$prot@data)))

metadata <- read.table(labels_file, header = T, sep = ',')
table(metadata$Assignment_CR, metadata$Assignment_simple)
metadata <- metadata[which(metadata$Assignment != 'Doublet'),]
table(metadata$Assignment_CR, metadata$Assignment_simple)

### Subset to the input batch
### Subset sdat
keep_cells <- colnames(sdat)[which(sdat@meta.data[, batch_name] == batch_value)]
sdat <- subset(sdat, cells = keep_cells)
### Subset background. This assumes that the second _ section of the colnames is the cellranger ID
covs <- covs[which(covs[, batch_name] == batch_value),]
background <- background[, colnames(background)[sapply(colnames(background), function(x) strsplit(x, '_')[[1]][2]) %in% unique(covs$CR_ID)]]

md <- sdat@meta.data

prot <- sdat@assays$prot@data
print(mean(colSums(prot)))
htos <- get_hashtag_names(hto_string, row.names(prot))

if(file.exists(panel_file)){
  panel_info <- read.table(panel_file, sep = ',', header = T, check.names = F)
  ### If we have input panel_file and no hto_groups, lets assume that HTOs being with 'HTO_', and create the hto_groups object
  hto_groups <- sapply(htos, function(y) colnames(panel_info)[which(panel_info[y,] == 1)])
  ### Remove panel prefixes from protein names in 'prot' and 'panel_info'
  panel_info$name <- sapply(row.names(panel_info), function(rn) gsub(paste0('^',panel_info[rn, 'panel_prefix'],'_'), '', rn))
  row.names(prot) <- panel_info[row.names(prot), 'name']
  row.names(panel_info) <- panel_info$name
  # Add panel column to md
  md$hto_groups <- hto_groups[gsub('-', '_', md$Assignment)]
} else{
  hto_groups <- rep(panel, length(htos))
  names(hto_groups) <- htos
  panel_info <- NA
  md$hto_groups <- panel
}

# sdat@assays$prot@counts <- sdat@assays$prot@data
# ### Creates a Seurat object and saves to RDS, out_file
# sdat@assays$prot@data <- DSB_once(prot, md, hto_groups, negative_mtx_rawprot = background, panel_info, panel)
### Save  object
# saveRDS(sdat, file = out_file)

### Save as data frame. When read should be added to Seurat assay "data" slot
dsb <- DSB_once(prot, md, hto_groups, negative_mtx_rawprot = background, panel_info, panel)
print(mean(colSums(dsb)))
saveRDS(dsb, file = out_file)

### create Seurat object with cell-containing drops (min.cells is a gene filter, not a cell filter)
# cseq <- Seurat::CreateSeuratObject(counts = rna[ , positive_cells], meta.data = metadata[positive_cells, , drop = F], assay = "RNA", min.cells = 20)
# 
# ### add DSB normalized "dsb_norm_prot" protein data to an assay called "CITE" created in step II 
# cseq[["prot"]] = Seurat::CreateAssayObject(data = dsb_norm_prot)
# ### Add the unnormalized data (but for only the cells that carried through by DSB) as 'counts'
# ### Seurat doesn't allow underscores in feature names
# row.names(prot) <- gsub('_', '-', row.names(prot)) 
# cseq@assays$prot@counts <- as.sparse(prot[row.names(cseq@assays$prot), colnames(cseq@assays$prot)])
# 
# ### Save seurat object as an RDS file
# if(!is.na(out_file)){
#   ### Get current date, to be used in the file name
#   #rds_file <- paste0(out_dir, panel, '_DSBnormalized_Seurat_', date, '.rds')
#   print(paste0('Saving file ', out_file))
#   saveRDS(cseq, file = out_file)
# }
