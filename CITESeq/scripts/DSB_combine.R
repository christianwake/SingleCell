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

source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

if(interactive()){
  # project <- '2021614_21-002'
  # qc_name <- '2023-05-30'
  # batch <- 'Sample_Name'
  # batches <- c('Su13_03_Innate', 'Su19_3_Innate', 'Su9_03_B_cells')

  # project <- '2023600_21-0012'
  # qc_name <- '2023-11-07'
  # batch <- 'CR_ID'
  # batches <- c('CR1-2', 'CR2-2')
  # isotype_controls <- 'C0090,C0091,C0092,C0095'
  
  project <- '2024605_Hillary_test'
  qc_name <- '2024-06-04'
  batch <- 'CR_ID'
  batches <- c('CR1')
  isotype_controls <- ''
  
  out_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', 
                     qc_name, '/DSB_normalized_data.RDS')
  pdf_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', 
                     qc_name, '/DSB_normalized_data.pdf')
  sdat_files <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, 
                       '/results/', qc_name, '/batches/', batches, '/Cell_filtered.RDS')
  
  mat_files <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project,
                      '/results/', qc_name, '/batches/', batches, '/DSB_normalized_data.RDS')
}else{
  args <- commandArgs(trailingOnly=TRUE)
  
  out_file <- args[1]
  pdf_file <- args[2]
  batch <- args[3]
  isotype_controls <- args[4]
  all_files <- args[5:length(args)]
  sdat_files <- all_files[1:(length(all_files)/2)]
  mat_files <- all_files[((length(all_files)/2)+1):length(all_files)]
}

isotype_controls <- strsplit(isotype_controls, ',')[[1]]
##### Standardize protein names
### space then parentheses <- parentheses
isotype_controls <- gsub(' \\(', '\\(', isotype_controls)
### spaces, underscores to '.'
isotype_controls <- make.names(isotype_controls, allow_ = 'F')

### RNA counts is integer values. RNA data is integer values (probably identical)
### Don't have @assays$prot@counts, has integer @assays$prot@data
sdats <- lapply(sdat_files, function(sdat_file) readRDS(sdat_file))

### NOTE - merge doesn't work on an assay if the counts slot is empty but the data slot is full. Then both will end up empty.
for(i in 1:length(sdats)){
  print('per batch')
  print(paste0('Mean data: ', mean(colSums(sdats[[i]]@assays$prot@data))))
  sdats[[i]]@assays$prot@counts <- sdats[[i]]@assays$prot@data
  ### Print mean nCount_prot
  print(paste0('Mean counts: ', mean(colSums(sdats[[i]]@assays$prot@counts))))
}

### merge.data = F will get rid of anything in 'data' and make them 'counts' instead. This is necessary if you have DSB normalized values in 'data'
### because it will try to normalize after the merge, which causes an error. Thus, Seurat objects normalized separately cannot be correctly merged.
if(length(sdats) == 1){
  sdat <- sdats[[1]]
} else{
  sdat <- merge(x = sdats[[1]], y = array(unlist(sdats[2:length(sdats)])), merge.data = F, project = 'SeuratProject')
}
DefaultAssay(sdat) <- 'prot'
print('Merged')
print(paste0('Mean counts: ', mean(colSums(sdats[[i]]@assays$prot@counts))))
print(paste0('Mean data: ', mean(colSums(sdats[[i]]@assays$prot@data))))
row.names(sdat@assays$prot@data) <- gsub('-totalseq', '', row.names(sdat@assays$prot))
row.names(sdat@assays$prot@counts) <- gsub('-totalseq', '', row.names(sdat@assays$prot))

### Remove isotype controls (they aren't in the DSB outputs)
to_keep <- row.names(sdat@assays$prot@data)[which(!row.names(sdat@assays$prot@data) %in% isotype_controls)]
sdat@assays$prot@data <- sdat@assays$prot@data[to_keep, ]
sdat@assays$prot@counts <- sdat@assays$prot@counts[to_keep, ]
#rm(sdats)

### Read DSB normalized counts files
mats <- lapply(mat_files, function(mat_file) readRDS(mat_file))
### Merge DSB counts
mat <- mats[[1]]
if(length(mats) > 1){
  for(i in 2:length(mats)){
    mat <- cbind(mat, mats[[i]])
  }
}

mat <- mat[row.names(sdat@assays$prot),]
print('matrix')
print(mean(colSums(mat)))
#rm(mats)
### Add DSB counts as 'data' to sdat
#sdat@assays$prot@counts <- sdat@assays$prot@data
sdat@assays$prot@data <- mat

print('Final?')
print(mean(colSums(sdat@assays$prot@counts)))
print(mean(colSums(sdat@assays$prot@data)))

print(range(sdat@assays$prot@counts))
print(range(sdat@assays$prot@data))
sdat$nFeature_prot <- colSums(sdat@assays$prot@counts > 0)
sdat$nCount_prot <- colSums(sdat@assays$prot@counts)
sdat$nFeature_dsb <- colSums(sdat@assays$prot@data > 0)
sdat$nCount_dsb <- colSums(sdat@assays$prot@data)
### Save sdat
saveRDS(sdat, file = out_file)

pdf(pdf_file)
VlnPlot(sdat, assay = 'prot', features = 'nCount_prot', group.by = batch, pt.size = 0) + theme(legend.position = "none")
VlnPlot(sdat, assay = 'prot', features = 'nCount_dsb', group.by = batch, pt.size = 0) + theme(legend.position = "none")
dev.off()

### create Seurat object with cell-containing drops (min.cells is a gene filter, not a cell filter)
# sdat <- Seurat::CreateSeuratObject(counts = rna[ , positive_cells], meta.data = metadata[positive_cells, , drop = F], assay = "RNA", min.cells = 20)
# 
# ### add DSB normalized "dsb_norm_prot" protein data to an assay called "CITE" created in step II 
# sdat[["prot"]] = Seurat::CreateAssayObject(data = dsb_norm_prot)
# ### Add the unnormalized data (but for only the cells that carried through by DSB) as 'counts'
# ### Seurat doesn't allow underscores in feature names
# row.names(prot) <- gsub('_', '-', row.names(prot)) 
# sdat@assays$prot@counts <- as.sparse(prot[row.names(sdat@assays$prot), colnames(sdat@assays$prot)])
# 
# ### Save seurat object as an RDS file
# if(!is.na(out_file)){
#   ### Get current date, to be used in the file name
#   #rds_file <- paste0(out_dir, panel, '_DSBnormalized_Seurat_', date, '.rds')
#   print(paste0('Saving file ', out_file))
#   saveRDS(sdat, file = out_file)
# }
