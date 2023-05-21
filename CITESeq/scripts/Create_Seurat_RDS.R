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
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R')

if(interactive()){
  project <- '2021614_21-002'
  covs_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/Sample_sheet.csv'
  # project <- '2022619_857.3b'
  # covs_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022619_857.3b/Sample_sheet.csv'
  runs_dir <- '/hpcdata/vrc/vrc1_data/douek_lab/Runs/'
  labels_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/Cell_data.csv')
  out_rds <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/All_data.RDS')
  out_pdf <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/All_data.pdf')
} else{
  args = commandArgs(trailingOnly=TRUE)
  
  runs_dir <- args[1]
  project <- args[2]
  covs_file <- args[3]
  labels_file <- args[4]
  out_rds <- args[5]
  out_pdf <- args[6]  
}

print(out_rds)
project_dir <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project)
cellranger_quant <- c('multi', 'count')

### Read covariates
dat <- read.csv(covs_file, stringsAsFactors = F, header = T,
                 check.names = F, 
                 colClasses = 'character',
                 sep = ',')
row.names(dat) <- dat$Sample_ID

### Add Read count info per sample from "multi_output/CRN/outs/per_sample_outs/CRN/metrics_summary.csv"
dat[, c('nRead_RNA', 'nRead_prot', 'nRead_VDJ', 'Sequencing_saturation')] <- NA

for(cr_id in unique(dat$CR_ID)){
  metrics_file <- paste0(project_dir, '/multi_output/', cr_id, '/outs/per_sample_outs/', cr_id, '/metrics_summary.csv')
  met <- read.table(metrics_file, sep = ',', header = T)
  nread <- met[which(met[, 'Library.Type'] == 'Gene Expression' & met[, 'Grouped.By'] == 'Physical library ID' &  met[, 'Metric.Name'] == 'Number of reads'), 'Metric.Value']
  nread <- as.numeric(gsub(',', '', nread))
  dat[which(dat$CR_ID == cr_id), 'nRead_RNA'] <- nread
  
  nread <- met[which(met[, 'Library.Type'] == 'Antibody Capture' & met[, 'Grouped.By'] == 'Physical library ID' &  met[, 'Metric.Name'] == 'Number of reads'), 'Metric.Value']
  nread <- as.numeric(gsub(',', '', nread))
  dat[which(dat$CR_ID == cr_id), 'nRead_prot'] <- nread
  
  nread <- met[which(met[, 'Library.Type'] == 'Antibody Capture' & met[, 'Grouped.By'] == 'Physical library ID' &  met[, 'Metric.Name'] == 'Number of reads'), 'Metric.Value']
  nread <- as.numeric(gsub(',', '', nread))
  dat[which(dat$CR_ID == cr_id), 'nRead_prot'] <- nread
  
  if('VDJ B' %in% unique(met[, 'Library.Type'])){
    nread <- met[which(met[, 'Library.Type'] == 'VDJ B' & met[, 'Grouped.By'] == 'Physical library ID' &  met[, 'Metric.Name'] == 'Number of reads'), 'Metric.Value']
    nread <- as.numeric(gsub(',', '', nread))
    dat[which(dat$CR_ID == cr_id), 'nRead_VDJ'] <- nread
  }
  
  ### Saturation
  sat <- met[which(met[, 'Library.Type'] == 'Gene Expression' & met[, 'Grouped.By'] == 'Physical library ID' &  met[, 'Metric.Name'] == 'Sequencing saturation'), 'Metric.Value']
  sat <- as.numeric(gsub('%', '', sat)) / 100
  dat[which(dat$CR_ID == cr_id), 'Sequencing_saturation'] <- sat
}

plotdat <- unique(dat[, c('Date_sort', 'nRead_RNA')])
p <- ggplot(plotdat, aes(x = Date_sort, y = nRead_RNA)) + 
  geom_dotplot(binaxis = 'y', stackdir = 'center') + ylim(c(0, (max(plotdat$nRead_RNA) * 1.1)))
p

plotdat <- unique(dat[, c('Date_sort', 'Sequencing_saturation')])
p <- ggplot(plotdat, aes(x = Date_sort, y = Sequencing_saturation)) + 
  geom_dotplot(binaxis = 'y', stackdir = 'center') + ylim(c(0, (max(plotdat$Sequencing_saturation) * 1.1)))
p

htos <- unique(dat$HTO_index_name)

labels <- read.table(labels_file, header = T, sep = ',')
row.names(labels) <- labels$cell_id

### If some of the sample IDs in labels are not in dat, select the cell_ids that are (in this project)
if(!all(labels[which(!labels$Assignment %in% c('Negative', 'Doublet')), 'Assignment'] %in% dat$Sample_ID)){
  print(paste0('Selecting cells from samples within the input sample sheet: ', toString(dat$Sample_ID)))
  project_cells <- labels[which(labels$Assignment %in% dat$Sample_ID), 'cell_id']    
} else{
  project_cells <- row.names(labels)
}

rlist <- as.list(rep(NA, length(unique(dat$CR_ID))))
names(rlist) <- unique(dat$CR_ID)
plist <- as.list(rep(NA, length(unique(dat$CR_ID))))
names(plist) <- unique(dat$CR_ID)
### For each cellranger count output, read and filter based on info in 'labels'
for(id in unique(dat$CR_ID)){
  print(id)
  ### Get desired directory of data
  data_dir <- my_read_10x(id, dat, runs_dir, project_dir, cellranger = cellranger_quant, filtered = 'raw')
  
  CR_dat <- Read10X(data.dir = data_dir)
  ### RNA data
  rna <- CR_dat[['Gene Expression']]
  ### Account for weird "multi_" prefacing cell_id if processed by cellranger multi
  colnames(rna) <- sapply(colnames(rna), function(cn) strsplit(cn, '_')[[1]][2])
  ### Remove the trailing "-1"
  colnames(rna) <- gsub('-1$', '', colnames(rna))
  ### Prepend CR ID
  colnames(rna) <- paste0(id, '_', colnames(rna))
  ### Protein data
  prot <- CR_dat[[names(CR_dat)[which(names(CR_dat) != 'Gene Expression')]]]
  ### Account for weird "multi_" prefacing cell_id if processed by cellranger multi
  colnames(prot) <- sapply(colnames(prot), function(cn) strsplit(cn, '_')[[1]][2])
  ### Remove the trailing "-1"
  colnames(prot) <- gsub('-1$', '', colnames(prot))
  ### Prepend CR id
  colnames(prot) <- paste0(id, '_', colnames(prot))
  row.names(prot) <- gsub('_totalseq', '', row.names(prot))
  rm(CR_dat)
  

  ### Keep cells by filtering on the dehashing summary column, Assignment_simple
  positive_cells <- labels[which(labels$Assignment_simple == 'Positive' & labels$cell_id %in% colnames(rna)), , drop = F]$cell_id
  if(!all(positive_cells %in% project_cells)){
    print(paste0('Removing ', sum(!(positive_cells %in% project_cells)), ' cells from samples not in this project.'))
    positive_cells <- positive_cells[which(positive_cells %in% project_cells)]
  }
  rna <- rna[, which(colnames(rna) %in% positive_cells)]
  prot <- prot[, which(colnames(prot) %in% positive_cells)]
  
  rlist[[id]] <- rna
  plist[[id]] <- prot
  rm(rna)
  rm(prot)
}
print('combining RNA')
### Combine RNA data. Conversion to dataframe crashes, so can't use rblindlist
rna <- rlist[[1]]
for(i in 2:length(rlist)){
  rna <- cbind(rna, rlist[[i]])
}
colnames(rna) <- unlist(lapply(rlist, function(x) colnames(x)))
print('combining protein')
### Combine prot data. 
markers <- unique(unlist(lapply(plist, function(prot) row.names(prot))))
addEmptyRows <- function(prot, markers){
  toAdd <- markers[which(!markers %in% row.names(prot))]
  myMatrix <- matrix(ncol = ncol(prot), nrow = length(toAdd))
  row.names(myMatrix) <- toAdd
  prot <- rbind(prot, myMatrix)
  return(prot)
}
### Add empty rows when fewer than maximum hashes
prot <- addEmptyRows(plist[[1]], markers)
for(i in 2:length(plist)){
  psub <- addEmptyRows(plist[[i]], markers)
  prot <- cbind(prot, psub)
}
colnames(prot) <- unlist(lapply(plist, function(x) colnames(x)))

print('editing labels')
### Subset labels
labels <- labels[which(labels$cell_id %in% colnames(rna)),]
colnames(labels)[which(colnames(labels) == 'Assignment')] <- 'Sample_ID'
### Add useful info from dat to labels (which will become Seurat metadata)
dat <- dat[, sapply(colnames(dat), function(x) length(unique(dat[, x])) > 1)]
dat <- dat[, !colnames(dat) %in% colnames(labels)]
labels[, colnames(dat)] <- NA
for(s in unique(labels$Sample_ID)){
  labels[which(labels$Sample_ID == s), colnames(dat)] <- dat[s, ]
}

print('making seurat')
### create Seurat object with cell-containing drops (min.cells is a gene filter, not a cell filter)
cseq <- Seurat::CreateSeuratObject(counts = rna, meta.data = labels, assay = "RNA", min.cells = 20, project = 'CITESeq')

print('adding prot')
#row.names(prot) <- gsub('-totalseq', '', row.names(prot))
### Remove htos
print(dim(prot))
prot <- prot[which(!(row.names(prot) %in% htos)),]
print(dim(prot))

cseq[["prot"]] <- Seurat::CreateAssayObject(data = prot)
print('saving')
saveRDS(cseq, file = out_rds)

pdf(out_pdf)
VlnPlot(cseq, features = "nFeature_RNA", combine = T, log = T, pt.size = 0)
VlnPlot(cseq, features = "nCount_RNA", combine = T, log = T, pt.size = 0)
dev.off()
