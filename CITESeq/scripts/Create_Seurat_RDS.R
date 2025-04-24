library('sys')
library('Seurat')
library('ggplot2')
library('dplyr')
library('viridis')
library('data.table')
library('patchwork')
#library('PKI')
library('tinytex')
library('dsb')
#library('tidyverse')
#library('pastecs')
library('reticulate')
#library('umap')
library('gridExtra')
library('cowplot')
#library('logspline')
library('readxl')
library('WriteXLS')
library('Matrix')

source('/data/vrc_his/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/sample_sheet_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/sc_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/Utility_functions.R')

if(interactive()){
  # project <- '2021614_21-002'
  # qc_name <- 'B_cells_2024-04-05'
  # cell_type <- 'B'
  # demux_method <- 'Custom'
  # covs_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, 
  #                     '/Sample_sheet_B_cells.csv')
  
  project <-'2024605_Hillary_test'
  qc_name <- '2024-06-04'
  cell_type <- 'B'
  demux_method <- 'Trough'
  covs_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, 
                      '/Sample_sheet.csv')
  
  runs_dir <- '/data/vrc_his/douek_lab/Runs/'
  labels_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, 
                        '/data/Dehash_CRint_calls_', demux_method, '.tsv')
  labels_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, 
                        '/results/', qc_name, '/Dehash/Dehash_CRint_calls_', demux_method, '.tsv')
  
  out_rds <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/All_data.RDS')
  out_pdf <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/All_data.pdf')
  out_pdf2 <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Filter_info.pdf')
  } else{
  args = commandArgs(trailingOnly=TRUE)
  
  runs_dir <- args[1]
  project <- args[2]
  cell_type <- args[3]
  covs_file <- args[4]
  labels_file <- args[5]
  out_rds <- args[6]
  out_pdf <- args[7]  
  out_pdf2 <- args[8] 
}

#print(out_rds)
project_dir <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project)
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
  
  ncells <- met[which(met[, 'Library.Type'] == 'Gene Expression' & met[, 'Grouped.By'] == 'Physical library ID' &  met[, 'Metric.Name'] == 'Estimated number of cells'), 'Metric.Value']
  ncells <- as.numeric(gsub(',', '', ncells))
  dat[which(dat$CR_ID == cr_id), 'nCells_RNA'] <- ncells
  
  ncells <- met[which(met[, 'Library.Type'] == 'Antibody Capture' & met[, 'Grouped.By'] == 'Physical library ID' &  met[, 'Metric.Name'] == 'Estimated number of cells'), 'Metric.Value']
  ncells <- as.numeric(gsub(',', '', ncells))
  dat[which(dat$CR_ID == cr_id), 'nCells_prot'] <- ncells
  
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
p1 <- ggplot(plotdat, aes(x = Date_sort, y = nRead_RNA)) + 
  geom_dotplot(binaxis = 'y', stackdir = 'center') + ylim(c(0, (max(plotdat$nRead_RNA) * 1.1)))
p1

plotdat <- unique(dat[, c('Date_sort', 'Sequencing_saturation')])
p2 <- ggplot(plotdat, aes(x = Date_sort, y = Sequencing_saturation)) + 
  geom_dotplot(binaxis = 'y', stackdir = 'center') + ylim(c(0, (max(plotdat$Sequencing_saturation) * 1.1)))
p2

htos <- unique(dat$HTO_index_name)
### space then parentheses -> parentheses
htos <- gsub(' \\(', '\\(', htos)
### spaces, underscores to '.'
htos <- make.names(htos, allow_ = 'F')

labels <- read.table(labels_file, header = T, sep = '\t')
row.names(labels) <- labels$cell_id

#####
### Visualizing Doublet/Singlet info to inform cell filters
#####

### Use a scoring system to chose a threshold, and visualize with line graphs
### Note: A ROC curve could be a good addition
find_thresh <- function(labels, cr, val = 'ngene', value_ratio = 3){
  dub <- labels[which(labels$Assignment_simple == 'Ambiguous' & labels$CR_ID == cr), ]
  pos <- labels[which(labels$Assignment_simple == 'Positive' & labels$CR_ID == cr),]
  ### Positive's mean and standard deviation
  u <- mean(pos[which(pos$CR_ID == cr), val])
  sdv <- sd(pos[which(pos$CR_ID == cr), val])
  
  min_thresh <- floor(u + sdv)
  max_thresh <- max(dub[, val])- 1
  dat <- as.data.frame(t(as.data.frame(sapply(min_thresh:max_thresh, function(thresh){
    c(length(which(dub[, val] >= thresh)), length(which(pos[, val] >= thresh)))
  }))))
  colnames(dat) <- c('Ambiguous', 'Postives')
  row.names(dat) <- min_thresh:max_thresh
  dat$threshold <- min_thresh:max_thresh
  dat$score <- dat$Ambiguous - (value_ratio * dat$Postives)
  thresh <- dat[which(dat$score == max(dat$score))[1], 'threshold']
  best_sd <- (thresh - u)/sdv
  out <- c(thresh, best_sd)
  
  plotdat <- melt(dat, id = 'threshold')
  colnames(plotdat) <- c('threshold', 'Group', 'N_failed_filter')
  p <- ggplot(plotdat, aes(x=threshold, y = N_failed_filter, color = Group)) +
    geom_line() + 
    geom_vline(xintercept = thresh) +
    ggtitle(paste0(cr, ' best threshold: ', thresh, '(', signif(best_sd, digits = 3), ' SD) of positives'),
            subtitle = paste0('Score: Each removed Positive is worth ', value_ratio, ' removed Doublet'))
  print(p)
  return(out)
}
crs <- unique(labels$CR_ID)
pdf(out_pdf2)
threshs <- lapply(crs, function(cr) find_thresh(labels, cr))
dev.off()
threshs <- as.data.frame(t(as.data.frame(threshs)))
row.names(threshs) <- crs
colnames(threshs) <- c('Auto_Thresh', 'SDs_from_mean')

### Doublet/Singlet by rna_size(), ngene, fraction_MT, prot_size
### rna_size is log10(Matrix::colSums(rna))
### ngene is Matrix::colSums(rna > 0)
plot_dat <- labels[which(labels$Assignment_simple != 'Negative'), ]
pos_dat <- plot_dat[which(plot_dat$Assignment_simple == 'Positive'),]
### RNA size
ms <- sapply(crs, function(cr) mean(pos_dat[which(pos_dat$CR_ID == cr), 'rna_size']))
sds <- sapply(crs, function(cr) sd(pos_dat[which(pos_dat$CR_ID == cr), 'rna_size']))
rna_size_dat <- data.frame(mean = ms, sd1 = ms + sds, sd2 = (ms + 2*sds))
rna_size_dat$CR_ID <- row.names(rna_size_dat)
rna_size_dat <- melt(rna_size_dat)
### ngene
ms <- sapply(crs, function(cr) mean(pos_dat[which(pos_dat$CR_ID == cr), 'ngene']))
sds <- sapply(crs, function(cr) sd(pos_dat[which(pos_dat$CR_ID == cr), 'ngene']))
ngene_dat <- data.frame(mean = ms, sd1 = ms + sds, sd2 = (ms + 2*sds), auto = threshs$Auto_Thresh)
ngene_dat$CR_ID <- row.names(ngene_dat)
ngene_dat <- melt(ngene_dat)
### prot size
ms <- sapply(crs, function(cr) mean(pos_dat[which(pos_dat$CR_ID == cr), 'prot_size']))
sds <- sapply(crs, function(cr) sd(pos_dat[which(pos_dat$CR_ID == cr), 'prot_size']))
prot_size_dat <- data.frame(mean = ms, sd1 = ms + sds, sd2 = (ms + 2*sds))
prot_size_dat$CR_ID <- row.names(prot_size_dat)
prot_size_dat <- melt(prot_size_dat)
### fraction_MT
ms <- sapply(crs, function(cr) mean(pos_dat[which(pos_dat$CR_ID == cr), 'fraction_MT']))
sds <- sapply(crs, function(cr) sd(pos_dat[which(pos_dat$CR_ID == cr), 'fraction_MT']))
cbind(ms, ms + sds)
fraction_MT_dat <- data.frame(mean = ms, sd1 = ms + sds, sd2 = (ms + 2*sds))
fraction_MT_dat$CR_ID <- row.names(fraction_MT_dat)
fraction_MT_dat <- melt(fraction_MT_dat)

#### Violin plots pr CR of Doublet/Singlet distributions
p3 <- ggplot(plot_dat, aes(fill = Assignment_simple, x = CR_ID, y = rna_size)) +
  geom_violin() + ylab('log10(N RNA reads)')
p3 <- p3 + geom_point(data = rna_size_dat, aes(CR_ID, value, color = variable), inherit.aes = FALSE, alpha = 0.5, size = 0.5) + 
  coord_flip()
  
p4 <- ggplot(plot_dat, aes(fill = Assignment_simple, x = CR_ID, y = ngene)) +
  geom_violin() + ylab('N genes w/ > 0 counts')
p4 <- p4 + geom_point(data = ngene_dat, shape = 23, size = 3,  aes(CR_ID, value, color = variable, fill = variable), inherit.aes = FALSE, alpha = 0.5, size = 0.5) + 
  coord_flip()

p5 <- ggplot(plot_dat, aes(fill = Assignment_simple, x = CR_ID, y = prot_size)) +
  geom_violin() + ylab('log10(N protein reads)')
p5 <- p5 + geom_point(data = prot_size_dat, aes(CR_ID, value, color = variable), inherit.aes = FALSE, alpha = 0.5, size = 0.5) +
  coord_flip()

p6 <- ggplot(plot_dat, aes(fill = Assignment_simple, x = CR_ID, y = fraction_MT)) +
  geom_violin()
p6 <- p6 + geom_point(data = fraction_MT_dat, aes(CR_ID, value, color = variable), inherit.aes = FALSE, alpha = 0.5, size = 0.5) +
  coord_flip()

### If some of the sample IDs in labels are not in dat, select the cell_ids that are (in this project)
if(!all(labels[which(!labels$Assignment %in% c('Negative', 'Doublet', 'Ambiguous')), 'Assignment'] %in% dat$Sample_ID)){
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
  data_dir <- my_read_10x(id, dat, runs_dir, project_dir, cellranger = cellranger_quant, 
                          filtered = 'raw')
  
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
  ##### 
  ### Account for the '.1' appended to protein names if they are also a gene name.....
  if(any(grepl("\\.1$", row.names(prot)))){
    sub <- row.names(prot)[grepl("\\.1$", row.names(prot))]
    key <- sapply(sub, function(x) gsub('\\.1$', '', x))
    names(key) <- sub
    ### Check that the ".1" is there because the version without is in RNA
    key <- key[(key %in% row.names(rna)) & !(key %in% row.names(prot))]
    ### Replace values in row.names that match names(key) with key
    row.names(prot)[which(row.names(prot) %in% names(key))] <- key[ row.names(prot)[which(row.names(prot) %in% names(key))]]
  }  
  ### space then parentheses -> parentheses
  row.names(prot) <- gsub(' \\(', '\\(', row.names(prot))
  ### spaces, underscores to '.'
  row.names(prot) <- make.names(row.names(prot), allow_ = 'F')
  
  rm(CR_dat)
  
  print(table(labels$Assignment_simple))

  ### Keep cells by filtering on the dehashing summary column, Assignment_simple
  positive_cells <- labels[which(labels$Assignment_simple == 'Positive' & 
                                   labels$cell_id %in% colnames(rna)), , drop = F]$cell_id
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
if(length(rlist) > 1){
  for(i in 2:length(rlist)){
    rna <- cbind(rna, rlist[[i]])
  }  
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
if(length(plist) > 1){
  for(i in 2:length(plist)){
    psub <- addEmptyRows(plist[[i]], markers)
    prot <- cbind(prot, psub)
  }  
}

colnames(prot) <- unlist(lapply(plist, function(x) colnames(x)))

print('editing labels')
### Subset labels
labels <- labels[which(labels$cell_id %in% colnames(rna)),]
colnames(labels)[which(colnames(labels) == 'Assignment')] <- 'Sample_ID'
### Add useful info from dat (sample sheet) to labels (which will become Seurat metadata)
dat <- dat[, sapply(colnames(dat), function(x) length(unique(dat[, x])) > 1)]
dat <- dat[, !colnames(dat) %in% colnames(labels)]
labels[, colnames(dat)] <- NA
for(s in unique(labels$Sample_ID)){
  labels[which(labels$Sample_ID == s), colnames(dat)] <- dat[s, ]
}
### Check for cell type and SONAR info to add to meta data



print('making seurat')
### create Seurat object with cell-containing drops (min.cells is a gene filter, not a cell filter)
cseq <- Seurat::CreateSeuratObject(counts = rna, meta.data = labels, assay = "RNA", min.cells = 20, project = 'CITESeq')

print('adding prot')
### Remove htos and save as separate assay
hashdat <- prot[which((row.names(prot) %in% htos)),]
cseq[["HTO"]] <- Seurat::CreateAssayObject(counts = hashdat)
prot <- prot[which(!(row.names(prot) %in% htos)),]
cseq[["prot"]] <- Seurat::CreateAssayObject(data = prot)

print('saving')
saveRDS(cseq, file = out_rds)

pdf(out_pdf)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
VlnPlot(cseq, features = "nFeature_RNA", combine = T, log = T, pt.size = 0)
VlnPlot(cseq, features = "nCount_RNA", combine = T, log = T, pt.size = 0)
dev.off()
