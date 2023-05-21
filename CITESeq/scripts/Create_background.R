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
library('VennDiagram')

source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/sample_sheet_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R')

# source('/Volumes/VRC1_DATA/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
# source('/Volumes/VRC1_DATA/douek_lab/wakecg/sample_sheet_functions.R')
# source('/Volumes/VRC1_DATA/douek_lab/snakemakes/sc_functions.R')
# source('/Volumes/VRC1_DATA/douek_lab/snakemakes/Utility_functions.R')

if(interactive()){
  project <- '2021614_21-002'
  base_dir <- '/hpcdata/vrc/vrc1_data/douek_lab/'
  qc_name <- '2023-Jan'
  runs_dir <- '/hpcdata/vrc/vrc1_data/douek_lab/Runs/'

    # project <- '2022619_857.3b'
  # qc_name <- '2023-Mar'
  # runs_dir <- '/hpcdata/vrc/vrc1_data/douek_lab/Runs/'
  # covs_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022619_857.3b/Sample_sheet.csv'
  # labels_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022619_857.3b/data/Cell_data.csv'
  # qc_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022619_857.3b/QC_steps/DSB_background.csv'
  # out_rds <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022619_857.3b/results/', qc_name, '/Background.RDS')
  # out_pdf <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022619_857.3b/results/', qc_name, '/Background.pdf')
  
  # base_dir <- '/Volumes/VRC1_DATA/douek_lab/'
  # #base_dir <- '/hpcdata/vrc/vrc1_data/douek_lab/'
  # project <- '2022619_857.3b'
  # qc_name <- 'QC_first_pass'
  # runs_dir <- '/Volumes/VRC1_DATA/douek_lab/Runs/'
  
  covs_file <- paste0(base_dir, '/projects/RNASeq/', project, '/Sample_sheet.csv')
  labels_file <-  paste0(base_dir, '/projects/RNASeq/', project, '/data/Cell_data.csv')
  qc_file <-  paste0(base_dir, '/projects/RNASeq/', project, '/QC_steps/DSB_background.csv')
  gtf_file <- paste0(base_dir, 'projects/RNASeq/', project, '/data/gtf.RDS')
  out_rds <- paste0(base_dir, '/projects/RNASeq/', project, '/results/', qc_name, '/Background.RDS')
  out_pdf <- paste0(base_dir, '/projects/RNASeq/', project, '/results/', qc_name, '/Background.pdf')
  # 
  # setwd(paste0(base_dir, '/projects/RNASeq/2022619_857.3b'))

  
} else{
  args = commandArgs(trailingOnly=TRUE)
  
  runs_dir <- args[1]
  project <- args[2]
  covs_file <- args[3]
  labels_file <- args[4]
  qc_file <- args[5]
  gtf_file <- args[6]
  out_rds <- args[7]
  out_pdf <- args[8]
}

cellranger_quant <- c('multi', 'count')
print(out_rds)
project_dir <- paste0("/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/", project)
# project_dir <- paste0("/Volumes/VRC1_DATA/douek_lab/projects/RNASeq/", project)

### Read covariates
dat <- read.csv(covs_file, stringsAsFactors = F, header = T,
                check.names = F, 
                colClasses = 'character',
                sep = ',')
row.names(dat) <- dat$Sample_ID

htos <- unique(dat$HTO_index_name)

### Dehash info
labels <- read.table(labels_file, header = T, sep = ',')
#labels$cell_id <- sapply(labels$cell_id, function(x) paste0(strsplit(x, '_')[[1]][2], '_', gsub('-1', '', strsplit(x, '_')[[1]][1])))

### If created by the snakemake pipleine, this should be T, but we want to account for the situation where it is F
hash_includes_bg <- ('Background' %in% names(table(labels$Assignment_CR)))

qc <- read.table(qc_file, header = T, sep = ',')

if(gtf_file == '' | is.na(gtf_file)){
  gtf <- NA
} else{
  if(grepl('\\.gtf', gtf_file)){
    gtf <- read_gtf(gtf_file, 'gene', c('gene_name', 'gene_id')) 
  } else{
    gtf <- readRDS(gtf_file)
  }
  row.names(gtf) <- gtf$gene_id
  ### Determine whether 'gene_name' or 'gene_id' better matches those in the seurat object
  gtf_cols <- c('gene_name', 'gene_id')
}

### MODIFCATION - Set defaults if it is empty or not input
### Filter values from DSB documentation
# prot_size_max <- 3
# prot_size_min <- 1.5
# ngene_max <- 100
prot_size_max <- qc[which(qc$feature == 'prot_size_max'), 'value']
prot_size_min <- qc[which(qc$feature == 'prot_size_min'), 'value']
ngene_max <- qc[which(qc$feature == 'ngene_max'), 'value']
print(paste0('ngene_max:', ngene_max))
### Subset to CRID
labels <- labels[which(labels$CR_ID %in% unique(dat$CR_ID)),]
dim(labels)

#labels$cell_id <- sapply(labels$X, function(x) paste0(strsplit(x, '_')[[1]][2], '_',strsplit(x, '_')[[1]][1]))
labels <- labels[, c('cell_id', 'CR_ID', 'Assignment_CR', 'Assignment_simple', 'rna_size', 'ngene', 'prot_size', 'fraction_MT')]
### Account for weird "multi_" prefacing cell_id if processed by cellranger multi
#labels$cell_id <- sapply(labels$cell_id, function(x) gsub('multi_', '', x))

plist <- as.list(rep(NA, length(unique(dat$CR_ID))))
names(plist) <- unique(dat$CR_ID)
### For each cellranger count output, read and filter based on info in 'labels'
for(id in unique(dat$CR_ID)){
  ### Get desired directory of data
  data_dir <- my_read_10x(id, dat, runs_dir, project_dir, cellranger = cellranger_quant, filtered = 'raw')

  CR_dat <- Read10X(data.dir = data_dir)
  ### Protein data
  prot <- CR_dat[[names(CR_dat)[which(names(CR_dat) != 'Gene Expression')]]]
  ### Account for weird "multi_" prefacing cell_id if processed by cellranger multi
  colnames(prot) <- sapply(colnames(prot), function(cn) strsplit(cn, '_')[[1]][2])
  ### Prepend Cellranger ID, so they can be combined outside the loop
  colnames(prot) <- paste0(id, '_', colnames(prot))
  ### Remove the trailing "-1"
  colnames(prot) <- gsub('-1$', '', colnames(prot))
  
  
  ### If dehashing labels don't contain the CellRanger Background calls (as they would if the snakemake pipeline was used for dehashing), 
  ### because the dehashing filtered them, we need that info back
  if(!hash_includes_bg){
    ### Cellranger droplets that ARENT in the dehash info of this CR_ID
    CR_background <- colnames(prot)[which(!colnames(prot) %in% labels[which(labels$CR_ID == id), 'cell_id'])]
    ### Instantiate background labels data frame
    labels_bg <- as.data.frame(matrix(ncol = length(colnames(labels)), nrow = length(CR_background)))
    colnames(labels_bg) <- colnames(labels)
    row.names(labels_bg) <- CR_background
    labels_bg[, 'cell_id'] <- CR_background
    labels_bg[, 'CR_ID'] <- id
    labels_bg[, 'Assignment_CR'] <- 'Background'
    labels_bg[, 'Assignment_simple'] <- 'Cellranger_Background'
    
    ### Calculate prot_size and ngene since we will filter on those
    rna <- CR_dat$`Gene Expression`
    ### Account for weird "multi_" prefacing cell_id if processed by cellranger multi
    colnames(rna) <- sapply(colnames(rna), function(cn) strsplit(cn, '_')[[1]][2])
    ### Prepend Cellranger ID, so they can be combined outside the loop
    colnames(rna) <- paste0(id, '_', colnames(rna))
    ### Remove the trailing "-1"
    colnames(rna) <- gsub('-1$', '', colnames(rna))
    
    # create metadata of droplet QC stats used in standard scRNAseq processing
    rna_size <- log10(Matrix::colSums(rna))
    ngene <- Matrix::colSums(rna > 0)
    prot_size <- log10(Matrix::colSums(prot) + 1)
    ### Read from gtf if it is input
    if(is.data.frame(gtf)){
      ### Get the names of genes on the input chromosome
      mtgene <- unlist(lapply(c('gene_id', 'gene_name'), function(id) unique(gtf[which(gtf$chr == 'MT'), id])))
      ### Subset to those in the mtx
      mtgene <- mtgene[mtgene %in% row.names(rna)]
    } else{
      mtgene <- grep(pattern = "^MT-", rownames(rna), value = TRUE)
    }
    propmt <- Matrix::colSums(rna[mtgene, ]) / Matrix::colSums(rna)
    rm(rna)
    ### Add to background_drops
    labels_bg[, 'rna_size'] <- rna_size[CR_background]
    labels_bg[, 'ngene'] <- ngene[CR_background]
    labels_bg[, 'prot_size'] <- prot_size[CR_background]
    labels_bg[, 'fraction_MT'] <- propmt[CR_background]
    ### Append to labels
    labels <- rbind(labels, labels_bg)
    ### Keep background droplets by filtering on the dehashing summary column, Assignment_simple, and only those in this CR sample
    background_drops <- labels[which((labels$Assignment_simple %in% c('Negative', 'Cellranger_Background') | labels$Assignment_CR == 'Background') & 
                                       labels$cell_id %in% colnames(prot)), , drop = F]
    
  } else{
    ### Keep background droplets by filtering on the dehashing summary column, Assignment_simple, and only those in this CR sample
    background_drops <- labels[which(labels$Assignment_simple %in% c('Negative', 'Cellranger_Background') & labels$Assignment_CR == 'Background' & 
                                       labels$cell_id %in% colnames(prot)), , drop = F]
  }
  if(dim(background_drops)[1] == 0){
    warning('No background drops to filter')
  }
  rm(CR_dat)
  
  #labels$cell_id <- paste0(labels$cell_id, '-1_', labels$CR_ID)

  ### Now filter on protein library size and RNA content
  background_drops <- background_drops[background_drops$prot_size > prot_size_min & background_drops$prot_size < prot_size_max & background_drops$ngene < ngene_max, ]
  prot <- prot[, which(colnames(prot) %in% background_drops$cell_id)]
  print(paste0(id, ' NAs: ', sum(is.na(prot))))
  
  plist[[id]] <- prot
  rm(prot)
}

### Combine prot data and add empty rows when fewer than max hashes  
markers <- unique(unlist(lapply(plist, function(prot) row.names(prot))))
addEmptyRows <- function(prot, markers){
  ### Markers not in this CR ID
  toAdd <- markers[which(!markers %in% row.names(prot))]
  myMatrix <- matrix(ncol = ncol(prot), nrow = length(toAdd))
  row.names(myMatrix) <- toAdd
  prot <- rbind(prot, myMatrix)
  ### Conform row order
  prot <- prot[markers, ]
  return(prot)
}
### adding empty rows 
prot <- addEmptyRows(plist[[1]], markers)
for(i in 2:length(plist)){
  psub <- addEmptyRows(plist[[i]], markers)
  prot <- cbind(prot,psub)
}
colnames(prot) <- unlist(lapply(plist, function(x) colnames(x)))
row.names(prot) <- gsub('-totalseq', '', row.names(prot))
row.names(prot) <- gsub('_totalseq', '', row.names(prot))

### Remove htos
print(dim(prot))

prot <- prot[which(!(row.names(prot) %in% htos)),]

if(sum(is.na(prot)) > 0 ){
  warning(paste0('There should not me any NAs, but there are ', sum(is.na(prot))))
}

print(dim(prot))
### create Seurat object with cell-containing drops (min.cells is a gene filter, not a cell filter)
#cseq <- Seurat::CreateSeuratObject(counts = rna, meta.data = labels, assay = "RNA", min.cells = 20)

#cseq[["prot"]] = Seurat::CreateAssayObject(data = prot)

saveRDS(prot, file = out_rds)

print(table(labels$Assignment_CR, labels$Assignment_simple))

#### Cell Ranger negatives
original <- length(labels[labels$Assignment_CR != 'Cells',]$cell_id)
### Cell Ranger negatives matching the input QC criteria
remaining <- length(labels[labels$Assignment_CR != 'Cells' & labels$prot_size > prot_size_min & labels$prot_size < prot_size_max & 
                             labels$ngene < ngene_max, ]$cell_id)
print(paste0(signif((original- remaining)/original, 3), ' of background droplets filtered'))

p1 <- ggplot(labels, aes(x = log10(ngene), y = prot_size )) +
  theme_bw() + 
  geom_bin2d(bins = 80) + 
  scale_fill_viridis_c(option = "C") + 
  facet_wrap(~Assignment_CR) +
  geom_hline(yintercept = prot_size_max, linetype = "dashed", color = "red", size = 0.8) +
  geom_hline(yintercept = prot_size_min, linetype = "dashed", color = "red", size = 0.8) + 
  geom_segment(aes(x = log10(ngene_max), xend = log10(ngene_max), y = prot_size_min, yend = prot_size_max), size = 0.8, colour = "red", linetype = "dashed") +
  labs(title = 'Cell Ranger calls', subtitle = paste0(signif((original- remaining)/original, 3), ' of background droplets filtered'))

#### Dehash negatives
original <- length(labels[labels$Assignment_simple == 'Negative',]$cell_id)
### Dehash negatives matching the input QC criteria
remaining <- length(labels[labels$Assignment_simple == 'Negative' & labels$prot_size > prot_size_min & 
                             labels$prot_size < prot_size_max & labels$ngene < ngene_max, ]$cell_id)

p2 <- ggplot(labels, aes(x = log10(ngene), y = prot_size )) +
  theme_bw() + 
  geom_bin2d(bins = 80) + 
  scale_fill_viridis_c(option = "C") + 
  facet_wrap(~Assignment_simple) +
  geom_hline(yintercept = prot_size_max, linetype = "dashed", color = "red", size = 0.8) +
  geom_hline(yintercept = prot_size_min, linetype = "dashed", color = "red", size = 0.8) + 
  geom_segment(aes(x = log10(ngene_max), xend = log10(ngene_max), y = prot_size_min, yend = prot_size_max), size = 0.8, colour = "red", linetype = "dashed") +
  #geom_text(x = 1, y = 0, label = signif((original- remaining)/original, 3)) +
  labs(title = 'Dehash calls', subtitle = paste0(signif((original- remaining)/original, 3), ' of negative droplets filtered'))

p3 <- ggplot(labels, aes(x = log10(ngene), y = prot_size, color = fraction_MT, size = 0.2)) +
  theme_bw() +
  geom_point() +
  scale_fill_viridis_c(option = "C") +
  facet_wrap(~Assignment_CR)# +

print('here1')
print('Assignment_CR' %in% colnames(labels))
print(labels$Assignment_CR)

forvenn <- c(lapply(unique(labels$Assignment_CR), function(x) labels[which(labels$Assignment_CR == x), 'cell_id']), 
             lapply(unique(labels$Assignment_simple), function(x) labels[which(labels$Assignment_simple == x), 'cell_id']))

print(length(forvenn))
print(forvenn[1])

pdf(out_pdf)
unique(labels$Assignment_CR)
unique(labels$Assignment_simple)
p4 <- my_venn(c(unique(labels$Assignment_CR), unique(labels$Assignment_simple)), forvenn)

print('here2')
p1
p2
#p3
dev.off()

png(gsub('.pdf', '.png', out_pdf))
p3
dev.off()
