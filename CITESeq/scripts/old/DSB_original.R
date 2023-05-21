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

# args = commandArgs(trailingOnly=TRUE)
# 
runs_dir <- arg[1]
covs_file <- args[2]
labels_file <- args[3]
#RDS_file <- args[4]
hto_string <- args[4]
### These next two will be the same if there was only 1 panel in the config file (string, rather than file name input)
panel_file <- args[5]
panel <- args[6]
out_file <- args[7]

runs_dir <- '/hpcdata/vrc/vrc1_data/douek_lab/Runs/211123_A00243_0139_BHLKCFDSX2/'
covs_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/All_covariates.csv'
labels_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/data/Demultiplex_labels.csv'
labels_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/data/labels.csv'
labels_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/data/labels_limitedbackground.csv'
#RDS_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/data/Unnormalized_data.RDS'
hto_string <- 'C0251,C0252,C0253,C0254,C0255,C0256,C0257,C0258,C0259,C0260'
panel_file <- ''
panel <- 'Andrews'
out_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/DSB_normalized_data.RDS'
#panel <- 'Biolegend'

quant_func <- 'count'
batch <- 'Date_Sort'
### Read covariates
covs <- read.csv(covs_file, stringsAsFactors = F, header = T,
                check.names = F, 
                colClasses = 'character',
                sep = ',')
row.names(covs) <- covs$Sample_ID

if(is.na(quant_func)) {quant_func <- 'count'}
runs_dir <- paste0(runs_dir, '/', quant_func, '_output/')
if(quant_func == 'multi'){
  runs_subdir <- '/outs/multi/count/raw_feature_bc_matrix/'
} else if(quant_func == 'count'){
  runs_subdir <- '/outs/raw_feature_bc_matrix/'
}

metadata <- read.table(labels_file, header = T, sep = ',')
table(metadata$CR_Assignment, metadata$Assignment0)
metadata <- metadata[which(metadata$Assignment != 'Doublet'),]
table(metadata$CR_Assignment, metadata$Assignment0)
batchs <- unique(metadata[, batch])
blist <- as.list(rep(NA, length(batchs)))
names(blist) <- batchs
for(bat in batchs){
  ### Filter metadata to those in this "batch"
  md <- metadata[which(metadata[, batch] == bat),]
  dat <- covs[which(covs[,batch] == bat), ]
  ### Determine threhsolds for acceptable background droplets
  
  ngmax <- mean(md[which(md$Assignment == 'Negative'), 'ngene']) + 3 * sd(md[which(md$Assignment == 'Negative'), 'ngene'])
  pmax <- mean(md[which(md$Assignment == 'Negative'), 'prot_size']) + 3 * sd(md[which(md$Assignment == 'Negative'), 'prot_size'])
  pmin <- mean(md[which(md$Assignment == 'Negative'), 'prot_size']) - 1.5 * sd(md[which(md$Assignment == 'Negative'), 'prot_size'])
  
  ### Plot log10(ngenes) v. prot_size, for selection of main peak in ADT distribution
  ggplot(md, aes(x = log10(ngene), y = prot_size )) +
    theme_bw() + 
    geom_bin2d(bins = 300) + 
    scale_fill_viridis_c(option = "C") + 
    facet_wrap(~CR_Assignment) +
    geom_hline(yintercept = pmax, linetype = "dashed", color = "red", size = 1) +
    geom_hline(yintercept = pmin, linetype = "dashed", color = "red", size = 1) + 
    geom_vline(xintercept = log10(ngmax), linetype = "dashed", color = "red", size = 1)
  
  plist <- as.list(rep(NA, length(unique(dat$CR_ID))))
  names(plist) <- unique(dat$CR_ID)
  ### For each cellranger count output, read and demultiplex
  for(id in unique(dat$CR_ID)){
    data_dir <- paste0(runs_dir,  id, runs_subdir)
    CR_raw <- Read10X(data.dir = data_dir)
    CR_filtered <- Read10X(data.dir = gsub('raw', 'filtered', data_dir))
    # define a vector of cell-containing barcodes and remove them from unfiltered data 
    stained_cells <- colnames(cells$`Gene Expression`)
    
    prot <- CR_raw[[names(CR_raw)[which(names(CR_raw) != 'Gene Expression')]]]
    rm(CR_raw)
    colnames(prot) <- paste0(colnames(prot), '_', id)
    row.names(prot) <- gsub('_totalseq', '', row.names(prot))
    
    ##### Filter cells
    ### Using Seurat's min.features
    #prot <- CreateSeuratObject(CR_raw[[names(CR_raw)[which(names(CR_raw) != 'Gene Expression')]]], min.features = 5, assay = 'prot')
    #prot <- GetAssayData(prot, assay = 'prot', slot = 'counts') %>% as.matrix()
    ### Based on presence in Demultiplex.R output, which removes those with 0 counts in either RNA or ADT data
    prot <- prot[, colnames(prot) %in% metadata$cell_id]
    ### Based on log10 number of RNA and ADT counts
    to_remove <- md[which(md$Assignment == 'Negative' & (md$prot_size <= pmin | md$prot_size >= pmax | md$ngene >= ngmax)), 'cell_id']
    md <- md[which(!md$cell_id %in% to_remove),]
    prot <- prot[, which(!colnames(prot) %in% to_remove)]
    plist[[id]] <- prot
    rm(prot)
  }
  
  ### Combine protein data
  prot <- plist[[1]]
  for(i in 2:length(plist)){
    prot <- cbind(prot, plist[[i]])
  }
  colnames(prot) <- unlist(lapply(plist, function(x) colnames(x)))
  rm(plist)
  
  
  if(any(grepl('-', row.names(prot)))) {
    row.names(prot) <- gsub('-', '_', row.names(prot))
  }
  if(any(grepl('-', md$Assignment))) {
    md$Assignment <- gsub('-', '_', md$Assignment)
  }
  ### htos can be a comma-delimited string or a regular expression matching row.names. If empty, we assume it is '^HTO_' or '^CO_'
  ### Defaults, pattern '^HTO' if that matches rows, or then '^CO_' if that does. 
  #htos <- get_hashtag_names(hto_string, row.names(prot))
  #htos <- names(table(md$Assignment))
  #htos <- htos[which(!htos %in% c('Negative', 'Doublet'))]
  
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
  
  ### Run DSB on each panel separately, but using the same negatives.
  ### Creates a Seurat object and saves to RDS, out_file
  dsb_norm_prot <- DSB_once(prot, md, hto_groups, panel_info, panel)
  blist[[bat]] <- dsb_norm_prot
}


saveRDS(dsb_norm_prot, file = out_file)

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
