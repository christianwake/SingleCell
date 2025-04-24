library('sys')
library('Seurat')
library('stringr')
library('pheatmap')
library('ggplot2')
#library('umap')
#library('textshape')
library('dplyr')
library('biomaRt')
library('grid')
library('scales')
library('data.table')
library('VennDiagram')
#library('scuttle')

source('/data/vrc_his/douek_lab/snakemakes/sc_functions.R')
source('/data/vrc_his/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/Utility_functions.R')

if(interactive()){
  project <- '2021600_kristin'
  qc_name <- 'Run2025-04-10'
  batch_name <- 'FC_ID'
  batch_value <- 'HF5V3DRXY'
  
  # project <- '2024615_Boswell'
  # qc_name <- 'Run2025-04-10'
  # batch_name <- 'FC_ID'
  # batch_value <- 'FC1'
 
  base_dir <- '/data/vrc_his/'
  in_file <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/data/fastq/trimmed/multiqc_input_files.txt')
  out_file <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/data/fastq/filtered/multiqc_input_files.txt')
  
  covs_file <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/data/Covariates_QC_metrics.csv')
  filter_file <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, 
                        '/QC_steps/Cell_filters.csv')
} else {
  args = commandArgs(trailingOnly=TRUE)
  
  in_file <- args[1]
  filter_file <- args[2]
  covs_file <- args[3]
  out_file1 <- args[4]
  out_file2 <- args[5]
}

print(filter_file)
print(covs_file)

intxt <- read.table(in_file)
intxt$SeqRep_ID<- sapply(intxt[,1], function(x) strsplit(x, '/')[[1]][5])

### Read in covariates
covs <- read.csv(covs_file, stringsAsFactors = F, header = T,
                        check.names = F, 
                        colClasses = 'character',
                        sep = ',')
row.names(covs) <- covs$Cell_ID

filters <- read.table(filter_file, header = T, sep = ',')
### Remove NA values
filters <- filters[which(!is.na(filters$value)),]

nFeature_RNA <- filters[which(filters$type == 'cell' & filters$feature == 'nFeature_RNA'), 'value']
nCount_RNA <- filters[which(filters$type == 'cell' & filters$feature == 'nCount_RNA'), 'value']
fraction_MT <- filters[which(filters$type == 'cell' & filters$feature == 'fraction_MT'), 'value']

if(length(nCount_RNA) == 0 & length(nFeature_RNA) == 0){
  nCount_RNA <- '100,999999(SD)'
  filters[length(row.names(filters) + 1), ] <- c('cell', 'nCount_RNA', '100,99999999(SD)')
}
if(length(fraction_MT) == 0){
  fraction_MT <- '0,0.4'
  filters[length(row.names(filters) + 1), ] <- c('cell', 'fraction_MT', '0,0.4')
}
cell_filters <- filters[which(filters$type == 'cell' & (!is.na(filters$value))), 'feature'] 

cell_filters <- cell_filters[which(cell_filters %in% colnames(covs))]


### Record raw number of cells that fail threshold
cell_list <- vector(mode = 'list', length = length(cell_filters))
names(cell_list) <- cell_filters

### Loop over each input cell filter
for(filter in cell_filters[which(cell_filters == 'Fastqc_Adapter_Content')]){
  ### Determine if filter is quantitative or qualitative
  quant <- T
  if(any(is.na(suppressWarnings(as.numeric(covs[, filter]))))){
    quant <- F
  }
  covs[, filter]
  if(quant){
    covs[, filter] <- as.numeric(covs[, filter])
    print(paste0(filter, " filter"))
    ### Get array of lower,upper bound thresholds from the comma-delimited string. Options are flat values, or if in values of standard deviation, 
    ### with a trailing (SD)
    thresh_str <- filters[which(filters$type == 'cell' & filters$feature == filter), 'value']
    threshs <- threshold_string_seurat(covs, thresh_str, filter, T)
    ### Lower bound
    cells1 <-  row.names(covs[which(covs[, filter] < threshs[1]), ])
    ### Upper bound
    cells2 <- row.names(covs[which(covs[, filter] > threshs[2]), ])
    cell_list[[cell_filters[which(cell_filters == filter)][1]]] <- c(cells1, cells2)
  }else{
    filt <- filters[which(filters$type == 'cell' & filters$feature == filter), 'value']
    exclusive <- ifelse(substr(filt, 1, 1) == '!', T, F)
    if(exclusive){
      filt <- substr(filt, 2, nchar(filt))
      filt <- strsplit(filt, ',')[[1]]
      cells <- row.names(covs[which(covs[, filter] %in% filt),])
    } else{
      filt <- strsplit(filt, ',')[[1]]
      cells <- row.names(covs[which(!covs[, filter] %in% filt),])
    }
    cell_list[[cell_filters[which(cell_filters == filter)][1]]] <- cells
  }
}
cells <- unlist(cell_list)
reps <- unlist(lapply(covs[cells, 'SeqRep_ID'], function(r) strsplit(r, ' & ')[[1]]))

outtxt <- intxt[which(!intxt$SeqRep_ID %in% reps), 1, drop = F]

write.table(outtxt1, file = out_file, row.names = F, col.names = F, quote = F)

### Loop over each input cell filter
for(filter in cell_filters){
  ### Determine if filter is quantitative or qualitative
  quant <- T
  if(any(is.na(suppressWarnings(as.numeric(covs[, filter]))))){
    quant <- F
  }
  covs[, filter]
  if(quant){
    covs[, filter] <- as.numeric(covs[, filter])
    print(paste0(filter, " filter"))
    ### Get array of lower,upper bound thresholds from the comma-delimited string. Options are flat values, or if in values of standard deviation, 
    ### with a trailing (SD)
    thresh_str <- filters[which(filters$type == 'cell' & filters$feature == filter), 'value']
    threshs <- threshold_string_seurat(covs, thresh_str, filter, T)
    ### Lower bound
    cells1 <-  row.names(covs[which(covs[, filter] < threshs[1]), ])
    ### Upper bound
    cells2 <- row.names(covs[which(covs[, filter] > threshs[2]), ])
    cell_list[[cell_filters[which(cell_filters == filter)][1]]] <- c(cells1, cells2)
  }else{
    filt <- filters[which(filters$type == 'cell' & filters$feature == filter), 'value']
    exclusive <- ifelse(substr(filt, 1, 1) == '!', T, F)
    if(exclusive){
      filt <- substr(filt, 2, nchar(filt))
      filt <- strsplit(filt, ',')[[1]]
      cells <- row.names(covs[which(covs[, filter] %in% filt),])
    } else{
      filt <- strsplit(filt, ',')[[1]]
      cells <- row.names(covs[which(!covs[, filter] %in% filt),])
    }
    cell_list[[cell_filters[which(cell_filters == filter)][1]]] <- cells
  }
}
cells <- unlist(cell_list)
reps <- unlist(lapply(covs[cells, 'SeqRep_ID'], function(r) strsplit(r, ' & ')[[1]]))

outtxt <- intxt[which(!intxt$SeqRep_ID %in% reps), 1, drop = F]

write.table(outtxt2, file = out_file, row.names = F, col.names = F, quote = F)
