library('sys')
library('Seurat')
library('stringr')
library('pheatmap')
library('ggplot2')
library('umap')
library('textshape')
library('dplyr')
library('biomaRt')
library('grid')
library('scales')
library('data.table')
library('VennDiagram')
library('scuttle')

source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R')


# source('/Volumes/VRC1_DATA/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
# source('/Volumes/VRC1_DATA/douek_lab/wakecg/sample_sheet_functions.R')
# source('/Volumes/VRC1_DATA/douek_lab/snakemakes/sc_functions.R')
# source('/Volumes/VRC1_DATA/douek_lab/snakemakes/Utility_functions.R')

if(interactive()){
  # project <- '2021614_21-002'
  # qc_name <- '2023-March'
  # qc_name <- '2023-Jan'
  # qc_name <- 'Both_celltypes'
  # sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/All_data.RDS')
  # covs_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/Sample_sheet.csv')
  # filter_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/QC_steps/RNA_filters.csv')
  # batch_name <- 'Date_sort'
  # batch_value <- '2021-11-09'
  # batch_value <- '2022-08-12'
  # #batch_value <- '2021-12-02'
  # base_dir <- '/hpcdata/vrc/vrc1_data/'
  
  base_dir <- '/Volumes/VRC1_DATA/'
  #base_dir <- '/hpcdata/vrc/vrc1_data/'
  project <- '2022619_857.3b'
  qc_name <- 'QC_SampleName_pass'
  sdat_file <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/data/All_data.RDS')
  covs_file <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/Sample_sheet.csv')
  filter_file <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/QC_steps/step3_cell_and_feature_filters.csv')
  batch_name <- 'Sample_Name'
  batch_value <- '15C225_w31_Probe+'
  gtf_file <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/data/gtf.RDS')
  out_rds <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/', batch_value, '/RNA_cell_filtered.RDS')
  out_pdf <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/', batch_value, '/MT_nFeatures.pdf')
  out_txt1 <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/', batch_value, '/Excluded_genes.txt')
  out_txt2 <- paste0(base_dir, '/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/', batch_value, '/Filtered_Ncells.txt')
  
  
  
  setwd(paste0(base_dir, 'douek_lab/projects/RNASeq/2022619_857.3b'))
} else {
  args = commandArgs(trailingOnly=TRUE)
  
  sdat_file <- args[1]
  #checkpt_file <- args[2]
  filter_file <- args[2]
  covs_file <- args[3]
  exclude_file <- args[4]
  gtf_file <- args[5]
  batch_value <- args[6]
  batch_name <- args[7]
  out_rds <- args[8]
  out_pdf <- args[9]
  out_txt2 <- args[10]
}

print(sdat_file)
print(filter_file)
print(covs_file)
print(gtf_file)
print(batch_value)
print(batch_name)
print(out_rds)
print(out_pdf)
print(out_txt2)

filters <- read.table(filter_file, header = T, sep = ',')
nFeature_RNA <- filters[which(filters$type == 'cell' & filters$feature == 'nFeature_RNA'), 'value']
nCount_RNA <- filters[which(filters$type == 'cell' & filters$feature == 'nCount_RNA'), 'value']
fraction_MT <- filters[which(filters$type == 'cell' & filters$feature == 'fraction_MT'), 'value']

cell_filters <- filters[which(filters$type == 'cell'), 'feature'] 
### Set defaults
# if(is.na(fraction_cells) | fraction_cells == ''){
#   fraction_cells <- '<0.03'
# }
if(length(nFeature_RNA) == 0){
  nFeature_RNA <- '20,3(SD)'
}
if(length(nCount_RNA) == 0){
  nCount_RNA <- '100,999999(SD)'
}
#if(is.na(fraction_MT) | fraction_MT == ''){
if(length(fraction_MT) == 0){
  fraction_MT <- '0,0.4'
}

sdat <- readRDS(sdat_file)
DefaultAssay(sdat) <- 'RNA'

#thresh <- threshold_string_seurat(sdat, '100,999999(SD)', 'nCount_RNA', T)
#VlnPlot(sdat, 'nCount_RNA', pt.size = 0, log = T)  + geom_hline(yintercept = thresh[1], color = 'red')

### Subset to the input batch
print(paste0('Subsetting to input batch, ', batch_name))
print(table(sdat@meta.data[, batch_name]), useNA = 'ifany')
keep_cells <- colnames(sdat)[which(sdat@meta.data[, batch_name] == batch_value)]
sdat <- subset(sdat, cells = keep_cells)
#saveRDS(sdat, file = gsub('/RNA_cell_filtered.RDS', '.RDS', out_rds))
#sdat <- readRDS(gsub('/RNA_cell_filtered.RDS', '.RDS', out_rds))
print(table(sdat@meta.data[, batch_name]), useNA = 'ifany')

#sdat <- readRDS(gsub('.RDS', '_prior.RDS', out_rds))
# saveRDS(sdat, file = gsub('.RDS', '_prior.RDS', out_rds))
### To create a test set for interactive use
# csub <- subset(x = sdat, downsample = 1000)
# print(dim(csub))
# saveRDS(csub, file = '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/Downsample_0.8/Test_subset.RDS')

print(paste0('Reading gtf file ', gtf_file))
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
  #id_type <- names(which.max(sapply(gtf_cols, function(col) sum(row.names(sdat) %in% gtf[, col]))))
  id_type <- get_id_name(row.names(sdat), gtf)
}

if(file.info(exclude_file)[, 'size'] > 0){
  exclude_gene_names <- read.table(exclude_file)[,1]
} else{
  exclude_gene_names <- c()
}
print('Determing which genes are mitochondiral')
### Determine which genes are mitochondrial by gtf if it is input. Otherwise, by gene name pattern '^MT'.
if(is.data.frame(gtf)){
  sdat[['MT_sum']] <- sum_chromosome(mtx = sdat@assays$RNA@counts, chr ='MT', gtf = gtf)
} else {
  mts <- row.names(sdat@assays$RNA)[grepl('^MT-', row.names(sdat@assays$RNA))]
  sdat[['MT_sum']] <- colSums(sdat@assays$RNA@counts[mts, ])
}
### Calculate mitochondial percentage
sdat[['fraction_MT']] <- sdat[['MT_sum']]/sdat[['nCount_RNA']]

# cor(sdat$nFeature_RNA, sdat$fraction_MT)
# plot(sdat$nFeature_RNA, sdat$fraction_MT, xlab = 'N genes', ylab = 'Sum MT TPM')
# #FeatureScatter(sdat, feature1 = "nFeature_RNA", feature2 = "MT_sum", group.by = "sample")
# 
# plot(sdat$nFeature_RNA, sdat$fraction_MT, xlab = 'N genes', ylab = 'Sum MT TPM') +
#   abline(v = NFeat_thresh[1], col = 'red') +
#   abline(v = NFeat_thresh[2], col = 'red') +
#   abline(h = MT_thresh[2], col = 'red')

print('Creating text file with summary of filter numbers')
#filter_names <- c('nFeature_RNA', 'nCount_RNA', 'fraction_MT')
filter_names <- rep(cell_filters, each = 2)
names(filter_names) <- paste(rep(cell_filters, each = 2), c('low', 'high'), sep = '_')

txt <- as.data.frame(matrix(ncol = 3, nrow = length(filter_names) + 2))
row.names(txt) <- c('original', names(filter_names), 'remaining')
colnames(txt) <- c('step', 'N_fail', 'N_by_step')
txt[, 'step'] <- 0:(length(filter_names) + 1)

sdat_full <- sdat
### Record raw number of cells that fail threshold
cell_list <- vector(mode = 'list', length = length(filter_names))
names(cell_list) <- names(filter_names)

sdat$Sample_ID <- factor(sdat$Sample_ID, levels = unique(sdat$Sample_ID)) 
sdat$Sample_Name <- factor(sdat$Sample_Name, levels = unique(sdat$Sample_Name))  

### N cells per sample
samp_ns <- as.character(sapply(levels(sdat$Sample_Name), function(b) length(which(sdat@meta.data[, 'Sample_Name'] == b))))
pdf(out_pdf)
### Loop over each input cell filter
for(filter in unique(filter_names)){
  print(paste0(filter, " filter"))
  ### Get array of lower,upper bound thresholds from the comma-delimited string. Options are flat values, or if in values of standard deviation, 
  ### with a trailing (SD)
  thresh_str <- filters[which(filters$type == 'cell' & filters$feature == filter), 'value']
  threshs <- threshold_string_seurat(sdat, thresh_str, filter, T)
  ### Lower bound
  cells1 <-  row.names(sdat@meta.data[which(sdat@meta.data[, filter] < threshs[1]), ])
  cell_list[[names(filter_names[which(filter_names == filter)])[1]]] <- cells1
  txt[names(filter_names[which(filter_names == filter)])[1], 'N_fail'] <- length(cells1)
  ### Upper bound
  cells2 <- row.names(sdat@meta.data[which(sdat@meta.data[, filter] > threshs[2]), ])
  cell_list[[names(filter_names[which(filter_names == filter)])[2]]] <- cells2
  txt[names(filter_names[which(filter_names == filter)])[2], 'N_fail'] <- length(cells2)
  
  if(all(sdat@meta.data[, filter] <= 1 & sdat@meta.data[, filter] >= 0)){
    do_log <- F
  } else{
    do_log <- T
  }
  #sids <- names(table(sdat$Sample_ID)[which(table(sdat$Sample_ID) < 50)])
  #too_few_cells <- row.names(sdat@meta.data[which(!sdat@meta.data$Sample_ID %in% sids),])
  #sdat1 <- subset(sdat, cells = too_few_cells)
  
  ### Violin plot, filter by Sample_Name 
  text_y <- max(sdat@meta.data[, filter]) * 0.92
  p <- VlnPlot(sdat, features = filter, group.by = 'Sample_Name', log = do_log, pt.size = 0) + theme(legend.position = 'none') + 
    annotate('text', x = 0.56, y = text_y, size = 3, label = 'N') +
    annotate('text', x = 1:length(unique(sdat@meta.data[, 'Sample_Name'])), y = text_y, label = samp_ns, size = 3, angle = 30) + 
    labs(subtitle = paste0(signif((length(cells1) + length(cells2)) / dim(sdat)[2], digits = 4), ' of cells fail')) +  
    coord_flip() 
  if(length(cells1) > 0){
    p <- p + geom_hline(yintercept = threshs[1], color = 'red') 
    p <- p + annotate('rect', ymin = min(sdat@meta.data[, filter]), ymax = threshs[1], xmin = 0, xmax = length(unique(sdat$Sample_Name)), alpha = 0.1, fill = 'red') 
  }
  if(length(cells2) > 0){
    p <- p + geom_hline(yintercept = threshs[2], color = 'red')
    p <- p + annotate('rect', ymin = threshs[2], ymax = max(sdat@meta.data[, filter]), xmin = 0, xmax = length(unique(sdat$Sample_Name)), alpha = 0.1, fill = 'red') 
  }
  print(p)
}

### Venn diagram showing overlap of multiple filters
filter_plot_names <- gsub('_high', ' upper bound', gsub('_low', ' lower bound', names(filter_names)))
### Remove those that removed no cells
filter_plot_names <- filter_plot_names[which(sapply(cell_list, function(x) length(x)) > 0)]
cell_list <- cell_list[which(sapply(cell_list, function(x) length(x)) > 0)]
if (length(filter_plot_names) > 1) {
  plot.new()
  my_venn(names = filter_plot_names, cell_list)
  dev.off()
}

print('Subsetting the seurat object')
### Recording the number of cells removed at each step
current_size <- dim(sdat_full)[2]
for(filter in unique(filter_names)){
  print(filter)
  thresh_str <- filters[which(filters$type == 'cell' & filters$feature == filter), 'value']
  threshs <- threshold_string_seurat(sdat, thresh_str, filter, T)
  ##### Do subsetting of each threshold, and record how many are removed at that step
  ### Lower bound
  print('lower')
  cells_to_keep <- row.names(sdat@meta.data[which(sdat@meta.data[, filter] >= threshs[1]), ])
  sdat <- subset(sdat, cells = cells_to_keep)
  txt[names(filter_names[which(filter_names == filter)])[1], 'N_by_step'] <- (current_size - dim(sdat)[2])
  print(paste0(filter, ' lower bound of ', threshs[1], ' removes ', (current_size - dim(sdat)[2]), ' cells'))
  current_size <- dim(sdat)[2]
  
  ### Upper bound
  print('upper')
  cells_to_keep <- row.names(sdat@meta.data[which(sdat@meta.data[, filter] <= threshs[2]), ])
  sdat <- subset(sdat, cells = cells_to_keep)
  txt[names(filter_names[which(filter_names == filter)])[2], 'N_by_step'] <- (current_size - dim(sdat)[2])
  print(paste0(filter, ' upper bound of ', threshs[2], ' removes ', (current_size - dim(sdat)[2]), ' cells'))
  current_size <- dim(sdat)[2]
}
#plot(sdat$nFeature_RNA, sdat$fraction_MT, xlab = 'N genes', ylab = 'Sum MT TPM')
#dev.off()

txt['original', 'N_fail'] <- dim(sdat_full)[2]
txt['original', 'N_by_step'] <- NA
txt['remaining', 'N_fail'] <- dim(sdat)[2]
txt['remaining', 'N_by_step'] <- NA
txt$freq_fail <- txt$N_fail /  txt[1,2]
txt$freq_by_step <- txt$N_by_step /  txt[1,2]
write.table(txt, out_txt2, sep = ',', quote = F, row.names = T, col.names = T)
#sdat_removed <- subset(sdat_full, cells = colnames(sdat_full)[!colnames(sdat_full) %in% colnames(sdat)])

### Fraction kept by each sample. The first two lines are necessary in the case that an entire sample is filtered
# temp <- table(sdat$sample)
# temp[unique(sdat_full$sample)[which(!unique(sdat_full$sample) %in% unique(sdat$sample))]] <- 0
# temp/table(sdat_full$sample)[names(temp)]
print('Saving new seurat object')
saveRDS(sdat, file = out_rds)

