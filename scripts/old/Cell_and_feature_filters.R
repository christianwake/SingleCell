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
library('vegan')
library('data.table')

source(paste0('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R') )
source(paste0('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R'))
source(paste0('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R'))


if(interactive()){

#  rootPath <- "/Volumes/VRC1_DATA/douek_lab"

  # project <- '2021617_mis-c'

  project <- '2021600_kristin'
  qc_name <- 'Run2023-05-14'
  sdat_step <- paste0('data/All_data.RDS')

  # project <- '2022620_857.1'
  # qc_name <- '2022-12-09'
  # sdat_step <- paste0('data/All_data.RDS')
  
  # project <- '2022619_857.3b'
  # qc_name <- 'QC_first_pass'
  # sdat_step <- paste0('data/All_data.RDS')

  sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/', sdat_step)
  filter_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/QC_steps/step3_cell_and_feature_filters.csv')
  gtf_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/gtf.RDS')
  out_rds <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/PostQC3.RDS')
  out_pdf <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Cell_filters.pdf')
  out_txt0 <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Filters_Nfeatures.txt')
  out_txt1 <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Excluded_genes.txt')
  out_txt2 <-  paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Filters_Ncells.txt')
 
 # sdat_file <- paste0(rootPath, '/projects/RNASeq/', project, '/', sdat_step)
 # filter_file <- paste0(rootPath,'/projects/RNASeq/', project, '/QC_steps/step3_cell_and_feature_filters.csv')
  # gtf_file <- paste0(rootPath,'/projects/RNASeq/', project, '/data/gtf.RDS')
  # out_rds <- paste0(rootPath,'/projects/RNASeq/', project, '/results/', qc_name, '/PostQC3.RDS')
  # out_pdf <- paste0(rootPath,'/projects/RNASeq/', project, '/results/', qc_name, '/Cell_filters.pdf')
  # out_txt0 <- paste0(rootPath,'/projects/RNASeq/', project, '/results/', qc_name, '/Filters_Nfeatures.txt')
  # out_txt1 <- paste0(rootPath,'/projects/RNASeq/', project, '/results/', qc_name, '/Excluded_genes.txt')
  # out_txt2 <-  paste0(rootPath,'/projects/RNASeq/', project, '/results/', qc_name, '/Filters_Ncells.txt')
} else{
  args = commandArgs(trailingOnly=TRUE)

  # rootPath <- "/hpcdata/vrc/vrc1_data/douek_lab"
  
  sdat_file <- args[1]
  filter_file <- args[2]
  gtf_file <- args[3]
  out_rds <- args[4]
  out_pdf <- args[5]
  out_txt0 <- args[6]
  out_txt1 <- args[7]
  out_txt2 <- args[8] 
}
batch <- 'Lane'
# source(paste0(rootPath,'/snakemakes/sc_functions.R') )
# source(paste0(rootPath,'/wakecg/CITESeq/CITESeq_functions.R'))
# source(paste0(rootPath,'/snakemakes/Utility_functions.R'))

print('Reading seurat object')
sdat <- readRDS(sdat_file)
DefaultAssay(sdat) <- 'RNA'

print('Reading filter file')
filters <- read.table(filter_file, header = T, sep = ',')
print(filters)
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

### ADD DOWNSAMPLE

##### FEATURE FILTERS
exclude_gene_names <- c()
feat_txt <- filters[which(filters$type == 'feature'),]
feat_txt$N <- NA

### Based on fraction of cells with at least one count
fraction_cells <- filters[which(filters$type == 'feature' & filters$feature == 'fraction_cells'), 'value']
if(length(fraction_cells) > 0){
  print('Feature filter based on fraction of cells with at least one count')
  fraction_cells <- as.numeric(gsub('<', '', fraction_cells))
  ### Original way, using rowSums of boolean table
  # a1 <- Sys.time()
  # genes <- row.names(sdat)[which(rowSums(sdat@assays$RNA@data > 0)/length(colnames(sdat)) <= fraction_cells)]
  # a2 <- Sys.time()
  # orignal_time <- a2-a1
  
   a1 <- Sys.time()
    a <- CreateSeuratObject(sdat@assays$RNA@data, assay = 'RNA', min.cells = length(colnames(sdat)) * fraction_cells)
    genes <- setdiff(row.names(sdat), row.names(a))
   a2 <- Sys.time()
   a2-a1
  
  exclude_gene_names <- c(exclude_gene_names, genes)
  feat_txt[which(feat_txt$feature == 'fraction_cells'), 'N'] <- length(genes)
  hist(rowSums(sdat@assays$RNA@data > 0)/length(row.names(sdat)), xlab = 'Fraction of cells with >0 counts', main = 'N features by fraction of cells')
}

### Based on minimum counts 
count_min <- gsub('<', '', gsub('=', '', filters[which(filters$type == 'feature' & filters$feature == 'counts'), 'value']))
if(length(count_min) > 0){
   print('Feature filter based on the minimum total counts')
   count_min <- as.numeric(count_min) 
   genes <- row.names(sdat)[which(rowSums(sdat@assays$RNA@data) <= count_min)]
   exclude_gene_names <- c(exclude_gene_names, genes)
   feat_txt[which(feat_txt$feature == 'counts'), 'N'] <- length(genes) 
}

### Based on input key words
key_words <- filters[which(filters$type == 'feature' & filters$feature == 'pattern'), 'value']
if(length(key_words) > 0){
  print('Feature filtering based on key-word')
  key_words <- strsplit(key_words, ',')[[1]]
  if(id_type == 'gene_id'){
    all_names <-  gtf$gene_name
    genes <- row.names(sdat)[sapply(row.names(sdat), function(id) any(sapply(key_words, function(y) grepl(pattern = y, x = gtf[id, 'gene_name'], perl = T))))]
  } else if(id_type == 'gene_name'){
    all_names <- row.names(sdat@assays$RNA@counts)
    genes <- all_names[sapply(all_names, function(x) any(sapply(key_words, function(y) grepl(pattern = y, x = x, perl = T))))]
  } else {
    ### Determine whether each is a gene ID or name
    all_names <- sapply(row.names(sdat@assays$RNA@counts), function(x) get_id_name(x, gtf))
    names(all_names) <- row.names(sdat@assays$RNA@counts)
    table(all_names)
    ### For those that match neither, try redoing after splitting on a period.
    is <- which(all_names == 'neither')
    all_names[is] <- sapply(names(all_names)[is], function(x) get_id_name(strsplit(x, '\\.')[[1]][1], gtf))
    ### First gene_names
    gene_names <- names(all_names)[all_names == 'gene_name']
    genes <- gene_names[sapply(gene_names, function(x) any(sapply(key_words, function(y) grepl(pattern = y, x = x, perl = T))))]
    ### Then gene_ids
    gene_ids <- names(all_names)[all_names == 'gene_id']
    genes <- c(genes, gene_ids[sapply(gene_ids, function(id) any(sapply(key_words, function(y) grepl(pattern = y, x = gtf[id, 'gene_name'], perl = T))))])
  }
  #print(all_names)

  exclude_gene_names <- c(exclude_gene_names, genes)
  feat_txt[which(feat_txt$feature == 'pattern'), 'N'] <- length(genes) 
}

exclude_gene_names <- sort(unique(exclude_gene_names))
#gtf[exclude_gene_names, 'gene_name']

feat_txt <- rbind(feat_txt, c('feature', 'total filtered', '', length(unique(exclude_gene_names))))
feat_txt <- rbind(feat_txt, c('feature', 'total remaining', '', length(row.names(sdat)) - length(unique(exclude_gene_names))))

write.table(feat_txt, out_txt0, quote = F, sep =',', row.names = F, col.names = T)
write.table(exclude_gene_names, out_txt1, quote = F, sep = ',', row.names = F, col.names = F)

print('Determing which genes are mitochondiral')
### Determine which genes are mitochondrial by gtf if it is input. Otherwise, by gene name pattern '^MT'.
if(is.data.frame(gtf)){
  sdat[['MT_sum']] <- sum_chromosome(sdat@assays$RNA@counts, chr ='MT', gtf)
} else { 
  mts <- row.names(sdat@assays$RNA)[grepl('^MT-', row.names(sdat@assays$RNA))]
  sdat[['MT_sum']] <- colSums(sdat@assays$RNA@counts[mts, ])
}

##### This should distinguish between results coming from SC_R pipeline (counts) and SmartSeq pipeline (TPM)
### Calculate mitochondial percentage
### TenX
if('nCount_RNA' %in% colnames(sdat@meta.data)){
  print('Proportion MT can be calcualted from nCount_RNA info')
  sdat[['fraction_MT']] <- sdat[['MT_sum']] / sdat[['nCount_RNA']]
} else{
  print('Proportion MT can be calculated from TPM data')
  sdat[['fraction_MT']] <- sdat[['MT_sum']] /colSums(sdat@assays$RNA@counts[, ])
}

### Cell filters
cell_filters <- filters[which(filters$type == 'cell'), 'feature']

### Shannon diversity calculation, if in cell_filters
if(('Shannon_diversity' %in% cell_filters) & (!'Shannon_diversity' %in% colnames(sdat@meta.data))){
 print('Cell filter based on Shannon diversity, calculation is needed') 
 div <- tryCatch({
   diversity(x = sdat@assays$RNA@counts, index = 'shannon', MARGIN = 2)
 }, error = function(e) {
   warning('Cannot calculate Shannon diversity, so removing this filter')
   NA
 }
 )
 if(is.na(div[1])){
   cell_filters <- cell_filters[which(cell_filters != 'Shannon_diversity')]
 } else{
   sdat[['Shannon_diversity']] <- NA
 }
} else if('Shanon_diversity' %in% colnames(sdat@meta.data)){
  print('Cell filter based on Shannon diversity') 
  
  sdat[['Shannon_diversity']] <- as.numeric(sdat@meta.data[, 'Shannon_diversity'])
}

### Determine if any input cell filter criteria aren't in the seurat object
absent <- cell_filters[!cell_filters %in% colnames(sdat@meta.data)]
if(length(absent) > 0){
  print(paste0('Input cell filter criteria is/are not present: ', toString(absent)))
}

print('value')
cell_filters <- cell_filters[cell_filters %in% colnames(sdat@meta.data)]
### Read values from filter file
values <- sapply(cell_filters, function(feature) filters[which(filters$feature == feature), 'value'])
names(values) <- cell_filters
### Set default values
default_values  <- c('<0.03', '1.5(SD),3(SD)', '1.5(SD),3(SD)', '-1,0.2', '3(SD),3(SD)')
names(default_values) <- c('fraction_cells', 'nFeature_RNA', 'nCount_RNA', 'fraction_MT', 'Shannon_diversity')
### Are in sdat
default_values <- default_values[which(names(default_values) %in% colnames(sdat@meta.data))]
### Are not in values already
default_values <- default_values[which(!names(default_values) %in% names(values))]
values <- c(values, default_values)
cell_filters <- c(cell_filters, names(default_values))

print('thresholds')
### Get array of lower,upper bound thresholds from the comma-delimited string. Options are flat values, or if in values of standard deviation, with a trailing (SD)
threshs <- lapply(1:length(values), function(i) threshold_string_seurat(sdat, values[i], names(values)[i]))
names(threshs) <- cell_filters
### Violin plots
print('plots1')

pre_plots <- vector(mode='list', length = length(cell_filters))
names(pre_plots) <- cell_filters
post_plots <- vector(mode='list', length = length(cell_filters))
names(post_plots) <- cell_filters
for(cf in cell_filters){
  sdat[[cf]] <- as.numeric(sdat@meta.data[, cf])
  # print(
  #   VlnPlot(sdat, features = cf, log = T, pt.size = 0) & 
  #     geom_hline(yintercept = threshs[[cf]][1]) & 
  #     geom_hline(yintercept = threshs[[cf]][2])
  #   )
  
  ### Make sure it has levels 
  bats <- unique(sdat@meta.data[, batch])[order(unique(sdat@meta.data[, batch]))]
  sdat@meta.data[, batch] <- factor(sdat@meta.data[, batch], levels = bats)
  ### Violin plots
  batch_ns <- as.character(sapply(levels(sdat@meta.data[, batch]), function(b) length(which(sdat@meta.data[, batch] == b))))
  
  text_y <- max(sdat@meta.data[, cf]) * 0.92
  p <- VlnPlot(sdat, group.by = batch, features = cf, pt.size = 0) + 
    theme(legend.position = "none") +
    ggtitle(paste0(cf, ' by ', batch), subtitle = paste0('Pre filter: ', dim(sdat)[2], ' cells total'))
  p <- p + annotate('text', x = 0.56, y = text_y, size = 3, label = 'N')
  p <- p + annotate('text', x = 1:length(unique(sdat@meta.data[, batch])), y = text_y, label = batch_ns, size = 3, angle = 30)
  #### Only print horizontal line if it filters anything
  ### Lower bound
  if(any(sdat@meta.data[which(!is.na(sdat@meta.data[, cf])), cf] < threshs[[cf]][1])){
    p <- p + geom_hline(yintercept = threshs[[cf]][1])
  }
  ### Upper bound
  if(any(sdat@meta.data[which(!is.na(sdat@meta.data[, cf])), cf] > threshs[[cf]][2])){
    p <- p + geom_hline(yintercept = threshs[[cf]][2])
  }
  pre_plots[[cf]] <- p
}

### Pairs of plots
print('filter')
### Instantiate txt dataframe to hold filter Ns
txt <- as.data.frame(matrix(ncol = 2, nrow = 1 + 2*(length(threshs))))
colnames(txt) <- c('filter', 'n')
txt[, 'filter'] <- c('Original', paste(rep(cell_filters,each = 2), rep(c('low', 'high'), length(cell_filters)), sep = '_'))
row.names(txt) <- txt[, 'filter']

### Subset seurat object for each filter, and record number of cells in txt
sdat_full <- sdat
txt['Original', 'n'] <- dim(sdat_full)[2]
for(i in 1:length(cell_filters)){
  print(cell_filters[i])
  ### Lower bound
  to_keep <- row.names(sdat@meta.data[which(sdat@meta.data[, cell_filters[i]] > threshs[[i]][1]), ])
  sdat <- subset(sdat, cells = to_keep) 
  txt[paste0(cell_filters[i], '_low'), 'n'] <- dim(sdat)[2]
  ### Upper bound
  to_keep <- row.names(sdat@meta.data[which(sdat@meta.data[, cell_filters[i]] < threshs[[i]][2]), ])
  sdat <- subset(sdat, cells = to_keep)
  txt[paste0(cell_filters[i], '_high'), 'n'] <- dim(sdat)[2]
}

### Violin plot for each cell filter
for(cf in cell_filters){
  sdat[[cf]] <- as.numeric(sdat@meta.data[, cf])
  # print(VlnPlot(sdat, features = cf, log = T, pt.size = 0) & 
  #         geom_hline(yintercept = threshs[[cf]][1]) & 
  #         geom_hline(yintercept = threshs[[cf]][2]))

  ### Make sure it has levels 
  bats <- unique(sdat@meta.data[, batch])[order(unique(sdat@meta.data[, batch]))]
  sdat@meta.data[, batch] <- factor(sdat@meta.data[, batch], levels = bats)
  ### Violin plots
  batch_ns <- as.character(sapply(levels(sdat@meta.data[, batch]), function(b) length(which(sdat@meta.data[, batch] == b))))
  
  text_y <- max(sdat@meta.data[, cf]) * 0.92
  p <- VlnPlot(sdat, group.by = batch, features = cf, pt.size = 0) + 
    theme(legend.position = "none") +
    ggtitle(paste0(cf, ' by ', batch), subtitle = paste0('Post filter: ', dim(sdat)[2], ' cells total'))
  p <- p + annotate('text', x = 0.56, y = text_y, size = 3, label = 'N')
  p <- p + annotate('text', x = 1:length(unique(sdat@meta.data[, batch])), y = text_y, label = batch_ns, size = 3, angle = 30)
  #### Only print horizontal line if it filters anything
  ### Lower bound
  if(any(sdat@meta.data[which(!is.na(sdat@meta.data[, cf])), cf] < threshs[[cf]][1])){
    p <- p + geom_hline(yintercept = threshs[[cf]][1])
  }
  ### Upper bound
  if(any(sdat@meta.data[which(!is.na(sdat@meta.data[, cf])), cf] > threshs[[cf]][2])){
    p <- p + geom_hline(yintercept = threshs[[cf]][2])
  }
  
  post_plots[[cf]] <- p
}


pdf(out_pdf)
for(cf in cell_filters){
  print(pre_plots[[cf]])
  print(post_plots[[cf]])
}
dev.off()

# cor(sdat$nFeature_RNA, sdat$fraction_MT)
# plot(sdat$nFeature_RNA, sdat$fraction_MT, xlab = 'N genes', ylab = 'MT proportion')
# #FeatureScatter(sdat, feature1 = "nFeature_RNA", feature2 = "MT_sum", group.by = "sample")
# plot(sdat$nFeature_RNA, sdat$fraction_MT, xlab = 'N genes', ylab = 'MT proportion') +
#   abline(v = NFeat_thresh[1], col = 'red') +
#   abline(v = NFeat_thresh[2], col = 'red') +
#   abline(h = MT_thresh[2], col = 'red')
# plot(sdat$nCount_RNA, sdat$Shannon_diversity, xlab = 'N reads', ylab = 'Shannon Diversity') +
#   abline(v = NCount_thresh[1], col = 'red') +
#   abline(v = NCount_thresh[2], col = 'red') +
#   abline(h = div_thresh[1], col = 'red') +
#   abline(h = div_thresh[2], col = 'red')
# 
# plot(sdat$nFeature_RNA, sdat$fraction_MT, xlab = 'N genes', ylab = 'MT proportion', main = 'post-filtering')
# plot(sdat$nCount_RNA, sdat$Shannon_diversity, xlab = 'N reads', ylab = 'Shannon Diversity', main = 'post-filtering')


#sdat_removed <- subset(sdat_full, cells = colnames(sdat_full)[!colnames(sdat_full) %in% colnames(sdat)])

### Fraction kept by each sample. The first two lines are necessary in the case that an entire sample is filtered
# temp <- table(sdat$sample)
# temp[unique(sdat_full$sample)[which(!unique(sdat_full$sample) %in% unique(sdat$sample))]] <- 0
# temp/table(sdat_full$sample)[names(temp)]

write.table(txt, out_txt2, sep = ',', quote = F, row.names = F, col.names = T)
saveRDS(sdat, file = out_rds)

