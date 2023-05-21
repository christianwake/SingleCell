
library('sys')
library('Seurat')
library('SeuratObject')
library('ggplot2')
library('dplyr')
library('viridis')
library('data.table')
library('PKI')
library('tinytex')
library('biomaRt')
library('fgsea')
library('GSEABase')
library('WriteXLS')
library('logistf')
library('readxl')

source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')

if(interactive()){
  project <- '2021600_kristin'
  qc_name <- 'Run2022-11-14'
  test <- 'TCR.status'
  stratification <- 'seurat_clusters-0'
  sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Filtered_clustered.RDS')
  
  # project <- '2022620_857.1'
  # qc_name <- '2023-01-09'
  # test <- 'timepoint'
  # stratification <- 'All'
  # stratification <- 'seurat_clusters-1+3'
  #sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Filtered.RDS')
  exclude_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Excluded_genes.txt')
  gtf_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/gtf.RDS')
  out_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/DE/', test, '/', stratification, '/DE_results.tsv')
  } else{
  args = commandArgs(trailingOnly=TRUE)
  sdat_file <- args[1]
  exclude_file <- args[2]
  gtf_file <- args[3]
  test <- args[4]
  stratification <- args[5]
  out_file <- args[6]
}

sdat <- readRDS(sdat_file)
### If the file exists and is not empty, read it
if(file.exists(exclude_file) & file.info(exclude_file)$size != 0){
  exclude_gene_names <- read.table(exclude_file)[,1]
} else{
  exclude_gene_names <- c()
}
features <- row.names(sdat)[which(!(row.names(sdat) %in% exclude_gene_names))]

##### Parse stratification input
### If stratification is input (not 'All', NA or '')
if(stratification != '' & !is.na(stratification) & stratification != 'All'){
  strat_info <- strsplit(stratification, '-')[[1]]
  strat <- strat_info[2]
  ### "+" separates two groups to be combined
  if(grepl('+', strat)){
    strat <- strsplit(strat, '\\+')[[1]]
  }
  strat_name <- strat_info[1]
  ### If the input stratification actually matches a column name
  if(strat_name %in% colnames(sdat@meta.data)){
    ### If this subset has fewer than 5 samples
    if(sum(table(sdat@meta.data[, strat_name])[strat]) < 10){
      stop('Cannot do stratification. Very few cells in subset.')
    } 
    ### Subset sdat to the stratification
    cell_ids<- row.names(sdat@meta.data[which(sdat@meta.data[, strat_name] %in% strat), ])
    sdat <- subset(sdat, cells = cell_ids) 
    print('Stratification successful')
  } else{
    stop('Cannot do stratification. Input does not match a covariate name')
  }
}
### Specific to Kristin's T cell data
if('TCR.status' %in% colnames(sdat@meta.data)){
  sdat$TCR.status[which(sdat$TCR.status == 'confirmed')] <- 'high confidence'
}

### DE of test, pseudo-bulk
log_file <- gsub('\\.tsv', '\\.log', out_file)
### Determine which value of the test variable will correspond to larger FC when it has higher expression
vals <- unique(sdat@meta.data[, test])
### First check if any are generic names of controls
control_names <- c('CONTROL', 'C', '0','BASELINE', 'WT', 'WILDTYPE')
### Try numeric
vals_num <- sapply(vals, function(x) gsub('[A-Z]', '', toupper(x)))
vals_num <- suppressWarnings(as.numeric(vals_num))

if(any(control_names %in% toupper(vals))){
  print('Determining which will be baseline by generic name of test value')
  positive_smaller <- vals[which(toupper(vals) %in% control_names)]
} else if(all(!is.na(vals_num))){ ### Otherwise, look for a number
  print('Determining which will be baseline by smaller numeric value')
  positive_smaller <- vals[which(vals_num == min(vals_num))]
} else { ### Otherwise, don't make a choice
  print('Nothing to determine which will be baseline')
  positive_smaller <- vals[1]
}
postive_larger <- vals[which(vals != positive_smaller)]

### Reorder
sdat@meta.data[, test] <- factor(sdat@meta.data[, test], levels = c(postive_larger, positive_smaller))

Idents(sdat) <- sdat[[test]]
idents <- levels(Idents(sdat))
if(length(idents) > 2){
  warning(paste('Too many values in the test variable. Can only compare 2. Using ', idents[1], ' and ', idents[2], '.'))
}
de <- FindMarkers(sdat, features = features, ident.1 = idents[1], ident.2 = idents[2], min.pct = 0, logfc.threshold = 0, min.cells.group = 0, test.use = 'LR')
print(paste0('Positive LogFC will be when ', idents[1], ' > ', idents[2]))

### Do direct matching from gtf
if(grepl('\\.gtf', gtf_file)){
  gtf <- read_gtf(gtf_file, feature_type = 'gene', atts_of_interest = c('gene_id', 'gene_name',  'gene_biotype'))
} else if(grepl('\\.RDS', gtf_file)){
  gtf <- readRDS(gtf_file)
} else{
  warning('No gtf file')
}
row.names(gtf) <- gtf$gene_id
### Determine whether 'gene_name' or 'gene_id' better matches those in the seurat object
gtf_cols <- c('gene_name', 'gene_id')
de$gene_name <- gtf[row.names(de), 'gene_name']
de <- de[, c('gene_name', colnames(de)[which(colnames(de) != 'gene_name')])]

write.table(de, out_file, sep = '\t', quote = F, row.names = T, col.names = T)

id <- row.names(de)[1]
ptext <- paste0('adjp: ', signif(de[id, 'p_val_adj'], digits = 3))
fctext <- paste0('FC: ',signif(de[id, 'avg_log2FC'], digits = 3))
ymax <- max(sdat@assays$RNA@data[id, ]) * 1.09
ytext1 <- max(sdat@assays$RNA@data[id, ]) * 1.01
ytext2 <- max(sdat@assays$RNA@data[id, ]) * 1.05

p <- VlnPlot(sdat, features = id, split.by = test, y.max = ymax) + ggtitle(id) + ylab('TPM') + 
  annotate("text", x = 1.5, y = ytext1, label = ptext) + 
  annotate("text", x = 1.5, y = ytext2, label = fctext)
p
