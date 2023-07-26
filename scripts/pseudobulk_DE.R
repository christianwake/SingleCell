
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
  qc_name <- 'Run2023-05-14'
  test_name <- 'Time'
  value1 <- '30'
  value2 <- '31'
  strat_name <- 'RNA_clusters'
  strat_values <- '0+1'

  sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Filtered_clustered.RDS')
  
  # project <- '2022620_857.1'
  # qc_name <- '2023-01-09'
  # test <- 'timepoint'
  # strat_name <- 'All'
  # strat_name <- 'seurat_clusters-1+3'
  #sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Filtered.RDS')
  exclude_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Excluded_genes.txt')
  gtf_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/gtf.RDS')
  out_tsv <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/DE/', test, '/', strat_name, '/DE_results.tsv')
  } else{
  args = commandArgs(trailingOnly=TRUE)
  sdat_file <- args[1]
  exclude_file <- args[2]
  gtf_file <- args[3]
  test_name <- args[4]
  value1 <- args[5]
  value2 <- args[6]
  strat_name <- args[7]
  strat_values <- args[8]
  out_tsv <- args[9]
  out_pdf <- args[10]
}
### Read Seurat object
sdat <- readRDS(sdat_file)
### If the gene-exclusion file exists and is not empty, read it
if(file.exists(exclude_file) & file.info(exclude_file)$size != 0){
  exclude_gene_names <- read.table(exclude_file)[,1]
} else{
  exclude_gene_names <- c()
}
### Features to include in the pseudobulk analysis
features <- row.names(sdat)[which(!(row.names(sdat) %in% exclude_gene_names))]

### Read gtf
if(grepl('\\.gtf', gtf_file)){
  gtf <- read_gtf(gtf_file, feature_type = 'gene', atts_of_interest = c('gene_id', 'gene_name',  'gene_biotype'))
} else if(grepl('\\.RDS', gtf_file)){
  gtf <- readRDS(gtf_file)
} else{
  warning('No gtf file')
}
row.names(gtf) <- gtf$gene_id

##### Parse stratification input
### If stratification is input (not 'All', NA or '')
if(strat_name != '' & !is.na(strat_name) & strat_name != 'All'){
  ### "+" separates two groups to be combined
  if(grepl('+', strat_values)){
    strat_values <- strsplit(strat_values, '\\+')[[1]]
  }
  ### If the input stratification actually matches a column name
  if(strat_name %in% colnames(sdat@meta.data)){
    ### If this subset has fewer than 5 samples
    if(sum(table(sdat@meta.data[, strat_name])[strat_values]) < 10){
      stop('Cannot do stratification. Very few cells in subset.')
    } 
    ### Subset sdat to the stratification
    cell_ids<- row.names(sdat@meta.data[which(sdat@meta.data[, strat_name] %in% strat_values), ])
    sdat <- subset(sdat, cells = cell_ids) 
    print('Stratification successful')
  } else{
    stop('Cannot do stratification. Input does not match a covariate name')
  }
}
### Specific to Kristin's T cell data
# if('TCR.status' %in% colnames(sdat@meta.data)){
#   sdat$TCR.status[which(sdat$TCR.status == 'confirmed')] <- 'high confidence'
# }

### Check that the input data is present
if(test_name %in% colnames(sdat@meta.data)){
  ### DE of test_name, pseudo-bulk
  log_file <- gsub('\\.tsv', '\\.log', out_tsv)
  ### Determine which value of the test variable will correspond to larger FC when it has higher expression
  vals <- unique(sdat@meta.data[, test_name])
  ### First check if any are generic names of controls
  control_names <- c('CONTROL', 'C', '0','BASELINE', 'WT', 'WILDTYPE')
  ### Try numeric
  vals_num <- sapply(vals, function(x) gsub('[A-Z]', '', toupper(x)))
  vals_num <- suppressWarnings(as.numeric(vals_num))
  ### Unless the name is a clusters...
  if((!test_name %in% c("RNA_clusters", 'prot_clusters'))){
    if(any(control_names %in% toupper(vals))){
      print('Determining which will be baseline by generic name of test value')
      positive_smaller <- vals[which(toupper(vals) %in% control_names)]
    } else if(all(!is.na(vals_num))){ ### Otherwise, look for a number
      print('Determining which will be baseline by smaller numeric value')
      vals <- vals[order(vals_num)]
      positive_smaller <- vals[1]
    } else { ### Otherwise, don't make a choice
      print('Nothing to determine which will be baseline')
      positive_smaller <- vals[1]
    }
    positive_larger <- vals[which(vals != positive_smaller)]
    
    ### Reorder the factor
    sdat@meta.data[, test_name] <- factor(sdat@meta.data[, test_name], levels = c(positive_smaller, positive_larger))
  }

  Idents(sdat) <- sdat[[test_name]]
  idents <- levels(Idents(sdat))
  ### First, account for "NA" as value1 and value2 if only the test_name is input because it is binary

  if(value1 == 'NA' | value2 == 'NA'){
    value1 <- idents[1]
    value2 <- idents[2]
    if(length(unique(sdat@meta.data[, test_name])) != 2){
      warning(paste('Too many values in the test variable. Can only compare 2. Using ', idents[1], ' and ', idents[2], '.'))
    }
  }
  if(all(value1 %in% vals) & all(value2 %in% vals)){
    de <- FindMarkers(sdat, features = features, ident.1 = value1, ident.2 = value2, min.pct = 0, logfc.threshold = 0, min.cells.group = 0, test.use = 'LR')
    print(paste0('Positive LogFC will be when ', value1, ' > ', value2))
    
    ### Determine whether 'gene_name' or 'gene_id' better matches those in the seurat object
    gtf_cols <- c('gene_name', 'gene_id')
    de$gene_name <- gtf[row.names(de), 'gene_name']
    de <- de[, c('gene_name', colnames(de)[which(colnames(de) != 'gene_name')])]
    
    write.table(de, out_tsv, sep = '\t', quote = F, row.names = T, col.names = T)
    
    id <- row.names(de)[1]
    ptext <- paste0('adjp: ', signif(de[id, 'p_val_adj'], digits = 3))
    fctext <- paste0('FC: ',signif(de[id, 'avg_log2FC'], digits = 3))
    ymax <- max(sdat@assays$RNA@data[id, ]) * 1.09
    ytext1 <- max(sdat@assays$RNA@data[id, ]) * 1.01
    ytext2 <- max(sdat@assays$RNA@data[id, ]) * 1.05
    
    pdf(out_pdf)
    p <- VlnPlot(sdat, features = id, split.by = test_name, y.max = ymax) + ggtitle(id) + ylab('TPM') + 
      annotate("text", x = 1.5, y = ytext1, label = ptext) + 
      annotate("text", x = 1.5, y = ytext2, label = fctext)
    print(p)
    dev.off()
  } else{
    stop('Cannot do test. Input test value(s) are not within the input factor')
  }
} else{
  stop('Cannot do test. Input factor does not match a covariate name')
}


