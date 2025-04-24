library('sys')
library('ggplot2')
library('stringr')
library('GenomicRanges')
library('dplyr')
library('readr')
library('effsize')
library('data.table')
library('biomaRt')
library('fgsea')
library('GSEABase')
library('readxl')
library('Seurat')

source('/data/vrc_his/douek_lab/snakemakes/Utility_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/DE_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/sc_functions.R')

if(interactive()){
  project <- '2021600_kristin'
  qc_name <- 'Run2022-11-14'
  test <- 'TCR.status'
  strats <- c('seurat_clusters-0', 'seurat_clusters-3')
  species <- 'hsapiens'
  sdat_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Filtered_clustered.RDS')

  de_files <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/DE/', test, '/', strats, '/DE_results.tsv')
  exclude_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Excluded_genes.txt')
  gtf_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/data/gtf.RDS')
  #gmt_file <- '/data/vrc_his/douek_lab/wakecg/genesets/c2.cp.v7.2.symbols.gmt'
  #gmt_file <- '/data/vrc_his/douek_lab/wakecg/genesets/h.all.v2022.1.Hs.symbols.gmt'
  gmt_file <- '/data/vrc_his/douek_lab/wakecg/genesets/c7.all.v2022.1.Hs.symbols.gmt'
  out_rds <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Classified.RDS')
} else{
  args = commandArgs(trailingOnly=TRUE)
  sdat_file <- args[1]
  res_file <- args[2]
  exclude_file <- args[3]
  gtf_file <- args[4]
  gmt_file <- args[5]
  species <- args[6]
  out_tsv <- args[7]
  #out_pdf <- args[8]
  custom_sets <- args[8]
}

# ### If the file exists and is not empty, read it
# if(file.exists(exclude_file) & file.info(exclude_file)$size != 0){
#   exclude_gene_names <- read.table(exclude_file)[,1]
# } else{
#   exclude_gene_names <- c()
# }

### Read counts
sdat <- readRDS(sdat_file)
### Specific to Kristin's T cell data
if('TCR.status' %in% colnames(sdat@meta.data)){
  sdat$TCR.status[which(sdat$TCR.status == 'confirmed')] <- 'high confidence'
}
#Idents(sdat) <- sdat[[test]]
Idents(sdat) <- test
idents <- levels(Idents(sdat))

# norm_counts <- as.data.frame(sdat@assays$RNA@data)
# norm_counts <- norm_counts[which(!(row.names(norm_counts) %in% exclude_gene_names)),]

### Convert to numeric
#norm_counts <- mutate_all(norm_counts, function(x) as.numeric(x))

### Do direct matching from gtf
if(grepl('\\.gtf', gtf_file)){
  gtf <- read_gtf(gtf_file, feature_type = 'gene', atts_of_interest = c('gene_id', 'gene_name',  'gene_biotype'))
} else if(grepl('\\.RDS', gtf_file)){
  gtf <- readRDS(gtf_file)
} else{
  warning('No gtf file')
}
row.names(gtf) <- gtf$gene_id

### Read DE results
des <- lapply(de_files, function(de_file) read.table(de_file, header = T, check.names = F, stringsAsFactors = F, na.strings = c("", "NA"), sep = '\t'))
names(des) <- strats

### Run biomaRt to get gene name information from ensembl IDs.
res$ensembl <- row.names(res)
#colnames(res) <- gsub('gene_name', 'gene_name_gtf', colnames(res))
if(species != 'hsapiens'){
  ### BiomaRt, for gene names
  id_key <- run_biomaRt(res, to_merge = F, species = species)
  ### (For now) Only keep one-to-one orthologs.
  id_key <- id_key[which(id_key$ortholog_type == 'ortholog_one2one'),]
  id_key[id_key == ''] <- NA
  row.names(id_key) <- id_key$ens_short
  colSums(is.na(id_key))
  ### Prioritize Gene names from biomart human orthologs
  res[row.names(id_key), 'gene_name'] <- id_key$gene_name
  ### But when that is NA, use gene name from the other species' gtf
  res[which(is.na(res$gene_name) & !is.na(res$gene_name_gtf)), 'gene_name'] <- res[which(is.na(res$gene_name) & !is.na(res$gene_name_gtf)), 'gene_name_gtf']
  res[, 'gene_name_plot'] <- res$gene_name
  res[which(is.na(res$gene_name)), 'gene_name_plot'] <- row.names(res[which(is.na(res$gene_name)), ])
  ### Adjust it if there are duplicates in 'gene_name_plot'
  dups <- names(table(res$gene_name_plot))[which(table(res$gene_name_plot) > 1)]
  res[which(res$gene_name_plot %in% dups), 'gene_name_plot'] <- 
    paste0(res[which(res$gene_name_plot %in% dups), 'gene_name_plot'], '(', row.names(res[which(res$gene_name_plot %in% dups), ]), ')')
} else{
  id_key <- gtf[, c('gene_id', 'gene_name')]
  colnames(id_key) <- c('gene_id', 'hsapiens_name')
  if('gene_name_gtf' %in% colnames(res)){
    res$gene_name <- res$gene_name_gtf
  } else{
    res$gene_name <- id_key[row.names(res), 'hsapiens_name']
    res$gene_name_plot <- res$gene_name
  }
}

TopN <- 1
sdat@meta.data[, paste0(test, '_classification')] <- NA
for(stratification in strats){
  ##### Parse stratification input
  strat_info <- strsplit(stratification, '-')[[1]]
  strat <- strat_info[2]
  ### "+" separates two groups to be combined
  if(grepl('+', strat)){
    strat <- strsplit(strat, '\\+')[[1]]
  }
  strat_name <- strat_info[1]
  ### If the input stratification actually matches a column name
  if(!(strat_name %in% colnames(sdat@meta.data))){
    warning('Stratification is not in the Seurat object. Was the wrong object input?')
  }
  ### Cell IDs from those in this stratification 
  cell_ids <- row.names(sdat@meta.data[which(sdat@meta.data[, strat_name] %in% strat), ])
  ### 
  de <- des[[stratification]]
  sig <- de[which(de$p_val_adj < 0.1),]
  if(TopN > length(row.names(sig))){
    topN <- length(row.names(sig))
  } else{
    topN <- TopN
  }
  sig <- sig[order(abs(sig$avg_log2FC), decreasing = T)[1:topN], ]
  nearer_dat <- as.data.frame(matrix(ncol = topN, nrow = length(cell_ids)))
  row.names(nearer_dat) <- cell_ids
  colnames(nearer_dat) <- row.names(sig)
  for(rn in row.names(sig)){
    avgs <- sapply(idents, function(id) mean(sdat@assays$RNA[rn, which(sdat[[strat_name]] == strat & sdat[[test]] == id)]))
    names(avgs) <- idents
    
    #sdat@assays$RNA[row.names(sig), cell_id]
    ### For each cell, select the ident whose mean of gene X is closest to that cells gene X value
    nearer_dat[cell_ids, rn] <- sapply(cell_ids, function(cell_id) names(avgs)[order(abs(avgs - as.vector(sdat@assays$RNA[rn, cell_id])), decreasing = F)[1]])
  }
  ### Bind all ident counts together
  scores <- rbind(sapply(idents, function(val) rowSums(nearer_dat == val)))
  ### Call each cell whichever ident has the highest count
  sdat@meta.data[cell_ids, paste0(test, '_classification')] <- sapply(cell_ids, function(cell_id) colnames(scores)[order(scores[cell_id, ], decreasing = T)[1]])
}

saveRDS(sdat, file = out_rds)

acc_tab <- table(sdat[[test]][[1]], sdat[[paste0(test, '_classification')]][[1]])
acc_tab
(acc_tab[1,2] + acc_tab[2,1]) / sum(acc_tab)

pdf('/data/vrc_his/douek_lab/projects/RNASeq/2021600_kristin/specificity/Classification_by_top1gene.pdf')
DimPlot(sdat, group.by = test, label = T, reduction = 'umap')
DimPlot(sdat, group.by = paste0(test, '_classification'), label = T, reduction = 'umap')
dev.off()

c0 <- sdat@meta.data[which(sdat$seurat_clusters == '0'),]
table(c0[,test], c0[,paste0(test, '_classification')])

c3 <- sdat@meta.data[which(sdat$seurat_clusters == '0'),]
table(c3[,test], c3[,paste0(test, '_classification')])
