library('sys')
library('Seurat')
library('SeuratWrappers')
library('harmony')
library('stringr')
library('pheatmap')
library('ggplot2')
#library('umap')
#library('textshape')
library('dplyr')
library('biomaRt')
library('grid')
library('scales')
library('vegan')
library('data.table')

source('/data/vrc_his/douek_lab/snakemakes/sc_functions.R')
source('/data/vrc_his/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/Utility_functions.R')

if(interactive()){
  #project <- '2021617_mis-c'
  project <- '2022620_857.1'
  qc_name <- 'Run2022'
  sdat_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/data/All_data.RDS')
  QC_input_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/QC_steps/step2_imputation.csv')
  out_rds <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/PostQC2.RDS')
  out_pdf <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Imputed.pdf')
  #out_txt1 <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/Excluded_genes.txt')
} else{
  args = commandArgs(trailingOnly=TRUE)
  sdat_file <- args[1]
  QC_input_file <- args[2]
  out_rds <- args[3]
  out_pdf <- args[4]
}

qc_input <- read.table(QC_input_file, sep = ',', header = T)
### Currently not doing anything with this because ALRA is the only option
imputation_method <- qc_input$method[1] 

sdat <- readRDS(sdat_file)
DefaultAssay(sdat) <- 'RNA'
sdat[['Shannon_diversity']] <- diversity(x = sdat@assays$RNA@layers$counts, index = 'shannon', MARGIN = 2)
#mtgene <- grep(pattern = "^MT-", rownames(sdat@assays$RNA@layers$counts), value = TRUE)
#sdat$fraction_MT <- Matrix::colSums(sdat@assays$RNA@layers$counts[mtgene, ]) / Matrix::colSums(sdat@assays$RNA@layers$counts)



### ALRA places the data in the 'data' slot of a new assay 'alra'
sdat <- RunALRA(sdat)
### Switch assay names
sdat@assays$RNA_orig <- sdat@assays$RNA
sdat@assays$RNA <- sdat@assays$alra
DefaultAssay(sdat) <- 'RNA'
sdat[['alra']] <- NULL
### Put the imputed counts into 'counts' slot
sdat@assays$RNA@layers$counts <- sdat@assays$RNA@layers$data

### But will need to recalculate the nFeature and nCount, etc.
cols <- c('nCount_RNA', 'nFeature_RNA', 'Shannon_diversity')
# if(length(mtgene) > 0){
#   cols <- c(cols, 'fraction_MT')
# }

cols <- cols[which(cols %in% colnames(sdat@meta.data))]
for(col in cols){
  sdat[[paste0(col, '_orig')]] <- sdat[[col]]
}

### nCount_RNA was removed if data is SmartSeq
if('nCount_RNA' %in% colnames(sdat@meta.data)){
  sdat$nCount_RNA <- colSums(sdat@assays$RNA@layers$counts)
}
sdat$nFeature_RNA <- colSums(sdat@assays$RNA@layers$counts > 0 )
sdat[['Shannon_diversity']] <- diversity(x = sdat@assays$RNA@layers$counts, index = 'shannon', MARGIN = 2)
mtgene <- grep(pattern = "^MT-", rownames(sdat@assays$RNA@layers$counts), value = TRUE)
#sdat$fraction_MT <- Matrix::colSums(sdat@assays$RNA@layers$counts[mtgene, ]) / Matrix::colSums(sdat@assays$RNA@layers$counts)
saveRDS(sdat, file = out_rds)

pdf(out_pdf)
for(col in cols){
  #print(VlnPlot(sdat, features = col) + labs(subtitle = 'After'))
  after <- sdat@meta.data[, col, drop = F]
  before <- sdat@meta.data[, paste0(col, '_orig'), drop = F]
  colnames(before) <- col
  plot_dat <- list(after, before)
  names(plot_dat) <- c('after', 'before')
  plot_dat <- rbindlist(plot_dat, idcol = T)
  colnames(plot_dat) <- c('Integration', col)
  
  p <- ggplot(plot_dat, aes(x = Integration, y = eval(parse(text = col)))) + 
    geom_violin() + ylab(col) 
  print(p)
}

dev.off()

