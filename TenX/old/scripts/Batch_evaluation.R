library('sys')
library('Seurat')
library('stringr')
library('pheatmap')
library('ggplot2')
#library('reticulate'); use_virtualenv("r-reticulate")
library('umap')
library('textshape')
library('dplyr')
library('biomaRt')
library('grid')
library('scales')
library('data.table')

source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/Utility_functions.R')

args = commandArgs(trailingOnly=TRUE)

sdat_file <- args[1]
candidates <- args[2]
method <- args[3]
plots_path <- args[4]
txt_out <- args[5]
pdf_out <- args[6]

# sdat_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021618_galt/results/Filtered.RDS'
# candidates <- c('Lane')
# method <- 'Harmony'
# #gtf_file <- '/hpcdata/vrc/vrc1_data/douek_lab/reference_sets/tenX/Homo_sapiens.GRCh38.93/GRCh38_protein_coding_only/genes/genes.gtf'
# plots_path <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021618_galt/plots/'
# txt_out <- "/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021618_galt/results/batch_evaluation.txt"
# pdf_out <- "/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021618_galt/results/batch_evaluation.pdf"

### Analysis only done in the 'else' statement, if the user has input some candidate batch variable to evaluate. If they haven't, only output a message recommending that they do.
if(candidates == ''){
  cat(paste0('Automated recommendation: Enter batch information in config file', ''), file = txt_out, sep = "\n", append = T)
  cat('Model batch: ', file = txt_out, sep = "\n", append = T)
  cat(paste0('Automated recommendation: Enter batch information in config file', ''), file = pdf_out, sep = "\n")
} else {
  ### Set default method
  if(method == ''){
    method <- 'Harmony'
  }
  sdat <- readRDS(sdat_file)
  candidates <- trimws(strsplit(candidates,',')[[1]])
  
  covs <- sdat@meta.data
  covs$cell_id <- row.names(covs)
  batch_comparisons <- c('nFeature_RNA', 'nCount_RNA', 'MT_sum', 'N_mapped', 'N_multimapped', 'Shannon_diversity')
  batch_comparisons <- batch_comparisons[which(batch_comparisons %in% colnames(covs))]
  ###  Determine type of value and convert
  for(c in colnames(covs)){
    covs[,c] <- determine_type(covs, c)
  }
  covs_numeric <- colnames(covs)[sapply(colnames(covs), function(x) typeof(covs[, x]) == 'double')]
  covs_factor <-  colnames(covs)[sapply(colnames(covs), function(x) typeof(covs[, x]) == 'integer')]

  ### Do pairwise association tests of input covariates
  cov_associations0 <- pairwise_associations(covs, pc = 'p', plots_path = plots_path,
                                             covs_numeric = covs_numeric, 
                                             covs_factor = covs_factor,
                                             covs_binary = c())
  adjps <- cov_associations0[['res']]
  test_names <- cov_associations0[['test']]
  
  DimPlot(sdat, label = T) + NoLegend()
  pdf(pdf_out)
  for(batch in candidates){
    ### Violin plots
    for(comparison in batch_comparisons){
      print(VlnPlot(sdat, group.by = batch, features = comparison) + theme(legend.position = "none") +
              ggtitle(paste0(comparison, ' by ', batch), subtitle = paste0(signif(adjps[batch, comparison]), ' (', test_names[batch, comparison], ' adjusted p)')))
    }
    #a <- table(sdat$seurat_clusters, sdat@meta.data[, batch])
    ### 857 - Plates which have more than 2x the Cluster 1 average
    #names(plate_match[which(plate_match %in% colnames(a)[which(a['1', ] > 2*mean(a['1', ]))])])
  }
  dev.off()
  
  ### If any of these features are significantly (0.05) associated with one of the features to potentially integrate over, do that integration
  a <- rowMeans(adjps[candidates, batch_comparisons])
  
  if(any(adjps[candidates, batch_comparisons] < 0.05)){
    batch_feature <- names(a[which(a == min((a)))])
  } else{
    batch_feature <- 'None'
  }
  cat(paste0('Automated recommendation for Harmony integration: ', batch_feature), file = txt_out, sep = "\n", append = T)
  cat(paste0('Harmony batch: ', batch_feature), file = txt_out, sep = "\n", append = T)
  
}