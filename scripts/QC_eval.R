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

source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R')

if(interactive()){
  project <- '2021614_21-002'
  #project <- '2022619_857.3b'
  qc_name <- 'Both_celltypes'
  candidates <- 'Date_sort'
  ### Before filtering
  sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/All_data.RDS')
  plots_path <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/batch_plots/')
  txt_out <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/Batch.txt')
  pdf_out <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/N_cells.pdf')
  ### After filtering
  sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/DSB_normalized_data.RDS')
  plots_path <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/post_DSB/')
  txt_out <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Batch.txt')
  pdf_out <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/N_cells.pdf')
  
  
  # sdat_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021618_galt/data/All_data.RDS'
  # QC_input_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021618_galt/QC_steps/QC_steps.csv'
  # candidates <- 'Lane'
  # gtf_file <- '/hpcdata/vrc/vrc1_data/douek_lab/reference_sets/tenX/Homo_sapiens.GRCh38.93/GRCh38_protein_coding_only/genes/genes.gtf'
  # plots_path <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021618_galt/plots/'
  # txt_out2 <- "/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021618_galt/results/Go1/batch_evaluation.txt"
  # pdf_out <- "/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021618_galt/results/Go1/batch_evaluation.pdf"
  
  # sdat_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021617_mis-c/data/All_data.RDS'
  # QC_input_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021617_mis-c/QC_steps/QC_steps.csv'
  # candidates <- 'Lane,CR_ID'
  # gtf_file <- '/hpcdata/vrc/vrc1_data/douek_lab/reference_sets/tenX/Homo_sapiens.GRCh38.93/GRCh38_protein_coding_only/genes/genes.gtf'
  # plots_path <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021617_mis-c/plots/'
  # txt_out1 <- "/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021617_mis-c/results/QC_steps.csv"
  # txt_out <- "/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021617_mis-c/results/batch_evaluation.txt"
  # pdf_out <- "/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021617_mis-c/results/batch_evaluation.pdf"
  
  #sdat_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021600_kristin/data/All_data.RDS'
  #QC_input_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021600_kristin/QC_steps.csv'
  #candidates <- 'Lane'
  #gtf_file <- '/hpcdata/vrc/vrc1_data/douek_lab/reference_sets/tenX/Homo_sapiens.GRCh38.93/GRCh38_protein_coding_only/genes/genes.gtf'
  #plots_path <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021600_kristin/plots/'
  #txt_out <- "/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021600_kristin/results/Run2022/batch_evaluation.txt"
  #pdf_out <- "/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021600_kristin/results/Run2022/batch_evaluation.pdf"

  #sdat_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022619_857.3b/data/All_data.RDS'
  #QC_input_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022619_857.3b/QC_steps.csv'
  #candidates <- 'Lane'
  #gtf_file <- '/hpcdata/vrc/vrc1_data/douek_lab/reference_sets/tenX/Mmul_10_proteincoding/genes/genes.gtf'
  #plots_path <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022619_857.3b/plots/'
  #txt_out <- "/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022619_857.3b/results/Run2023/batch_evaluation.txt"
  #pdf_out <- "/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022619_857.3b/results/Run2023/batch_evaluation.pdf"
}else{
  args = commandArgs(trailingOnly = TRUE)
  
  sdat_file <- args[1]
  #QC_input_file <- args[2]
  candidates <- args[2]  ### CiteSeq snakemake inputs just 1 batch definiteively, not candidates
  plots_path <- args[3] 
  #txt_out1 <- args[5]
  txt_out <- args[4]
  pdf_out <- args[5]
}

# qc_input <- read.table(QC_input_file, sep = ',', header = T)
# qc_output <- as.data.frame(matrix(nrow = 4, ncol = 3))
# colnames(qc_output) <- c('step_num', 'step_name', 'skip')
# qc_output$step_name <- c('Sample filter', 'Imputation', 'Cell filter', 'Integration')
# qc_output$step_num <- paste0('step', 1:4)
# qc_output[, 'skip'] <- T
# qc_output[qc_input[which(!is.na(qc_input[, 'file'])), 'step'], 'skip'] <- F
# ### Add step 0 (not QC)
# qc_output['0', ] <- c('step0', 'Base_data', F)
# qc_output <- qc_output[as.character(0:4), ]
# write.table(qc_output, txt_out1, quote =F, sep = ',', row.names = T, col.names = T)
sdat <- readRDS(sdat_file)

print('Batch:')
print(candidates)
print(table(sdat@meta.data[, candidates]))

pdf(pdf_out)
if('Cell_Count' %in% colnames(sdat@meta.data)){
  print(sum(is.na(sdat$Cell_Count)))
  print('Making Cell_Count plots')
  sdat$Cell_Count <- as.numeric(sdat$Cell_Count)
  ### Confirm that Cell_count is uniform within samples
  if(all(sapply(unique(sdat$Sample_ID), function(sid) length(unique(sdat@meta.data[which(sdat$Sample_ID == sid), 'Cell_Count'])) == 1))){
    print('Cell_Count is uniform within Sample_ID')
  }
  ### Get the Cell_Count (per Sample_ID)
  nCells_actual <- sapply(unique(sdat$Sample_ID), function(sid) unique(sdat@meta.data[which(sdat$Sample_ID == sid), 'Cell_Count']))
  ### Get the Seurat objects number of cells per Sample_ID
  nCells_observed <- table(sdat$Sample_ID)
  nCells_observed <- nCells_observed[names(nCells_actual)]
  ### Get the batch value for each Sample_ID
  Batch <- sapply(names(nCells_actual), function(sid) unique(sdat@meta.data[which(sdat$Sample_ID == sid), candidates[1]]))
  ### Put those three info into a data frame
  Ncells <- as.data.frame(cbind(nCells_observed, nCells_actual, Batch))
  Ncells$nCells_observed <- as.numeric(Ncells$nCells_observed)
  Ncells$nCells_actual <- as.numeric(Ncells$nCells_actual)
  
  
  maxvalue <- max(Ncells[, c('nCells_observed', 'nCells_actual')]) * 1.05
  p1 <- ggplot(Ncells, aes(x = nCells_actual, y = nCells_observed, color = Batch, label = row.names(Ncells))) + 
    geom_text(size = 3) + geom_abline(slope = 1) + ylim(0, maxvalue) + xlim(0, maxvalue)
  
  Ncells$Sample_ID <- factor(row.names(Ncells), levels = row.names(Ncells[order(Ncells$nCells_actual),]))
  #Ncells <- Ncells[order(Ncells$Ncells_actual),]
  p2 <- ggplot(data = Ncells, aes(x = Sample_ID, y = nCells_actual, color = Batch)) +
    geom_bar(stat = "identity") + coord_flip()
  
  
  #pdf(gsub('.pdf', '_ncells.pdf', pdf_out))
  print(p1)
  print(p2)
  #dev.off()
}


### string to array
candidates <- trimws(strsplit(candidates, ',')[[1]])
### only those in seurat object
candidates <- candidates[candidates %in% colnames(sdat@meta.data)]
### Analysis only done in the 'else' statement, if the user has input some candidate batch variable to evaluate. If they haven't, only output a message recommending that they do.
if(length(candidates) == 0){
  cat(paste0('Automated recommendation: Enter batch information in config file', ''), file = txt_out, sep = "\n", append = T)
  cat('Model batch: ', file = txt_out, sep = "\n", append = T)
  cat(paste0('Automated recommendation: Enter batch information in config file', ''), file = pdf_out, sep = "\n")
} else {
  covs <- sdat@meta.data
  covs$cell_id <- row.names(covs)
  ### Things to compare to batch?
  batch_comparisons <- c('nFeature_RNA', 'nCount_RNA', 'MT_sum', 'MT_percent', 'fraction_MT', 'N_mapped', 'N_multimapped', 'Shannon_diversity', 'diversity',
                         'nFeature_prot', 'nCount_prot', 'nFeature_dsb', 'nCount_dsb')
  batch_comparisons <- batch_comparisons[which(batch_comparisons %in% colnames(covs))]
  ### Reduce covs to those we are using
  covs <- covs[, c(candidates, batch_comparisons)]
  ###  Determine type of value and convert
  for(c in colnames(covs)){
    covs[,c] <- determine_type(covs, c)
    sdat@meta.data[, c] <- determine_type(covs, c)
    print(paste0(c, ' is ', typeof(covs[, c])))
    
  }
  covs_numeric <- colnames(covs)[sapply(colnames(covs), function(x) typeof(covs[, x]) == 'double')]
  covs_factor <-  colnames(covs)[sapply(colnames(covs), function(x) typeof(covs[, x]) == 'integer')]

  ### Do pairwise association tests of input covariates
  cov_associations0 <- pairwise_associations(covs, pc = 'p', plots_path = plots_path,
                                             covs_numeric = covs_numeric, 
                                             covs_factor = covs_factor,
                                             covs_binary = c(),
                                             verbose = T)
  adjps <- cov_associations0[['res']]
  test_names <- cov_associations0[['test']]
  
  ### If there are dimentionsality reductions, lets look at one
  if(length(names(sdat@reductions)) > 0){
    DimPlot(sdat, label = T) + NoLegend()
  }
  for(batch in candidates){
    ### Make sure it has levels 
    bats <- unique(sdat@meta.data[, batch])[order(unique(sdat@meta.data[, batch]))]
    sdat@meta.data[, batch] <- factor(sdat@meta.data[, batch], levels = bats)
    ### Violin plots
    batch_ns <- as.character(sapply(levels(sdat@meta.data[, batch]), function(b) length(which(sdat@meta.data[, batch] == b))))
    for(comparison in batch_comparisons){
      print(comparison)
      text_y <- max(sdat@meta.data[, comparison]) * 0.92
      p <- VlnPlot(sdat, group.by = batch, features = comparison, pt.size = 0) + theme(legend.position = "none") +
              ggtitle(paste0(comparison, ' by ', batch), subtitle = paste0(signif(adjps[batch, comparison]), ' (', test_names[batch, comparison], ' adjusted p)'))
      p <- p + annotate('text', x = 0.56, y = text_y, size = 3, label = 'N')
      p <- p + annotate('text', x = 1:length(unique(sdat@meta.data[, batch])), y = text_y, label = batch_ns, size = 3, angle = 30)
      print(p)
    }
    #a <- table(sdat$seurat_clusters, sdat@meta.data[, batch])
    ### 857 - Plates which have more than 2x the Cluster 1 average
    #names(plate_match[which(plate_match %in% colnames(a)[which(a['1', ] > 2*mean(a['1', ]))])])
  }
  print('violins complete')
  ### If any of these features are significantly (0.05) associated with one of the features to potentially integrate over...
  a <- rowMeans(adjps[candidates, batch_comparisons])
  
  if(any(adjps[candidates, batch_comparisons] < 0.05)){
    batch_feature <- names(a[which(a == min((a)))])
  } else{
    batch_feature <- 'None'
  }
  ### Print determination to txt_out
  cat(paste0('Automated recommendation for Harmony integration: ', batch_feature), file = txt_out, sep = "\n", append = T)
  cat(paste0('Harmony batch: ', batch_feature), file = txt_out, sep = "\n", append = T)
}
dev.off()
