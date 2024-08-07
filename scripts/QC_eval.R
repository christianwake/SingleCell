library('sys')
library('Seurat')
library('stringr')
library('pheatmap')
library('ggplot2')
#library('umap')
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
  # project <- '2021614_21-002'
  # qc_name <- '2024-01-20'
  # candidates <- 'Sample_Name'
  # tests <- 'Arm,Visit'
  # strats <- 'Cell_subset-B_cells,Cell_subset-Innate'
  # QC_input_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/Test_comparisons.csv')

  project <- '2024605_Hillary_test'
  qc_name <- '2024-06-04'
  candidates <- 'CR_ID'
  QC_input_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/Test_comparisons.csv')
  ### Before filtering
  sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/All_data.RDS')
  plots_path <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/batch_plots/')
  txt_out <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/Batch.txt')
  pdf_out <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/Initial_QC.pdf')
  ### After filtering
  # sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/DSB_normalized_data.RDS')
  # plots_path <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/post_DSB/')
  # txt_out <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Batch.txt')
  # pdf_out <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Post_QC.pdf')
  
 }else{
  args = commandArgs(trailingOnly = TRUE)
  
  sdat_file <- args[1]
  QC_input_file <- args[2]
  candidates <- args[3]  ### CiteSeq snakemake inputs just 1 batch definitively, not candidates
  #tests <- args[3]
  #strats <- args[4]
  plots_path <- args[4] 
  #txt_out1 <- args[5]
  txt_out <- args[5]
  pdf_out <- args[6]
}

qc_input <- read.table(QC_input_file, sep = ',', header = T)
tests <- c(qc_input[, 'test_name'], qc_input[, 'strat_name1'], qc_input[, 'strat_name2'])
tests <- unique(tests[which(!is.na(tests) & tests != '')])

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

downsample_var <- 0.3
if(!interactive()){
  ### Read Seurat data
  sdat <- readRDS(sdat_file)
  print('Done reading Seurat object')
  ### Downsample so I can work interactiveley for testing
  ssub <- DietSeurat(sdat, counts = F, assays = c('RNA', 'prot'))
  cells <- sample(x = colnames(ssub), size = (length(colnames(ssub)) * downsample_var), replace = F)
  
  ssub <- subset(ssub, cells = cells)
  saveRDS(ssub, file = gsub('.RDS', paste0('_DownSampledTo', downsample_var, '.RDS'), sdat_file))
  print('saved downsampled version for testing')
} else{
  sdat <- readRDS(gsub('.RDS', paste0('_DownSampledTo', downsample_var, '.RDS'), sdat_file))
}

print('Batch:')
print(candidates)
if(!all(candidates %in% colnames(sdat@meta.data))){
  warning("An input batch candidate is not in the Seurat object")
}
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
  CR_ID <- sapply(names(nCells_actual), function(sid) unique(sdat@meta.data[which(sdat$Sample_ID == sid), 'CR_ID']))
  ### Put those three info into a data frame
  Ncells <- as.data.frame(cbind(nCells_observed, nCells_actual, Batch, CR_ID))
  Ncells$nCells_observed <- as.numeric(Ncells$nCells_observed)
  Ncells$nCells_actual <- as.numeric(Ncells$nCells_actual)
  ### Add tests and strat vars
  for(cov in tests){
    col <- sapply(names(nCells_actual), function(sid) unique(sdat@meta.data[which(sdat$Sample_ID == sid), cov]))
    if(all(sapply(col, function(x) length(x)) == 1)){
      Ncells[, cov] <- col
    }
  }
  maxvalue <- max(Ncells[, c('nCells_observed', 'nCells_actual')]) * 1.05
  axis_size <- heatmap_N_to_size(length(Batch), size_range = c(10, 1.5), n_range = c(10, 100))
  #### If the color legend will be illegible, change it
  if(length(Batch) > 10){
    color <- colnames(Ncells)[length(colnames(Ncells))]
  } else{
    color <- 'Batch'
  }
  Ncells$Sample_ID <- factor(row.names(Ncells), levels = row.names(Ncells[order(Ncells$nCells_actual),]))
  
  # p1 <- ggplot(Ncells, aes(x = nCells_actual, y = nCells_observed, color = color, label = row.names(Ncells))) + 
  #   geom_text(size = 3) + geom_abline(slope = 1) + ylim(0, maxvalue) + xlim(0, maxvalue)
  p1 <- ggplot(Ncells, aes_string(x = 'nCells_actual', y = 'nCells_observed', color = color, 
                                  label = 'row.names(Ncells)')) + 
    geom_text(size = 3) + geom_abline(slope = 1) + ylim(0, maxvalue) + xlim(0, maxvalue)
  
  #Ncells <- Ncells[order(Ncells$Ncells_actual),]
  p2 <- ggplot(data = Ncells, aes_string(x = 'Sample_ID', y = 'nCells_actual', color = color)) +
    geom_bar(stat = "identity") + coord_flip() +
    theme(axis.text.y = element_text(size = axis_size))
  p3 <- ggplot(data = Ncells, aes_string(x = 'Sample_ID', y = 'nCells_actual', color = 'CR_ID')) +
    geom_bar(stat = "identity") + coord_flip() +
    theme(axis.text.y = element_text(size = axis_size))
  
  ### Key to convert from Sample ID to CR ID
  ### Doesn't work if there are NAs, as the == returns character(0)
  # sids <- unique(sdat$Sample_ID)
  # StoCR <- lapply(sids, function(sid) unique(sdat@meta.data[which(sdat@meta.data$Sample_ID == sid), 'CR_ID']))
  # all(sapply(StoCR, function(x) length(x)) == 1)
  # StoCR <- unlist(StoCR)
  # names(StoCR) <- sids
  
  sdat@meta.data[, 'nCells_RNA']
  ### Confirm that there is only 1 nCells_RNA value per CR_ID
  all(sapply(unique(sdat$CR_ID), function(cr) length(unique(sdat@meta.data[which(sdat@meta.data[, 'CR_ID'] == cr), 'nCells_RNA']))) == 1)
  CellRanger <- sapply(unique(sdat$CR_ID), function(cr) unique(sdat@meta.data[which(sdat@meta.data[, 'CR_ID'] == cr), 'nCells_RNA']))
  Observed <- sapply(unique(Ncells$CR_ID), function(cr) sum(Ncells[which(Ncells$CR_ID == cr), 'nCells_observed']))
  Actual <- sapply(unique(Ncells$CR_ID), function(cr) sum(Ncells[which(Ncells$CR_ID == cr), 'nCells_actual']))
  
  ### Put those three info into a data frame
  Ncells_CR <- as.data.frame(cbind(CellRanger, Observed, Actual))
  ### Add tests and strat vars
  for(cov in tests){
    col <- sapply(names(CellRanger), function(cr) length(unique(sdat@meta.data[which(sdat$CR_ID == cr), cov])))
    if(all(col == 1)){
      Ncells_CR[, cov] <- sapply(names(CellRanger), function(cr) unique(sdat@meta.data[which(sdat$CR_ID == cr), cov]))
    }
  }
  maxvalue <- max(Ncells_CR[, c('Observed', 'Actual', 'CellRanger')]) * 1.05
  p4 <- ggplot(Ncells_CR, aes_string(x = 'Actual', y = 'Observed', label = 'row.names(Ncells_CR)')) + 
    geom_text(size = 3) + geom_abline(slope = 1) + ylim(0, maxvalue) + xlim(0, maxvalue)
  p5 <- ggplot(Ncells_CR, aes_string(x = 'Actual', y = 'CellRanger', label = 'row.names(Ncells_CR)')) + 
    geom_text(size = 3) + geom_abline(slope = 1) + ylim(0, maxvalue) + xlim(0, maxvalue)
  p6 <- ggplot(Ncells_CR, aes_string(x = 'Observed', y = 'CellRanger', label = 'row.names(Ncells_CR)')) + 
    geom_text(size = 3) + geom_abline(slope = 1) + ylim(0, maxvalue) + xlim(0, maxvalue)
  
  #pdf(gsub('.pdf', '_ncells.pdf', pdf_out))
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p5)
  print(p6)
  #dev.off()
}

### string to array
candidates <- trimws(strsplit(candidates, ',')[[1]])
### Try adding these as well, as likely options for some technical effects
candidates <- unique(c(candidates, 'Date_sort', 'Lane'))
### only those in seurat object
candidates <- candidates[candidates %in% colnames(sdat@meta.data)]
### Analysis only done in the 'else' statement, if the user has input some candidate batch variable to evaluate. If they haven't, only output a message recommending that they do.
print(paste0('Batch candidates: ', toString(candidates)))

### Add test and strat variables
candidates <- c(candidates, tests)
candidates <- candidates[which(candidates %in% colnames(sdat@meta.data))]
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
  if(length(covs_numeric) == 0){
    covs_numeric <- c()
  }
  if(length(covs_factor) == 0){
    covs_factor <- c()
  }
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
    if(any(is.na(sdat@meta.data[, batch]))){
      warning(paste0('There are NAs in the batch candiate ', batch))
    }
    ### Make sure it has levels 
    bats <- unique(sdat@meta.data[, batch])[order(unique(sdat@meta.data[, batch]))]
    sdat@meta.data[, batch] <- factor(sdat@meta.data[, batch], levels = bats)
    ##### Violin plots
    ### N cells for annotation
    batch_ns <- as.character(sapply(levels(sdat@meta.data[, batch]), function(b) length(which(sdat@meta.data[, batch] == b))))
    ### Dynamic text size for legibility
    axis_size <- heatmap_N_to_size(length(batch_ns), size_range = c(10,1), n_range = c(10,100))
    annot_size <- heatmap_N_to_size(length(batch_ns), size_range = c(6,1.5), n_range = c(10,100))

    for(comparison in batch_comparisons){
      print(paste0('Violin plot of ', comparison, ' by ', batch))
      ### y-axis location to put extra text
      text_y <- max(sdat@meta.data[, comparison]) * 0.92
      p <- VlnPlot(sdat, group.by = batch, features = comparison, pt.size = 0) + 
        theme(legend.position = "none", axis.text.x = element_text(size = axis_size)) + 
        ggtitle(paste0(comparison, ' by ', batch), subtitle = paste0(signif(adjps[batch, comparison]), ' (', test_names[batch, comparison], ' adjusted p)'))
      p <- p + annotate('text', x = 0.56, y = text_y, size = annot_size, label = 'N')
      p <- p + annotate('text', x = 1:length(unique(sdat@meta.data[, batch])), y = text_y, 
                        label = batch_ns, size = annot_size, angle = 30)
      print(p)
      
      ### Print to separate pdf if too large
      if(length(batch_ns) > 30){
        pdf(gsub('.pdf', paste0('_', batch, '_', comparison, '.pdf'), pdf_out),
            height = heatmap_N_to_pdfsize(length(batch_ns)))
        print(p + coord_flip() + theme(axis.text.x = element_text(size = 8)))
        dev.off()
      }
    }
    #a <- table(sdat$seurat_clusters, sdat@meta.data[, batch])
    ### 857 - Plates which have more than 2x the Cluster 1 average
    #names(plate_match[which(plate_match %in% colnames(a)[which(a['1', ] > 2*mean(a['1', ]))])])
  }
  print('violins complete')
  ### If any of these features are significantly (0.05) associated with one of the features to potentially integrate over...
  a <- rowMeans(adjps[candidates, batch_comparisons])
  
  if(any(candidates %in% row.names(adjps))){
    if(any(adjps[candidates, batch_comparisons] < 0.05)){
      batch_feature <- names(a[which(a == min((a)))])
    } else{
      batch_feature <- 'None'
    }
  } else{
    batch_feature <- 'None'
  }

  ### Print determination to txt_out
  cat(paste0('Automated recommendation for Harmony integration: ', batch_feature), file = txt_out, sep = "\n", append = T)
  cat(paste0('Harmony batch: ', batch_feature), file = txt_out, sep = "\n", append = T)
}
dev.off()
