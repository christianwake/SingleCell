
my_vln_plot <- function(sdat, de, gene_names, genes = NA, groupby = 'seurat_clusters', lone_cluster = NA, log = F){
  ### Allow for 'gene_names' input to be name or ID, and any length
  ### If the input is not in de$gene, assume it is the gene name, so get its ID
  if('ens_short' %in% colnames(genes)){
    ids <- sapply(gene_names, function(x) ifelse(x %in% de$gene, x, genes[which(genes$gene_name == x), 'ens_short']))  
  } else{
    ids <- gene_names
  }
  
  subtitle <- ''
  ### If DE was done on a lone cluster use that info in subtitles
  if(!is.na(lone_cluster)){
    if(ids %in% row.names(de)){
      subtitle <- paste0('DE ', lone_cluster, ' v. others adjp: ', signif(de[ids, 'p_val_adj']), ', logFC: ', signif(de[ids, 'avg_log2FC']))
    } 
  }
  if(log){
    vplot <- VlnPlot(sdat, group.by = groupby, features = ids) + 
      ggtitle(gene_names, subtitle = subtitle) +  theme(legend.position = "none") + scale_y_log10()
  } else{
    vplot <- VlnPlot(sdat, group.by = groupby, features = ids[[1]]) + 
      ggtitle(gene_names, subtitle = subtitle) +  theme(legend.position = "none")
  }
  return(vplot)
}


marker_cell_filter <- function(cseq, output = 'Seurat', negatives = c('CD3'), positives = c('CD19')){
  ### Reduce to input markers actually matching those in the seurat object
  negatives <- negatives[which(negatives %in% row.names(cseq))]
  positives <- positives[which(positives %in% row.names(cseq))]
  if(length(c(negatives, positives)) == 0){
    print('No valid markers input')
  }
  ncells1 <- dim(cseq)[2]
  
  ### Filter cells based on counts of the negative markers
  neg_remove <- c()
  for(marker in negatives){
    to_remove <- colnames(cseq)[which(cseq@assays$prot@counts[marker, ] > 0)]
    print(paste0(length(to_remove), ' cells with > 0 counts of ', marker))
    neg_remove <- c(neg_remove, to_remove)
  }
  
  ### Filter cells based on absence of counts of the positive markers
  pos_remove <- c()
  for(marker in positives){
    to_remove <- colnames(cseq)[which(cseq@assays$prot@counts[marker, ] == 0)]
    print(paste0(length(to_remove), ' cells with 0 counts of ', marker))
    pos_remove <- c(pos_remove, to_remove)
  }
  to_remove <- unique(c(neg_remove, pos_remove))
  out <- to_remove
  if(output == 'Seurat'){
    ### Filter
    to_keep <- colnames(cseq)[which(!colnames(cseq) %in% to_remove)]
    cseq <- subset(cseq, cells = to_keep)
    print(paste0('N cells removed by marker in total: ', (ncells1 - dim(cseq)[2]), ' of ', ncells1))
    out <- cseq
  }
  return(out)
}


threshold_string_seurat <- function(cseq, string, var, positive_integer = F){
  threshs <- c()
  mt <- strsplit(string, ',')[[1]]
  sd <- grepl('(SD)', mt)
  ### Lower bound
  if(sd[[1]] == F){
    threshs <- c(threshs, as.numeric(mt[1]))
  } else{
    threshs <- c(threshs, mean(cseq[[var]][,1]) - as.numeric(gsub('\\(SD\\)', '', mt[1])) * sd(cseq[[var]][,1]))
  }
  if(positive_integer & threshs[1] < 0){
    threshs[1] <- 0
  }
  ### Upperbound
  if(sd[[2]] == F){
    threshs <- c(threshs, as.numeric(mt[2]))
  } else{
    threshs <- c(threshs, mean(cseq[[var]][,1]) + as.numeric(gsub('\\(SD\\)', '', mt[2])) * sd(cseq[[var]][,1]))
  }
  return(as.numeric(threshs))
}

### Returns array of cell_ids to keep
cell_filter <- function(cseq, MT_sum = NA, nFeature_RNA = NA, out_pdf = NA, out_txt = NA, out_xls = NA){
  txt <- as.data.frame(matrix(ncol = 2, nrow = 5))
  colnames(txt) <- c('filter', 'n')
  txt[, 'filter'] <- c('nFeature_low', 'nFeature_high', 'MT_low', 'MT_high', 'remaining')
  row.names(nums) <- nums$Assignment
  nums[, 'N_post_lower_nFeature'] <- NA
  nums[, 'N_post_upper_nFeature'] <- NA
  nums[, 'N_post_upper_MTsum'] <- NA
  
  cseq_full <- cseq
  if(!is.na(nFeature_RNA)){
    print('filtering nFeature_RNA')
    ### Get array of lower,upper bound thresholds from the comma-delimited string. Options are flat values, or if in values of standard deviation, with a trailing (SD)
    NFeat_thresh <- threshold_string_seurat(cseq, nFeature_RNA, 'nFeature_RNA')
    cseq <- subset(cseq, subset = nFeature_RNA >  NFeat_thresh[1])
    txt[1, 'n'] <- (dim(cseq_full)[2] - dim(cseq)[2])
    nums[names(table(cseq$Assignment)), 'N_post_lower_nFeature'] <- table(cseq$Assignment)
    cseq <- subset(cseq, subset = nFeature_RNA <  NFeat_thresh[2])
    txt[2, 'n'] <- (dim(cseq_full)[2] - txt[1, 'n'] - dim(cseq)[2])
    nums[names(table(cseq$Assignment)), 'N_post_upper_nFeature'] <- table(cseq$Assignment)
  }
  if(!is.na(MT_sum)){
    print('filtering MT_sum')
    ### Get array of lower,upper bound thresholds from the comma-delimited string. Options are flat values, or if in values of standard deviation, with a trailing (SD)
    MT_thresh <- threshold_string_seurat(cseq, MT_sum, 'MT_sum')
    cseq <- subset(cseq, subset = MT_sum >  MT_thresh[1])
    txt[3, 'n'] <- (dim(cseq_full)[2] - txt[1, 'n'] - txt[2, 'n'] - dim(cseq)[2])
    cseq <- subset(cseq, subset = MT_sum < MT_thresh[2])
    txt[4, 'n'] <- (dim(cseq_full)[2] - txt[1, 'n'] - txt[2, 'n'] - txt[3, 'n']  - dim(cseq)[2])
    nums[names(table(cseq$Assignment)), 'N_post_upper_MTsum'] <- table(cseq$Assignment)
    txt[which(txt$filter == 'remaining'), 'n'] <- dim(cseq)[2]
    
  }
  
  if(!is.na(out_pdf)){
    pdf(out_pdf)
    VlnPlot(cseq_full, features = "nFeature_RNA")
    VlnPlot(cseq_full, features = 'MT_sum')
    
    cor(cseq_full$nFeature_RNA, cseq_full$MT_sum)
    plot(cseq_full$nFeature_RNA, cseq_full$MT_sum, xlab = 'N genes', ylab = 'Sum MT TPM')
    
    plot(cseq_full$nFeature_RNA, cseq_full$MT_sum, xlab = 'N genes', ylab = 'Sum MT TPM') +
      abline(v = NFeat_thresh[1], col = 'red') +
      abline(v = NFeat_thresh[2], col = 'red') +
      abline(h = MT_thresh[2], col = 'red')
    VlnPlot(cseq_full, group.by = 'Assignment', features = 'nFeature_RNA')
    VlnPlot(cseq_full, group.by = 'Assignment', features = 'nCount_RNA')
    
    plot(cseq$nFeature_RNA, cseq$MT_sum, xlab = 'N genes', ylab = 'Sum MT TPM')
    
    dev.off()
  }
  
  if(!is.na(out_txt)){
    write.table(txt, out_txt, sep = ',', quote = F, row.names =F, col.names = T)
  }
  if(!is.na(out_xls)){
    WriteXLS(c('txt', 'nums'), ExcelFileName = out_xls, SheetNames = c('Summary', 'Sample'))
  }
  return(colnames(cseq))
}
### Run FindAllMarkers() for differential expression of data type 'assay' (e.g. 'RNA') of each cluster of 'cluster_type' (e.g. 'RNA_clusters') versus others.
### If features is not input, default of FindAllMarkers is to use all features.
run_fam <- function(cseq, assay, cluster_type, topn, features = NULL, pdf_file = NA){
  Idents(cseq) <- cseq[[cluster_type]]
  DefaultAssay(cseq) <- assay
  
  diffs <- FindAllMarkers(cseq, logfc.threshold = 0, return.thresh = 1, assay = assay, features = features)
  print('FAM complete')
  diffs$diff <- diffs$pct.2 - diffs$pct.1
  
  if(!is.na(pdf_file)){
    pdf(pdf_file)
    topn <- 5
    ### By p-value
    #top_genes <- unique(unlist(lapply(unique(diffs$cluster), function(c) diffs[which(diffs$cluster == c), ][order(diffs[which(diffs$cluster == c), ][, 'p_val'])[1:topn], 'gene'])))
    ### By difference
    top_genes <- unique(unlist(lapply(unique(diffs$cluster), function(c) diffs[which(diffs$cluster == c), ][order(diffs[which(diffs$cluster == c), ][, 'diff'], decreasing = T)[1:topn], 'gene'])))
    ### The default scaling doesn't include all features, to save time and space
    #cseq <- FindVariableFeatures(cseq, assay = 'RNA', selection.method = "vst")
    #cseq <- ScaleData(cseq, assay = "RNA", features = row.names(cseq))
    
    DoHeatmap(cseq, features = top_genes) + ggtitle(paste0(assay, ' expression of top ', topn , ' DE features per ', cluster_type, ' cluster'))
    dev.off()
  }
  return(diffs)
}

### The following are for cluster filters
too_few_cells <- function(cseq, cluster_name, min_cells){
  ### Filter any largely empty clusters
  cluster_sizes <- table(cseq[[cluster_name]])
  return(cluster_sizes[which(cluster_sizes < as.numeric(min_cells))])
}

too_few_reads <- function(cseq, datatype, cluster_list, min_expression){
  nCount <- paste0('nCount_', datatype)
  cluster_means <- sapply(names(cluster_list), function(cluster) mean(cseq[[nCount]][cluster_list[[cluster]],]))
  return(cluster_means[which(cluster_means < as.numeric(min_expression))])
}

too_few_reads_relative <- function(cseq, datatype, cluster_list, min_expression){
  nCount <- paste0('nCount_', datatype)
  cluster_means <- sapply(names(cluster_list), function(cluster) mean(cseq[[nCount]][cluster_list[[cluster]],]))
  sans_cluster_means <- sapply(names(cluster_list), function(cluster) mean(cseq[[nCount]][unique(unlist(cluster_list[names(cluster_list) != cluster])),])) 
  ratio <- cluster_means/sans_cluster_means
  return(ratio[which(ratio < as.numeric(min_expression))])
}

too_many_negatives <- function(cseq, datatype, cluster_list, negative_markers, max_expression){
  ### Split on commas
  negative_markers <- trimws(strsplit(negative_markers, ',')[[1]])
  ### Only keep those actually in the data
  negative_markers <- negative_markers[which(negative_markers %in% row.names(cseq@assays[[datatype]]@counts))]
  
  clusters <- c()
  ### Loop over input features
  for(negative in negative_markers){
    cluster_means <- sapply(cluster_list, function(cluster) mean(as.numeric(cseq@assays[[datatype]][negative, cluster])))
    clusters <- c(clusters, cluster_means[which(cluster_means > max_expression)])
  }
  return(clusters)
}

too_few_positives <- function(cseq, datatype, cluster_list, positive_markers, min_expression){
  ### Split on commas
  positive_markers <- trimws(strsplit(positive_markers, ',')[[1]])
  ### Only keep those actually in the data
  positive_markers <- positive_markers[which(positive_markers %in% row.names(cseq@assays[[datatype]]@counts))]
  
  clusters <- c()
  ### Loop over input features
  for(negative in positive_markers){
    cluster_means <- sapply(cluster_list, function(cluster) mean(as.numeric(cseq@assays[[datatype]][negative, cluster])))
    clusters <- c(clusters, cluster_means[which(cluster_means < min_expression)])
  }
  clusters <- unique(names(clusters))
  return(clusters)
}

filter_clusters <- function(cseq, datatype, filter_type, features, value, cluster_name = NA){
  if(is.na(cluster_name)){
    cluster_name <- paste0(datatype, '_clusters')
  }
  nCount <- paste0('nCount_', datatype)
  filter_type <- gsub(' ', '', gsub('_', '', tolower(filter_type)))
  ### Get the cell_ids  for each cluster
  cluster_names <- levels(cseq[[cluster_name]][,1])
  ### Cell membership of each cluster as a list
  cluster_list <- lapply(cluster_names, function(cluster) row.names(cseq[[cluster_name]][which(cseq[[cluster_name]][,1] == cluster),, drop = F]))
  names(cluster_list) <- cluster_names
  ### All possible filter types
  if(!filter_type %in% c('mincells', 'expressionflat', 'expressionratio', 'negative', 'positive')){
    warning(paste0('Input filter type ', filter_type,' is not among the possible options. Skipping'))
    break
  }
  if(filter_type == 'mincells'){
    out <- too_few_cells(cseq, cluster_name, value)
  } else if(filter_type == 'expressionflat'){
    out <- too_few_reads(cseq, datatype, cluster_list, value)
  } else if (filter_type == 'expressionratio'){
    out <- too_few_reads_relative(cseq, datatype, cluster_list, value)
  } else if(filter_type == 'negative'){
    out<- too_many_negatives(cseq, datatype, cluster_list, features, value)
  } else if (filter_type == 'positive'){
    out <- too_few_positives(cseq, datatype, cluster_list, features, value)
  }
  ### Account for character(0) being output if no clusters match criteria. ifelse statment selects only the 1st element if there are more than one cluster.
  if(length(out) == 0){
    clusters <- ''
  } else{
    clusters <- names(out)
  }
  cell_ids <- unlist(lapply(names(out), function(cluster) cluster_list[[cluster]]))
  ncells <- length(cell_ids)
  res <- list('clusters' = clusters, 'N' = ncells, 'cell_ids' = cell_ids)
  return(res)
}

make_test_rds <- function(sdat_file){
  sdat <- readRDS(sdat_file)
  csub <- subset(x = sdat, downsample = 1000)
  print(dim(csub))
  saveRDS(csub, file = gsub('.RDS', '_TempTestSubset.RDS', sdat_file))
}

read_test_rds <- function(sdat_file){
  sdat <- readRDS(gsub('.RDS', '_TempTestSubset.RDS', sdat_file))
  return(sdat)
}
