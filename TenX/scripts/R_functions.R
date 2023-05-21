
### Creates an MTX object from the TPM values in the input gene_abundances.tab files
tab2mtx <- function(tabs, cell_ids = NA){
  ### the column names have not been specified, attempt to pull cell_ids from the file paths (animal, timepoint, plate, lane, cell)
  if(is.na(cell_ids[1])){
    cns <- sapply(tabs, function(tab) id_from_path(tab, plate = T))
  } else{
    cns <- cell_ids
  }
  ### Read from file
  mat <- read.csv(tabs[1], sep = '\t', stringsAsFactors = F)
  ### For the few cases of duplicate gene IDs, append a 'v1', 'v2', ...
  mat <- make_unique_GeneID(mat)
  ### Select only Gene.ID and TPM and rename columns
  mat0 <- mat[, c('Gene.ID', 'TPM'), drop = F]
  colnames(mat0) <- c('Gene.ID', cns[1])
  for(i in 2:length(tabs)){
    #for(tab in tabs[2:length(tabs)]){
    mat <- read.csv(tabs[i], sep = '\t', stringsAsFactors = F)
    mat <- make_unique_GeneID(mat)
    mat <- mat[, c('Gene.ID', 'TPM'), drop = F]
    colnames(mat) <- c('Gene.ID', cns[i])
    mat0 <- merge(mat0, mat, by = 'Gene.ID', all = T)
  }
  
  ### Put Gene.ID as row.names and remove its column
  row.names(mat0) <- mat0[, 'Gene.ID']
  mat0 <- mat0[, colnames(mat0)[which(!colnames(mat0) %in% c("Gene.ID"))]]
  #mat0 <- mat0[complete.cases(mat0), ]
  mat0[is.na(mat0)] <- 0
  return(mat0)
}


star_files_to_df <- function(star_files){
  stars <- lapply(star_files, function(star_file) read_star_summary(star_file))
  output <- as.data.frame(matrix(ncol = 3, nrow = length(stars)))
  colnames(output) <- c('N_reads', 'N_mapped', 'N_multimapped')
  output[, 'N_reads'] <- as.numeric(sapply(stars, function(s) s['Number of input reads', 'value']))
  output[, 'N_mapped'] <- as.numeric(sapply(stars, function(s) s['Uniquely mapped reads number', 'value']))
  output[, 'N_multimapped'] <- as.numeric(sapply(stars, function(s) s['Number of reads mapped to multiple loci', 'value']))
  return(output)
}

### Read values from trimmomatic log file
read_trim_summary <- function(trim_file, cnames = c("Input_Read_Pairs","Both_Surviving","Forward_Only_Surviving","Reverse_Only_Surviving","Dropped")){
  #print(trim_file)
  if(!file.exists(trim_file)){
    out <- rep('Absent', length(cnames))
    names(out) <- cnames
  } else{
    txt <- readLines(trim_file, skipNul = T, encoding = 'UTF-8')
    ### If error
    if(any(grepl('Error', txt))){
      out <- rep('Error', length(cnames))
      names(out) <- cnames
    } else{
      line <- strsplit(txt[grep('Input Read Pairs', txt)], '\t')[[1]]
      if(any(grepl('NaN', line))){
        out <- rep('Empty', length(cnames))
        names(out) <- cnames
      } else {
        line <- gsub(': ', '', line)
        
        pieces <- stri_split_regex(line, pattern = '\\([^A-Z]*\\%\\)')[[1]]
        pieces <- pieces[pieces != '']
        out <- unlist(lapply(pieces, function(x) split_num_letter_pairs(x, char = T)))
      }
    }
  }
  return(out)
}


mtx2seurat <- function(mtx, sample_info, stars, gtf, project, covs = c('sample', 'animal', 'timepoint', 'plate', 'flowcell')){
  my_cells <- CreateSeuratObject(mtx, project = project, assay = "RNA", min.cells = 2)
  covs <- covs[which(covs %in% colnames(sample_info))]
  for(cov in covs){
    my_cells[[cov]] <- sample_info[[cov]]
  }
  my_cells[['MT_sum']] <- sum_chromosome(mtx, chr ='MT', gtf)
  # my_cells[['N_reads']] <- stars[, 'N_reads']
  # my_cells[['N_mapped']] <- stars[, 'N_mapped']
  # my_cells[['N_multimapped']] <- stars[, 'N_multimapped']
  return(my_cells)
}


### Pairwise covariate association tests
pairwise_associations <- function(cov_data, pc, symetric = T, plots_path, ID = 'cell_id', covs_numeric = c('nCount_RNA', 'nFeature_RNA'), covs_factor = c('batch', 'seurat_clusters'), covs_binary = c()){
  ### Remove those that are not actually in cov_data
  covs_numeric <- covs_numeric[which(covs_numeric %in% colnames(cov_data))]
  covs_factor <- covs_factor[which(covs_factor %in% colnames(cov_data))]
  covs_binary <- covs_binary[which(covs_binary %in% colnames(cov_data))]
  ### Remove those that only have 1 level
  covs_numeric <- covs_numeric[sapply(covs_numeric, function(x) length(unique(cov_data[, x])) > 1)]
  covs_factor <- covs_factor[sapply(covs_factor, function(x) length(unique(cov_data[, x])) > 1)]
  covs_binary <- covs_binary[sapply(covs_binary, function(x) length(unique(cov_data[, x])) > 1)]
  
  covs_all = c(covs_numeric, covs_factor, covs_binary)
  cov_combos = combn(covs_all, 2, simplify = F)
  
  ### Define the access to the test results p-value or other 
  if(pc == 'p'){
    result_code = c('p.value', 'p.value', 'prob', 'Pr(>F)')
    names(result_code) = c('Chisq', 'Spearman', 'logistf', 'anova')
  } else if(pc == 'c'){
    result_code = c('statistic', 'estimate', 'coefficients', 'F value')
    names(result_code) = c('Chisq', 'Spearman', 'logistf', 'anova')
  } else(print('pc should be either "p" or "c"'))
  ### Create the first empty data frame
  res <- data.frame(matrix(NA, nrow = length(covs_all), ncol = length(covs_all)))
  colnames(res) = covs_all
  row.names(res) = covs_all
  test_names <- res
  print(paste0('There are ', length(cov_combos), ' pairwise combinations of covariates.'))
  for(cc in 1:length(cov_combos)){
    cov_a = cov_combos[[cc]][1]
    cov_b = cov_combos[[cc]][2]
    ### Remove samples where either covariate is NA
    #print(paste0(cov_b,' and ',cov_a))
    cov_data_sub <- cov_data[complete.cases(cov_data[, c(cov_a,cov_b)]), ]
    #print(paste0(cov_b,' and ',cov_a, ' exclude ', length(cov_data[, ID]) - length(cov_data_sub[, ID]), ' samples due to NAs.'))
    
    ### Chisq if both are factors or binary
    if(cov_a %in% c(covs_factor, covs_binary) & cov_b %in% c(covs_factor, covs_binary)){
      #print('Chi-square')
      r <- chisq.test(table(cov_data_sub[,cov_a],cov_data_sub[,cov_b]))
      res[cov_a,cov_b] <- r[result_code['Chisq']]
      test_names[cov_a, cov_b] <- 'Chisq'
    } else if(cov_a %in% covs_numeric & cov_b %in% covs_numeric){
      ###  Spearman correlation if both are numeric
      #print('Spearman correlation')
      r <- cor.test(x = cov_data_sub[,c(cov_a)], y = cov_data_sub[,c(cov_b)],method='spearman')
      res[cov_a,cov_b] <- r[result_code['Spearman']]
      test_names[cov_a, cov_b] <- 'Spearman'
      png(paste0(plots_path,cov_a,'_',cov_b,'.png'))
      plot(x = cov_data_sub[,c(cov_a)], y = cov_data_sub[,c(cov_b)], main = paste0(cov_a,' v. ', cov_b), xlab = cov_a, ylab = cov_b) 
      dev.off()
    } else if(cov_a %in% covs_binary & cov_b %in% covs_numeric){
      #print('Logistic regression')
      ### Firth logistic regression if cov_a is binary and cov_b is numeric
      r <- logistf(data = cov_data_sub, formula = eval(parse(text=cov_a)) ~ eval(parse(text=cov_b)))
      res[cov_a,cov_b] <- r[result_code['logistf']][[1]][2]
      test_names[cov_a, cov_b] <- 'Logistic Regression'
      png(paste0(plots_path,cov_b,'_',cov_a,'.png'))
      boxplot(cov_data_sub[,c(cov_b)] ~ cov_data_sub[,c(cov_a)], main = paste0(cov_b,' v. ', cov_a), xlab = cov_a, ylab = cov_b)
      dev.off()
    } else if(cov_b %in% covs_binary & cov_a %in% covs_numeric){
      #print('Logistic regression')
      ### Firth logistic regression if cov_b is binary and cov_a is numeric
      r <- logistf(data = cov_data_sub, formula = eval(parse(text=cov_b)) ~ eval(parse(text=cov_a)))
      res[cov_a, cov_b] <- r[result_code['logistf']][[1]][2]
      test_names[cov_a, cov_b] <- 'Logistic Regression'
      png(paste0(plots_path,cov_a,'_',cov_b,'.png'))
      boxplot(cov_data_sub[,c(cov_a)] ~ cov_data_sub[,c(cov_b)], main = paste0(cov_a,' v. ', cov_b), xlab = cov_b, ylab = cov_a) 
      dev.off()
    } else if(cov_a %in% covs_factor & cov_b %in% covs_numeric){
      ### ANOVA if cov_a is factor and cov_b is numeric
      #print('ANOVA')
      lm_res <- lm(formula = eval(parse(text=cov_b)) ~ eval(parse(text=cov_a)), data = cov_data_sub)
      anova_res <- anova(lm_res)
      res[cov_a, cov_b] <- anova_res[result_code['anova']][[1]][1]
      test_names[cov_a, cov_b] <- 'ANOVA'
      #tukey_res <- TukeyHSD(x=aov(lm_res), conf.level=0.95)
      
      png(paste0(plots_path, cov_b, '_', cov_a, '.png'))
      boxplot(cov_data_sub[ ,c(cov_b)] ~ cov_data_sub[, c(cov_a)], main = paste0(cov_b, ' v. ', cov_a), xlab = cov_a, ylab = cov_b)
      dev.off()
    } else if(cov_b %in% covs_factor & cov_a %in% covs_numeric){
      #print('ANOVA')
      ### ANOVA if cov_b is factor and cov_a is numeric
      lm_res <- lm(formula = eval(parse(text=cov_a)) ~ eval(parse(text=cov_b)), data = cov_data_sub)
      anova_res <- anova(lm_res)
      res[cov_a, cov_b] <- anova_res[result_code['anova']][[1]][1]
      test_names[cov_a, cov_b] <- 'ANOVA'
      #tukey_res <- TukeyHSD(x=aov(lm_res), conf.level=0.95)    
      png(paste0(plots_path, cov_a,'_', cov_b, '.png'))
      boxplot(cov_data_sub[, c(cov_a)] ~ cov_data_sub[, c(cov_b)], main = paste0(cov_a,' v. ', cov_b), xlab = cov_b, ylab = cov_a) 
      dev.off()
    }
    if(symetric == 'N'){
      res[cov_b, cov_a] <- length(cov_data_sub[, ID])
      test_names[cov_b, cov_a] <- length(cov_data_sub[, ID])
    } else if(symetric == T){
      res[cov_b, cov_a] <- res[cov_a, cov_b]
      test_names[cov_b, cov_a] <- test_names[cov_a, cov_b]
    }
  }
  return(list('res' = res, test = test_names))
}


run_ilisi <- function(sdat, reduction, batch_feature){
  ilisi <- compute_lisi(sdat@reductions[[reduction]]@cell.embeddings, sdat@meta.data[, batch_feature, drop = F], c(batch_feature))
  ilisi$Integrated <- reduction
  return(ilisi)
}

plot_ilisi <- function(sdat, reductions = c('umap', 'umap0'), batch_feature = 'plate'){
  ### cehck that reductions are actually in sdat
  if(!all(reductions %in% names(sdat@reductions))){
    warning("Some input reductions are not in the input seurat object.")
  }
  ### Reduce to those that are present
  reductions <- reductions[reductions %in% names(sdat@reductions)]
  ilisi <- lapply(reductions, function(reduction) 
    run_ilisi(sdat, reduction,batch_feature))
  
  
  ### Bind the two together for plotting
  ilisi <- rbindlist(ilisi)
  colnames(ilisi) <- c('iLISI', 'Integrated')
  ### Calculate medians of each
  med <- ddply(ilisi, "Integrated", summarise, grp.median = median(iLISI))
  ### Density plot
  p <- ggplot(ilisi, aes(x = iLISI, color = Integrated)) +
    geom_density()+
    geom_vline(data = med, aes(xintercept = grp.median, color = Integrated),
               linetype = "dashed")
  return(p)
}


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


bcl2fastq_count <- function(cfile, sample_info){
  count <- read.table(cfile, sep = ',', header = T)
  sample_sub <- sample_info[which(sample_info$count_files == cfile), ]
  row.names(sample_sub) <- sample_sub$UDP_ID
  ### Order by count, so they match
  sample_sub <- sample_sub[count$key,]
  ### Get Cell_ID in count
  count$Cell_ID <- sample_sub[count$key, 'Cell_ID']
  ### Drop any that aren't in sampleinfo
  count <- count[!is.na(count$Cell_ID),]
  count$FCID <- sample_sub[count$key, 'FCID']
  count$Lane <- sample_sub[count$key, 'Lane']
  count$Sample_Name <- sample_sub[count$key, 'Sample_Name']
  row.names(count) <- count$Cell_ID
  count <- count[, c('Cell_ID', 'FCID', 'Lane', 'key', 'count')]
  colnames(count) <- c('Cell_ID', 'FCID', 'Lane', 'UDP_ID', 'bcl2fastq_count')
  return(count)
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

