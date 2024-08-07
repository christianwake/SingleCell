library('sys')
library('ggplot2')
library('stringr')
library('GenomicRanges')
library('dplyr')
library('readr')
library('effsize')
library('data.table')
library('biomaRt')
library('GSEABase')
library('readxl')
library('viridis')
library('Seurat')
library('circlize')
library('ComplexHeatmap')

source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/DE_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')

if(interactive()){
  project <- '2021600_kristin'
  qc_name <- 'Run2023-05-14'
  de_files <- Sys.glob(file.path(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', 
                                           qc_name, '/DE/*/*/*DE*.tsv')))
  sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, 
                      '/Mapped.RDS')
  
  gtf_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/gtf.RDS')
  out_pdf1 <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, 
                     '/DE/Bubble_plots.pdf') 
  out_pdf1 <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, 
                     '/DE/DE_Heatmap_within_strat.pdf') 
  out_pdf2 <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, 
                     '/DE/DE_Heatmap_within_test.pdf') 
  
} else{
  args = commandArgs(trailingOnly=TRUE)
  out_pdf1 <- args[1]
  out_pdf2 <- args[2]
  sdat_file <- args[3]
  gtf_file <- args[4]
  de_files <- args[5:length(args)]
}

### Read counts
sdat <- readRDS(sdat_file)

### Do direct matching from gtf
if(grepl('\\.gtf', gtf_file)){
  gtf <- read_gtf(gtf_file, feature_type = 'gene', atts_of_interest = c('gene_id', 'gene_name',  'gene_biotype'))
} else if(grepl('\\.RDS', gtf_file)){
  gtf <- readRDS(gtf_file)
} else{
  warning('No gtf file')
}
row.names(gtf) <- gtf$gene_id


de_list <- lapply(de_files, function(de_file) read.table(de_file, header = T, check.names = F, stringsAsFactors = F, na.strings = c("", "NA"), sep = '\t'))
names(de_list) <- gsub('.tsv', '', basename(de_files))

test_chunk <- sapply(names(de_list), function(x) strsplit(x, '_DE')[[1]][1])
test_name <- sapply(test_chunk, function(x) strsplit(x, '-')[[1]][1])
test_value1 <- sapply(test_chunk, function(x) strsplit(x, '-')[[1]][2])
test_value2 <- sapply(test_chunk, function(x) strsplit(x, '-')[[1]][3])
### For those whose values are NA (because it was binary) figure out which is 1 and which is 2
is <- which(test_value1 == 'NA' | test_value2 == 'NA')
for(i in is){
  Idents(sdat) <- sdat[[test_name[1]]]
  idents <- levels(Idents(sdat))
  test_value1[i] <- idents[1]
  test_value2[i] <- idents[2]
}

strat1_chunk <- sapply(names(de_list), function(x) strsplit(strsplit(x, 'Strat1-')[[1]][2], '_Strat2')[[1]][1])
strat_name1 <- sapply(strat1_chunk, function(x) strsplit(x, '-')[[1]][1])
strat_value1A <- sapply(strat1_chunk, function(x) strsplit(x, '-')[[1]][2])
strat_value1B <- sapply(strat1_chunk, function(x) strsplit(x, '-')[[1]][3])

strat2_chunk <- sapply(names(de_list), function(x) strsplit(x, '_Strat2-')[[1]][2])
strat_name2 <- sapply(strat2_chunk, function(x) strsplit(x, '-')[[1]][1])
strat_value2A <- sapply(strat2_chunk, function(x) strsplit(x, '-')[[1]][2])
strat_value2B <- sapply(strat2_chunk, function(x) strsplit(x, '-')[[1]][3])

### Select those for which the test name is the same, all stratifcations are the same, but the test values are different.
### Unique combinations of test name and all stratifications
test_strat_group <- paste0(test_name, '_', strat1_chunk, '_', strat2_chunk)
test_strat_group <- sapply(unique(test_strat_group), function(x) which(test_strat_group == x))
### Heatmap can't handle just 1 column, so remove those with only 1
test_strat_group <- test_strat_group[which(sapply(test_strat_group, function(x) length(x) > 1))]
print('Same stratifications, different tests values')
pdf(out_pdf1)
### For each unique combination of test name and all strats,
for(test_strat in names(test_strat_group)){
  print(test_strat)
  is <- unlist(test_strat_group[test_strat])
  des <- de_list[is]
  i <- is[1]
  
  ### Reduce names but keep lost info in a plot description for title
  names(des) <- sapply(names(des), function(x) gsub(paste0('_DE_Strat1-', strat1_chunk[i]), '', x))
  names(des) <- sapply(names(des), function(x) gsub(paste0('_Strat2-', strat2_chunk[i]), '', x))
  names(is) <- names(des)

  plot_title <- paste0(test_name[i], ' tests')
  ### Add strat info to plot title
  if(strat_name1[i] != 'All'){
    value <- ifelse(strat_value1A[i] == strat_value1B[i], 
                    strat_value1A[i], 
                    paste0(strat_value1A[i], '-', strat_value1B[i]))
    plot_title <- paste0(plot_title, ' in ', strat_name1[i], '-', value)
  }
  if(strat_name2[i] != 'All'){
    value <- ifelse(strat_value2A[i] == strat_value2B[i], 
                    strat_value2A[i], 
                    paste0(strat_value2A[i], '-', strat_value2B[i]))
    plot_title <- paste0(plot_title, ' in ', strat_name2[i], '-', value)
  }
  
  ### Set and adjust p value threshold
  fc_threshold <- 0
  pthresh <- 0.05
  gene_ids <- unique(unlist(lapply(des, function(de) row.names(de[which(de$p_val_adj < pthresh), ]))))
  while(length(gene_ids) > 50){
    pthresh <- pthresh/10
    gene_ids <- unique(unlist(lapply(des, function(de) row.names(de[which(de$p_val_adj < pthresh), ]))))
  }
  ### Correct if threshold reduced too much
  if(length(gene_ids) == 0){
    pthresh <- pthresh * 10
    gene_ids <- unique(unlist(lapply(des, function(de) row.names(de[which(de$p_val_adj < pthresh), ]))))
  }
  
  genes <- gtf[gene_ids, 'gene_name']
  names(genes) <- gene_ids
  
  ### Remove test name from names(des) (which will be the heatmap column names)
  names(des) <- gsub(paste0(test_name[i], '-'), '', names(des))
  
  labs <- names(des)
  ### If any test values are combinations of multiple groups (has '+')
  if(any(grepl('\\+', names(des)))){
    ### Check if they are equal to all the groups without a +
    withplus <- unique(labs)[grepl('\\+', unique(labs))]
    sansplus <- unique(labs)[!grepl('\\+', unique(labs))]
    withplus <- lapply(withplus, function(x) strsplit(x, '-')[[1]])
    ### Just value1
    sansplus1 <- lapply(sansplus, function(x) strsplit(x, '-')[[1]][1])
    sansplus1 <- unique(unlist(sansplus))
    
    sets1 <-lapply(withplus, function(x) strsplit(x[1], '\\+')[[1]])
    sets2 <-lapply(withplus, function(x) strsplit(x[2], '\\+')[[1]])
    for(ii in 1:length(sets1)){
      a <- sansplus1[which(!sansplus1 %in% sets1[[ii]])]
      if(all(union(sets2[[ii]], a) == intersect(sets2[[ii]],a))){
        labs[which(labs == paste0(withplus[[ii]][1], '-', withplus[[ii]][2]))] <-
          paste0(withplus[[ii]][1], '-All')
      }
    }
    names(des) <- labs
  }
  
  fontsize <- heatmap_N_to_size(length(genes))
  row_title <- paste0('DE padj < ', signif(pthresh, digits = 3), ', logFC > ', fc_threshold)
  #column_labels <- c('testB', 'testC')
  ### genes is an array of gene_names (for plotting) with names as IDs as appear in de_list
  p <- my_heatmap(des, genes, categories = NA, fontsize = fontsize, plot_title, row_title, 
                      pheatmap = F, col_clust = F, sigs = F) 
  
  #print(p)
}
dev.off()

### first check if no stratification (strat_name == 'All'), then if the values are equal
reduce_strat_names <- function(sname, i){
  if(strat_name1[i] == 'All'){
    sname <- gsub('Strat1-All-All-All', 'No_strat', sname)
  } else if(strat_value1A[i] == strat_value1B[i]){
    sname <- gsub(paste0('Strat1-', strat1_chunk[i]), paste0('Strat1-', strat_name1[i], '-', strat_value1A[i]), sname)
  }
  if(strat_name2[i] == 'All'){
    sname <- gsub('_Strat2-All-All-All', '', sname)
  } else if(strat_value2A[i] == strat_value2B[i]){
    sname <- gsub(paste0('Strat2-', strat2_chunk[i]), paste0('Strat2-', strat_name2[i], '-', strat_value2A[i]), sname)
  }
  return(sname)
}

### Select those for which the test name is the same, test values are the same, but all stratifcations are different, 
### Unique combinations of test name and test_values
test_val_group <- paste0(test_name, '_', test_value1, '_', test_value2)
test_val_group <- sapply(unique(test_val_group), function(x) which(test_val_group == x))
### Heatmap can't handle just 1 column, so remove those with only 1
test_val_group <- test_val_group[which(sapply(test_val_group, function(x) length(x) > 1))]

print('Same tests, different stratifications')
pdf(out_pdf2)
### For each unique combination of test name and test values,
for(test_val in names(test_val_group)){
  print(test_val)
  is <- unlist(test_val_group[test_val])
  des <- de_list[is]
  
  ##### Reduce names
  ### Remove beginning test chunk
  names(des) <- sapply(names(des), function(x) gsub(paste0('^', gsub('+', '\\+', test_chunk[is[1]], fixed = T), '_DE_'), '', x))
  ### first check if no stratification (strat_name == 'All'), then if the values are equal
  names(des) <- sapply(1:length(is), function(j) reduce_strat_names(names(des[j]), is[j]))
  names(is) <- names(des)
  
  plot_title <- paste0(test_chunk[is[1]], ' tests')
  
  
  ### Set and adjust p value threshold
  pthresh <- 0.05
  gene_ids <- unique(unlist(lapply(des, function(de) row.names(de[which(de$p_val_adj < pthresh), ]))))
  while(length(gene_ids) > 50){
    pthresh <- pthresh/10
    gene_ids <- unique(unlist(lapply(des, function(de) row.names(de[which(de$p_val_adj < pthresh), ]))))
  }
  ### Correct if threshold reduced too much
  if(length(gene_ids) < 2){
    pthresh <- pthresh * 10
    gene_ids <- unique(unlist(lapply(des, function(de) row.names(de[which(de$p_val_adj < pthresh),]))))
  }
  genes <- gtf[gene_ids, 'gene_name']
  names(genes) <- gene_ids
  
  labs <- names(des)
  ### If any test values are combinations of multiple groups (has '+')
  if(any(grepl('\\+', labs))){
    ### Check if they are equal to all the groups without a +
    withplus <- unique(labs)[grepl('\\+', unique(labs))]
    sansplus <- unique(labs)[!grepl('\\+', unique(labs))]
    withplus <- lapply(withplus, function(x) strsplit(x, '-')[[1]])
    ### Just value1
    sansplus1 <- lapply(sansplus, function(x) strsplit(x, '-')[[1]][1])
    sansplus1 <- unique(unlist(sansplus))
    
    sets1 <-lapply(withplus, function(x) strsplit(x[1], '\\+')[[1]])
    sets2 <-lapply(withplus, function(x) strsplit(x[2], '\\+')[[1]])
    for(ii in 1:length(sets1)){
      a <- sansplus1[which(!sansplus1 %in% sets1[[ii]])]
      if(all(union(sets2[[ii]], a) == intersect(sets2[[ii]],a))){
        labs[which(labs == paste0(withplus[[ii]][1], '-', withplus[[ii]][2]))] <-
          paste0(withplus[[ii]][1], '-All')
      }
    }
    names(des) <- labs
  }
  
  fontsize <- heatmap_N_to_size(length(genes))
  row_title <- paste0('DE padj < ', signif(pthresh, digits = 3), ', logFC > ', fc_threshold)
  column_labels <- c('testB', 'testC')
  ### genes is an array of gene_names (for plotting) with names as IDs as appear in de_list
  p <- my_heatmap(des, genes, categories = NA, fontsize = fontsize, plot_title, row_title, 
                  pheatmap = F, col_clust = F, sigs = F) 
  #print(p)
} 

dev.off()
