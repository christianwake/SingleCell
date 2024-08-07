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
library('viridis')
library('Seurat')

source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/DE_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')

if(interactive()){
  project <- '2021600_kristin'
  qc_name <- 'Run2023-05-14'
  fgsea_files <- Sys.glob(file.path(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', 
                                           qc_name, '/DE/*/*/*GSEA*.tsv')))
  sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, 
                      '/Mapped.RDS')
  
  gtf_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/gtf.RDS')
  out_pdf1 <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, 
                    '/DE/Bubble_plots.pdf') 
  out_pdf1 <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, 
                     '/DE/GSEA_Bubble_within_strat.pdf') 
  out_pdf2 <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, 
                     '/DE/GSEA_Bubble_within_test.pdf') 
  
} else{
  args = commandArgs(trailingOnly=TRUE)
  out_pdf1 <- args[1]
  out_pdf2 <- args[2]
  sdat_file <- args[3]
  gtf_file <- args[4]
  fgsea_files <- args[5:length(args)]
}

font <- 'Times'
fig_theme = theme_bw() + 
  theme(text = element_text(family = font, color = "black"),
        axis.title = element_text(family = font, size = 8, color = "black"),
        axis.text  = element_text(family = font, size = 8, color = "black"),
        legend.title = element_text(family = font, size = 7, color = "black"),
        legend.text  = element_text(family = font, size = 7, color = "black"),
        legend.spacing = unit(5,'points'), 
        legend.key.size = unit(8, 'points'),
        legend.box.spacing = unit(0, 'points'),
        axis.ticks.length = unit(1, 'points'))

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


fgseas_list <- lapply(fgsea_files, function(fgsea_file) read.table(fgsea_file, header = T, 
                check.names = F, stringsAsFactors = F, na.strings = c("", "NA"), sep = '\t'))
names(fgseas_list) <- gsub('.tsv', '', basename(fgsea_files))

test_chunk <- sapply(names(fgseas_list), function(x) strsplit(x, '_GSEA')[[1]][1])
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


strat1_chunk <- sapply(names(fgseas_list), function(x) strsplit(strsplit(x, 'Strat1-')[[1]][2], '_Strat2')[[1]][1])
strat_name1 <- sapply(strat1_chunk, function(x) strsplit(x, '-')[[1]][1])
strat_value1A <- sapply(strat1_chunk, function(x) strsplit(x, '-')[[1]][2])
strat_value1B <- sapply(strat1_chunk, function(x) strsplit(x, '-')[[1]][3])

strat2_chunk <- sapply(names(fgseas_list), function(x) strsplit(x, '_Strat2-')[[1]][2])
strat_name2 <- sapply(strat2_chunk, function(x) strsplit(x, '-')[[1]][1])
strat_value2A <- sapply(strat2_chunk, function(x) strsplit(x, '-')[[1]][2])
strat_value2B <- sapply(strat2_chunk, function(x) strsplit(x, '-')[[1]][3])

### Select those for which the test name is the same, all stratifcations are the same, but the test values are different.
### Unique combinations of test name and all stratifications
test_strat_group <- paste0(test_name, '_', strat1_chunk, '_', strat2_chunk)
test_strat_group <- sapply(unique(test_strat_group), function(x) which(test_strat_group == x))

print('Same stratifications, different tests values')
pdf(out_pdf1)
### For each unique combination of test name and all strats,
for(test_strat in names(test_strat_group)){
  print(test_strat)
  is <- unlist(test_strat_group[test_strat])
  fgseas <- fgseas_list[is]
  i <- is[1]
  
  ### Reduce names but keep lost info in a plot description for title
  names(fgseas) <- sapply(names(fgseas), function(x) gsub(paste0('_GSEA_Strat1-', strat1_chunk[i]), '', x))
  names(fgseas) <- sapply(names(fgseas), function(x) gsub(paste0('_Strat2-', strat2_chunk[i]), '', x))
  names(is) <- names(fgseas)
  #names(fgseas) <- sapply(names(fgseas), function(x) )

  plot_title <- paste0(test_name[i], ' tests')
  ### Add strat info
  if(strat_name1[i] != 'All'){
    value <- ifelse(strat_value1A[i] == strat_value1B[i], strat_value1A[i], paste0(strat_value1A[i], '+', strat_value1B[i]))
    plot_title <- paste0(plot_title, ' in ', strat_name1[i], '-', value)
  }
  if(strat_name2[i] != 'All'){
    value <- ifelse(strat_value2A[i] == strat_value2B[i], strat_value2A[i], paste0(strat_value2A[i], '+', strat_value2B[i]))
    plot_title <- paste0(plot_title, ' in ', strat_name2[i], '-', value)
  }
  
  ### Set and adjust p value threshold
  pthresh <- 0.05
  pathways <- unique(unlist(lapply(fgseas, function(fgsea) fgsea[which(fgsea$padj < pthresh), 'pathway'])))
  while(length(pathways) > 50){
    pthresh <- pthresh/10
    pathways <- unique(unlist(lapply(fgseas, function(fgsea) fgsea[which(fgsea$padj < pthresh), 'pathway'])))
  }
  ### Correct if threshold reduced too much
  if(length(pathways) == 0){
    pthresh <- pthresh * 10
    pathways <- unique(unlist(lapply(fgseas, function(fgsea) fgsea[which(fgsea$padj < pthresh), 'pathway'])))
  }
  plot_dat <- rbindlist(fgseas, idcol = 'fgsea_res')
  plot_dat <- plot_dat[which(plot_dat$pathway %in% pathways),]
  plot_dat$fgsea_res <- factor(plot_dat$fgsea_res)
  
  ### Remove insignificant
  #plot_dat <- plot_dat[which(plot_dat$padj < 0.05), ]

  
  y_title <- paste0('pathways, adj p < ', pthresh)
  ### Find if there is a common prefix
  if(length(unique(sapply(pathways, function(x) substr(x, 1, 1)))) == 1){
    n <- 1
    while(length(unique(sapply(pathways, function(x) substr(x, 1, n)))) == 1){
      n <- n + 1
    }
    prefix <- unique(sapply(pathways, function(x) substr(x, 1, n - 1)))
    y_title <- paste0(prefix, y_title)
    ###  Remove refix in pathways array
    pathways <- sapply(pathways, function(x) gsub(paste0('^', prefix), '', x))
    ### Remove prefix in plot data
    plot_dat$pathway <- sapply(plot_dat$pathway, function(x) gsub(paste0('^', prefix), '', x))
  }
  
  ### Handle long pathway names
  ### manually eyeball gene names in the pathway names that aren't directly matched in the gtf file
  gene_names <- c('TH17', 'RANKLRANK', 'IL12', 'TNFR2', 'NTHI', 'AP1', 'TAK1', 'IKK', 'SASP', 'TXA2',
                  'STAT5', 'MTORC1', 'UV', 'TNFA', 'NFKB')
  plot_dat$pathway <- sapply(plot_dat$pathway, function(x) abbrev_pathway_name(x, gtf, gene_names))
  plot_dat <- as.data.frame(plot_dat)
  ### Order pathways
  #a<-hclust(dist(plot_dat$NES))
  
  pathways <- unique(plot_dat$pathway)
  pathways <- pathways[order(sapply(pathways, function(pathway) mean(abs(plot_dat[which(plot_dat$pathway == pathway), 'NES']))), decreasing = T)]
  
  plot_dat$pathway <- factor(plot_dat$pathway, levels = pathways)
  ### Remove test name from the display
  plot_dat$fgsea_res <- gsub(paste0(test_name[i], '-'), '', plot_dat$fgsea_res)
  
  legend_size <-7
  x_size <- 10
  y_size <- 9
  ### If any test values are combinations of multiple groups (has '+')
  if(any(grepl('\\+', plot_dat$fgsea_res))){
    ### Check if they are equal to all the groups without a +
    withplus <- unique(plot_dat$fgsea_res)[grepl('\\+', unique(plot_dat$fgsea_res))]
    sansplus <- unique(plot_dat$fgsea_res)[!grepl('\\+', unique(plot_dat$fgsea_res))]
    withplus <- lapply(withplus, function(x) strsplit(x, '-')[[1]])
    ### Just value1
    sansplus1 <- lapply(sansplus, function(x) strsplit(x, '-')[[1]][1])
    sansplus1 <- unique(unlist(sansplus))
    
    sets1 <-lapply(withplus, function(x) strsplit(x[1], '\\+')[[1]])
    sets2 <-lapply(withplus, function(x) strsplit(x[2], '\\+')[[1]])
    for(ii in 1:length(sets1)){
      a <- sansplus1[which(!sansplus1 %in% sets1[[ii]])]
      if(all(union(sets2[[ii]], a) == intersect(sets2[[ii]],a))){
        plot_dat[which(plot_dat$fgsea_res == paste0(withplus[[ii]][1], '-', withplus[[ii]][2])), 'fgsea_res'] <-
          paste0(withplus[[ii]][1], '-All')
      }
    }
  }
  
  p <- ggplot(data = plot_dat, aes(y = pathway, x = fgsea_res, size = -log10(padj), fill = NES)) +
    ggtitle(plot_title) + 
    labs(y = y_title) + 
    geom_point(shape = 21) + 
    scale_size_continuous(range = c(0,4)) +
    scale_fill_gradientn(colours = viridis(100)) +
    fig_theme +
    scale_x_discrete(position = "top") + 
    scale_y_discrete(position = "right") + 
    theme(axis.title.x = element_blank()) +
    labs(size = '-log10 adj. p value', fill = 'Normalized\nEnrichment Score') +
    theme(legend.position = "bottom",
          legend.title = element_text(colour = "black"),#, size = legend_size),
          legend.spacing = unit(0, 'points'), legend.key.size = unit(8, 'points'), 
          legend.direction = 'horizontal', legend.box = 'horizontal',
          legend.box.spacing = unit(1, 'points'),
          legend.box.margin = margin(t = 0, r = -10, b = 0, l = 10, unit = "point")) +
 
    guides(shape = guide_legend(override.aes = list(size = 8))) + 
    guides(color = guide_legend(override.aes = list(size = 8))) + 
    theme(panel.spacing = unit(1, "lines")) +
    guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5), size = guide_legend(title.position = "top", title.hjust = 0.5))
  
  print(p)
#   ### Print one example Vln plot of mean expressions, for NES direction reference
#   # max_pathway <- 'Inflammatory Response'
#   # max_test <- plot_dat[which(plot_dat$pathway == max_pathway), 'fgsea_res'][1]
#   max_pathway <- plot_dat[which(abs(plot_dat$NES) == max(abs(plot_dat$NES))), 'pathway']
#   max_test <- plot_dat[which(abs(plot_dat$NES) == max(abs(plot_dat$NES))), 'fgsea_res'][1]
#   gene_names <- strsplit(plot_dat[which(plot_dat$pathway == max_pathway), 'leadingEdgeStr'], split = ',')[[1]]
#   gene_ids <- gtf[which(gtf$gene_name %in% gene_names), 'gene_id']
# #  iii <- is[max_test]
#   test <- test_name[i]
#   ### Select the cells in the stratifications
#   if(strat_name1[i] == 'All'){
#     strat1 <- row.names(sdat@meta.data)
#   }else{
#     strat1 <- row.names(sdat@meta.data[
#                     which(sdat@meta.data[, strat_name1[i]] %in% c(strat_value1A[i], strat_value1B[i]))
#               ,])
#   }  
#   ### Select the cells in the stratifications
#   if(strat_name2[i] == 'All'){
#     strat2 <- row.names(sdat@meta.data)
#   }else{
#     strat2 <- row.names(sdat@meta.data[
#       which(sdat@meta.data[, strat_name2[i]] %in% c(strat_value2A[i], strat_value2B[i]))
#       ,])
#   }
#   
#   cell_ids <- intersect(strat1, strat2)
#   #### Subsest to those cells matching the test value
#   group1 <- row.names(sdat@meta.data[which(sdat@meta.data[, test_name[i]]  %in% test_value1[i] & 
#                                              sdat@meta.data[, 'Cell_ID'] %in% cell_ids), ])
#   group2 <- row.names(sdat@meta.data[which(sdat@meta.data[, test_name[i]]  %in% test_value2[i] & 
#                                              sdat@meta.data[, 'Cell_ID'] %in% cell_ids), ])
#   means1 <- as.data.frame(sapply(gene_ids, function(x) mean(sdat@assays$RNA@data[x, group1])))
#   means2 <- as.data.frame(sapply(gene_ids, function(x) mean(sdat@assays$RNA@data[x, group2])))
#   means1[, 2] <- row.names(means1)
#   means2[, 2] <- row.names(means2)
#   means1[, 3] <- test_value1[i]
#   means2[, 3] <- test_value2[i]
#   colnames(means1) <- c('means', 'gene', test_name[i])
#   colnames(means2) <- c('means', 'gene', test_name[i])
#   plot_dat2 <- rbind(means1, means2, make.row.names = F)
#   p <- ggplot(plot_dat2, aes(x = eval(parse(text = test_name[i])), y = means)) + 
#     geom_violin() + 
#     #geom_dotplot(binaxis='y', stackdir='center', dotsize =0.5) +
#     geom_boxplot(width=0.1) +
#     ggtitle(paste0(max_pathway, ' leading edge genes')) +
#     xlab(test_name[i])
#   print(p)
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

print('Same tests, different stratifications')
pdf(out_pdf2)
### For each unique combination of test name and test values,
for(test_val in names(test_val_group)){
  print(test_val)
  is <- unlist(test_val_group[test_val])
  fgseas <- fgseas_list[is]

  ##### Reduce names
  ### Remove beginning test chunk
  names(fgseas) <- sapply(names(fgseas), function(x) gsub(paste0('^', gsub('+', '\\+', test_chunk[is[1]], fixed = T), '_GSEA_'), '', x))
  ### first check if no stratification (strat_name == 'All'), then if the values are equal
  names(fgseas) <- sapply(1:length(is), function(j) reduce_strat_names(names(fgseas[j]), is[j]))


  names(is) <- names(fgseas)

  plot_title <- paste0(test_chunk[is[1]], ' tests')

  
  ### Set and adjust p value threshold
  pthresh <- 0.05
  pathways <- unique(unlist(lapply(fgseas, function(fgsea) fgsea[which(fgsea$padj < pthresh), 'pathway'])))
  while(length(pathways) > 50){
    pthresh <- pthresh/10
    pathways <- unique(unlist(lapply(fgseas, function(fgsea) fgsea[which(fgsea$padj < pthresh), 'pathway'])))
  }
  ### Correct if threshold reduced too much
  if(length(pathways) < 2){
    pthresh <- pthresh * 10
    pathways <- unique(unlist(lapply(fgseas, function(fgsea) fgsea[which(fgsea$padj < pthresh), 'pathway'])))
  }
  plot_dat <- rbindlist(fgseas, idcol = 'fgsea_res')
  plot_dat <- plot_dat[which(plot_dat$pathway %in% pathways),]
  plot_dat$fgsea_res <- factor(plot_dat$fgsea_res)
  
  ### Remove insignificant
  plot_dat <- plot_dat[which(plot_dat$padj < 0.05), ]
  
  
  y_title <- paste0('pathways, adj p < ', pthresh)
  ### Find if there is a common prefix
  if(length(unique(sapply(pathways, function(x) substr(x, 1, 1)))) == 1){
    n <- 1
    while(length(unique(sapply(pathways, function(x) substr(x, 1, n)))) == 1){
      n <- n + 1
    }
    prefix <- unique(sapply(pathways, function(x) substr(x, 1, n - 1)))
    y_title <- paste0(prefix, y_title)
    ###  Remove refix in pathways array
    pathways <- sapply(pathways, function(x) gsub(paste0('^', prefix), '', x))
    ### Remove prefix in plot data
    plot_dat$pathway <- sapply(plot_dat$pathway, function(x) gsub(paste0('^', prefix), '', x))
  }
  
  ### Handle long pathway names
  ### manually eyeball gene names in the pathway names that aren't directly matched in the gtf file
  gene_names <- c('TH17', 'RANKLRANK', 'IL12', 'TNFR2', 'NTHI', 'AP1', 'TAK1', 'IKK', 'SASP', 'TXA2',
                  'STAT5', 'MTORC1', 'UV', 'TNFA', 'NFKB')
  plot_dat$pathway <- sapply(plot_dat$pathway, function(x) abbrev_pathway_name(x, gtf, gene_names))
  plot_dat <- as.data.frame(plot_dat)
  ### Order pathways
  #a<-hclust(dist(plot_dat$NES))
  
  pathways <- unique(plot_dat$pathway)
  pathways <- pathways[order(sapply(pathways, function(pathway) mean(abs(plot_dat[which(plot_dat$pathway == pathway), 'NES']))), decreasing = T)]
  
  plot_dat$pathway <- factor(plot_dat$pathway, levels = pathways)

  legend_size <-7
  x_size <- 10
  y_size <- 9
  ### If any test values are combinations of multiple groups (has '+')
  if(any(grepl('\\+', plot_dat$fgsea_res))){
    ### Check if they are equal to all the groups without a +
    withplus <- unique(plot_dat$fgsea_res)[grepl('\\+', unique(plot_dat$fgsea_res))]
    sansplus <- unique(plot_dat$fgsea_res)[!grepl('\\+', unique(plot_dat$fgsea_res))]
    withplus <- lapply(withplus, function(x) strsplit(x, '-')[[1]])
    ### Just value1
    sansplus1 <- lapply(sansplus, function(x) strsplit(x, '-')[[1]][1])
    sansplus1 <- unique(unlist(sansplus))
    
    sets1 <-lapply(withplus, function(x) strsplit(x[1], '\\+')[[1]])
    sets2 <-lapply(withplus, function(x) strsplit(x[2], '\\+')[[1]])
    for(ii in 1:length(sets1)){
      a <- sansplus1[which(!sansplus1 %in% sets1[[ii]])]
      if(all(union(sets2[[ii]], a) == intersect(sets2[[ii]],a))){
        plot_dat[which(plot_dat$fgsea_res == paste0(withplus[[ii]][1], '-', withplus[[ii]][2])), 'fgsea_res'] <-
          paste0(withplus[[ii]][1], '-All')
      }
    }
  }
  
  p <- ggplot(data = plot_dat, aes(y = pathway, x = fgsea_res, size = -log10(padj), fill = NES)) +
    ggtitle(plot_title) + 
    labs(y = y_title) + 
    geom_point(shape = 21) + 
    scale_size_continuous(range = c(0,4)) +
    scale_fill_gradientn(colours = viridis(100)) +
    fig_theme +
    scale_x_discrete(position = "top") + 
    scale_y_discrete(position = "right") + 
    theme(axis.title.x = element_blank()) +
    labs(size = '-log10 adj. p value', fill = 'Normalized\nEnrichment Score') +
    theme(legend.position = "bottom",
          legend.title = element_text(colour = "black"),#, size = legend_size),
          legend.spacing = unit(0, 'points'), legend.key.size = unit(8, 'points'), 
          legend.direction = 'horizontal', legend.box = 'horizontal',
          legend.box.spacing = unit(1, 'points'),
          legend.box.margin = margin(t = 0, r = -10, b = 0, l = 10, unit = "point")) +
    
    guides(shape = guide_legend(override.aes = list(size = 8))) + 
    guides(color = guide_legend(override.aes = list(size = 8))) + 
    theme(panel.spacing = unit(1, "lines")) +
    guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5), size = guide_legend(title.position = "top", title.hjust = 0.5))
  
  print(p)
} 

dev.off()
