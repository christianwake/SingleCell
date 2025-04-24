
library('readr')
library('RColorBrewer')
library('cowplot')
library('grid')
library('scales')
library('ggnewscale')
library('sys')
library('stringr')
library('Seurat')
library('ggplot2')
library('dplyr')
library('viridis')
library('data.table')
library('patchwork')
#library('PKI')
library('tinytex')
library('dsb')
#library('tidyverse')
#library('pastecs')
library('reticulate')
#library('umap')
library('gridExtra')
#library('logspline')
library('readxl')
library('WriteXLS')
library('Matrix')
library('pastecs')

source('/data/vrc_his/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/sample_sheet_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/sc_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/dehash_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/Utility_functions.R')

if(interactive()){
  project <- '2021614_21-002'

  runs_dir <- '/data/vrc_his/douek_lab/Runs/'
  covs_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, 
                      '/Sample_sheet_Innate.csv')
  out_pdf1 <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/data/Dehash_comparisons_calls2.pdf')
  out_pdf2 <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/data/Dehash_comparisons_threshs2.pdf')
  out_threshs <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/data/Dehash_comparisons2.csv')
  demux_methods <- 'MULTIseqDemux,Trough,Custom'
  ### Method 1 of input: specifying Trough. But it won't be exactly reproducible due to randomness.
  do_trough <- F
  ### Method 2 of input: a csv of threshold values
  thresh_files <- sapply(strsplit(demux_methods, ',')[[1]], function(x) 
    paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/data/Dehash_threshs_', x, '.csv'))
  ### Method 3 of input: csv of cell dehashing calls
  calls_files <- sapply(strsplit(demux_methods, ',')[[1]], function(x) 
    paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/data/Dehash_calls_', x, '.tsv'))
}else{
  args <- commandArgs(trailingOnly=TRUE)

  runs_dir <- args[1]
  project <- args[2]
  demux_methods <- args[3]
  out_pdf1 <- args[4]
  out_pdf2 <- args[5]
  out_threshs <- args[6]
  covs_file <- args[7]
  the_rest <- args[8:length(args)]
  calls_files <- the_rest[grepl('Dehash_calls', the_rest)]
  thresh_files <- the_rest[grepl('Dehash_threshs', the_rest)]
  
  do_trough <- F
}

if(project == '2023105_21-0012'){
  CR_version = '1.1.1'
} else{
  CR_version = '7.1.1'
}
project_dir <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project)

### Read covariates
dat <- read.csv(covs_file, stringsAsFactors = F, header = T,
                check.names = F, 
                colClasses = 'character',
                sep = ',')

row.names(dat) <- dat$Sample_ID
##### Standardize protein names
### space then parentheses <- parentheses
dat$HTO_index_name <- gsub(' \\(', '\\(', dat$HTO_index_name)
### spaces, underscores to '.'
dat$HTO_index_name <- make.names(dat$HTO_index_name, allow_ = 'F')


method_options <- c('MULTIseqDemux', 'Trough', 'Custom')
demux_methods <- strsplit(demux_methods, ',')[[1]]
### Reduce inputs to those in options
demux_methods <- demux_methods[which(toupper(demux_methods) %in% toupper(method_options))]
if(length(demux_methods) < 1){
  warning(paste0('Input demultiplex method(s) are not among the options (', toString(method_options), ')'))
}
### Convert inputs to the exact case of the method_options
demux_methods <- sapply(demux_methods, function(x) method_options[which(toupper(method_options) == toupper(x))])

# names(thresh_files) <- demux_methods
# names(calls_files) <- demux_methods

names(calls_files) <- sapply(calls_files, function(x) gsub('Dehash_calls_', '', gsub('.tsv', '', basename(x))))
names(thresh_files) <- sapply(thresh_files, function(x) gsub('Dehash_threshs_', '', gsub('.csv', '', basename(x))))

# if(do_trough == T){
#   demux_methods <- c('Trough')
#   calls_files <- calls_files[which(!names(calls_files) == 'Trough')]
# } else{
#   demux_methods <- c()
# }


read_thresh <- function(thresh_files, demux_method){
  dat <- read.table(thresh_files[demux_method], sep = ',', header = F, quote = '')
  colnames(dat) <- c('CR_ID', 'hash_id', demux_method)
  return(dat)
}

### If the file is input, read list of threshs
if(length(thresh_files) > 0){
  if('Custom' %in% names(thresh_files)){
    thresh_dat <- read.table(thresh_files['Custom'], sep = ',', header = T, quote = '')
  }else{
    if(length(thresh_files) == 1){
      thresh_dat <- read_thresh(thresh_files, names(thresh_files))
    } else {
      ### Read dehash call files
      thresh_list <- lapply(names(thresh_files), function(x) read_thresh(thresh_files, x))
      ### Merge all into 1
      thresh_dat <- thresh_list[[1]]
      for(i in 2:length(thresh_list)){
        thresh_dat <- merge(thresh_dat, thresh_list[[i]], by = c('CR_ID', 'hash_id'))
      }
    }
  }

  ### To use the threshold csv file it must have info for all the hashes in the sample sheet
  if(all(paste0(dat$CR_ID, '_', dat$HTO_index_name) %in% paste0(thresh_dat$CR_ID, '_', thresh_dat$hash_id))){
    ### Add them to demux_methods if they aren't already
    newcols <- colnames(thresh_dat)[which(!colnames(thresh_dat) %in% c('CR_ID', 'hash_id', demux_methods))]
    demux_methods <- c(demux_methods, newcols)
  }
} else{
  thresh_dat <- NA
}
### Read dehash call files
labels_list <- lapply(calls_files, function(x) read.table(x, sep = '\t', header = T, quote = ''))
demux_methods <- unique(c(demux_methods, names(labels_list)))

pdf(out_pdf1)
for(i in names(labels_list)){
  png_file <- gsub('\\.pdf', paste0('_', i, '.png'), out_pdf1)
  labels <- labels_list[[i]]
  freqtable <- table(labels$Assignment_simple, labels$Assignment_CR)
  ### log10(RNA count) x log10(protein count) by Dehash call, split by CellRanger call
  p3 <- ggplot(labels, aes(x = rna_size, y = prot_size, col = Assignment_simple)) + 
    labs(title = 'All data (by cellranger call)', y='log10(protein count)', x = 'log10(RNA count)', color = paste0(i, '\nDehash call')) +
    geom_point(size =0.5 ) +
    facet_wrap(~Assignment_CR) +
    annotation_custom(tableGrob(freqtable), xmin = 1, xmax = 4, ymin = 0, ymax = 2)
  png(png_file)
  print(p3)
  dev.off()
  print(p3)
}
dev.off()

line_options <- c('dashed', 'dotted', 'dotdash', 'twodash', 'longdash', 
                  'dashed', 'dotted', 'dotdash', 'twodash', 'longdash')
line_types <- line_options[1:length(demux_methods)]

### Define a color for each dehash
plot_colors <- brewer.pal(n = length(demux_methods), 'Set2')
names(plot_colors) <- demux_methods

### Remove NAs
plot_colors <- plot_colors[which(!is.na(names(plot_colors)))]

csv_out_dat <- as.data.frame(matrix(nrow = 0, ncol = length(demux_methods) + 2))
colnames(csv_out_dat) <- c('CR_ID', 'hash_id', demux_methods)

pdf(out_pdf2)
### For each cellranger count output, read and demultiplex
for(id in unique(dat$CR_ID)){
#for(id in c('CR1')){
  print(id)
  
  hashtags <- as.data.frame(dat[which(dat$CR_ID == id),])[,'HTO_index_name']
  ### Sample names rather than CellRanger IDs
  specimens <- as.data.frame(dat[which(dat$CR_ID == id),])[, 'Sample_ID']
  names(specimens) <- hashtags
  specimens <- c(specimens, 'Negative', 'Doublet')
  names(specimens) <- c(hashtags, 'Negative', 'Doublet')
  ### Get desired directory of data
  data_dir <- my_read_10x(id, dat, runs_dir, project_dir, cellranger = c('multi', 'count'), filtered = 'filtered', CR_version)
  ### Read CellRanger-filtered data
  CR_filtered <- Read10X(data.dir = data_dir)
  stained_cells <- colnames(CR_filtered$`Gene Expression`)
  ### Account for weird "multi_" prefacing cell_id if processed by cellranger multi
  stained_cells <- sapply(stained_cells, function(cn) strsplit(cn, '_')[[1]][2])
  ### Remove the trailing "-1"
  stained_cells <- gsub('-1$', '', stained_cells)
  
  rm(CR_filtered)
  ### Get desired directory of data
  data_dir <- my_read_10x(id, dat, runs_dir, project_dir, cellranger = c('multi', 'count'), filtered = 'raw', CR_version)
  CR_raw <- Read10X(data.dir = data_dir)
  
  prot <- CR_raw[[names(CR_raw)[which(names(CR_raw) != 'Gene Expression')]]]
  ### Account for weird "multi_" prefacing cell_id if processed by cellranger multi
  colnames(prot) <- sapply(colnames(prot), function(cn) strsplit(cn, '_')[[1]][2])
  ### Remove the trailing "-1"
  colnames(prot) <- gsub('-1$', '', colnames(prot))
  
  #row.names(prot) <- gsub('_totalseq', '', row.names(prot))
  ### Account for the '.1' appended to protein names if they are also a gene name.....
  if(any(grepl("\\.1$", row.names(prot)))){
    sub <- row.names(prot)[grepl("\\.1$", row.names(prot))]
    key <- sapply(sub, function(x) gsub('\\.1$', '', x))
    names(key) <- sub
    ### Check that the ".1" is there because the version without is in RNA
    key <- key[(key %in% row.names(CR_raw[['Gene Expression']])) & !(key %in% row.names(prot))]
    ### Replace values in row.names that match names(key) with key
    row.names(prot)[which(row.names(prot) %in% names(key))] <- key[ row.names(prot)[which(row.names(prot) %in% names(key))]]
  }
  ### space then parentheses -> parentheses
  row.names(prot) <- gsub(' \\(', '\\(', row.names(prot))
  ### spaces, underscores to '.'
  row.names(prot) <- make.names(row.names(prot), allow_ = 'F')
  
  ### Save space by removing it now that we saved the items separately
  rm(CR_raw)
  
  hash <- t(as.data.frame(prot[hashtags, ]))
  ### CLR transformation
  datNorm <- as.data.frame(apply(hash, 2, clr))
  ### Use density with bandwidth (max - min)/50
  densities <- lapply(1:length(colnames(datNorm)), function(i) get_density(datNorm[,i], 75))
  names(densities) <- colnames(datNorm)
  
  ### Subset the labels to the CR_ID
  labels <- lapply(labels_list, function(x) x[which(x$CR_ID == id),])
  #hash_ids <- unique(unlist(lapply(labels, function(x) unique(x$Hashtag))))
  hash_ids <- unique(dat[which(dat$CR_ID == id), 'HTO_index_name'])
  hash_ids <- hash_ids[which(!hash_ids %in% c('Negative','Doublet', 'Ambiguous'))]
  ### Order them
  hash_ids <- hash_ids[order(parse_number(hash_ids))]
  ### Plot list
  plot_list <- vector(mode = 'list', length = length(hash_ids))
  names(plot_list) <- hash_ids
  maxy <- 0
  miny <- 1
  maxx <- 0
  minx <- 1
  for(i in 1:length(hash_ids)){
    hash_id <- hash_ids[i]
    print(hash_id)
 
    ### Determine the threshold used by the input dehash results
    #threshs <- sapply(labels, function(lab) determine_threshold(lab, datNorm, hash_id))
    threshs <- sapply(demux_methods, function(demux_method) determine_or_fetch_threshold(
                       demux_methods, demux_method, thresh_dat, labels, datNorm, id, hash_id)) 
    if(length(threshs[1]) > 1){
      warning(paste0('More than one threshold listed for ', id, ', ', hash_id))
    }
    names(threshs) <- demux_methods
    print(threshs)
    #print(N_thresholds[hash_id])
    
    thresh_plot_dat <- data.frame(xintercept = threshs, Method_threshold = names(threshs), color = plot_colors, stringsAsFactors = F)
    plot_dat <- data.frame(size = densities[[hash_id]]$x, density = densities[[hash_id]]$y)
    ### cells x the 1 HTO (counts)
    cr_dat <- datNorm[, hash_id, drop = F]
    ### Add column for boolean presence in stained_cells or not
    cr_dat$CR <- ifelse(row.names(cr_dat) %in% stained_cells, '+', '-')
    ### Add column with bin information based on the first column and bins based on the plots x axis
    cr_dat <- cr_dat %>% mutate(new_bin = cut(eval(parse(text = hash_id)), breaks = plot_dat$size))
    ### table of - and + per bin
    freqs <- as.data.frame.matrix(table(cr_dat$new_bin, cr_dat$CR))
    ### Cumulative + as bins progress rightward, or - leftward
    freqs$cumul_pos <- cumsum(freqs[, '+'])
    freqs$cumul_neg <- rev(cumsum(rev(freqs[, '-'])))
    ### Cumulative CR+ (append last value in order to match lengths)
    plot_dat[, 'cumulative_CR'] <- c(freqs[,'cumul_neg'], freqs[length(freqs$cumul_neg),'cumul_neg'])
    max_cr <- max(plot_dat[, 'cumulative_CR'])
    #which(plot_dat[, 'cumulative_CR'] > 2^14)
    miny <- min(miny, min(plot_dat$density))
    maxy <- max(maxy, max(plot_dat$density))
    minx <- min(minx, min(plot_dat$size))
    maxx <- max(maxx, max(plot_dat$size))
    ### Plot without thresholds yet
    print(max_cr)
    p <- ggplot(plot_dat) + 
         ggtitle(hash_id) + theme(plot.title = element_text(size = 10, face = "bold"), 
                            panel.background = element_rect(fill = "lightgrey", colour = "darkgrey",
                            size = 0.5, linetype = "solid")) + 
     aes(x = size, y = density) + 
     geom_line(aes(color = cumulative_CR), size = 1.5) +
     scale_colour_gradient2(name = 'cumulative\nCR- calls', trans = 'log2',
                            limits = c(1, max_cr/1000),
                            midpoint = log2(max_cr/5000),
                            high = 'red', low = 'green', mid = 'white', na.value = 'black') +
     scale_y_continuous(trans = 'log2', labels = scientific) + 
     new_scale_color() +
     geom_vline(aes(xintercept = xintercept, color = Method_threshold), thresh_plot_dat, size = 0.5, linetype = line_types) +
      theme(axis.title.y=element_blank()) + labs(color='Dehash method\nMin threshold')
    ### For those on the right-hand side (even), turn off the y axis tick labels
    if((i %% 2) == 0){
      p <- p + theme(axis.text.y=element_blank())
    }
    ### For those not on the bottom row, turn off the x axis label
    if(!(i %in% c(length(hash_ids), (length(hash_ids) - 1)))){
      p <- p + theme(axis.title.x=element_blank(), axis.text.x = element_blank())
    } else{
      p <- p + xlab('size')
    }
    plot_list[[hash_id]] <- p

    subdat <- as.data.frame(matrix(nrow = 0, ncol = length(threshs) + 2))
    subdat[1, ] <- c(id, hash_id, threshs)
    colnames(subdat) <- c('CR_ID', 'hash_id',  names(threshs))
    csv_out_dat <- rbind(csv_out_dat, subdat)
  }
  ### Cowplot
  legend <- get_legend(plot_list[[1]] + theme(legend.box.margin = margin(0, 0, 0, 12)))
  ## Turn off legends
  plot_list <- lapply(plot_list, function(x)  x + theme(legend.position = "none"))
  ### Set maxy upper bound
  maxy <- min(2.5, maxy) 
  ### Set miny lower bound
  miny <- max(0.0000025)
  ### Set axis ranges uniformly
  plot_list <- lapply(plot_list, function(x)  x + scale_y_continuous(trans = 'log2', labels = scientific, limits = c(miny, maxy)))
  plot_list <- lapply(plot_list, function(x)  x + xlim(minx, maxx))
  
  rel_heights <- c(rep(1, length(plot_list)/2-1), 1.3)
  p0 <- plot_grid(plotlist = plot_list, ncol = 2, rel_heights = rel_heights)
  p1 <- plot_grid(p0, legend, rel_widths = c(3, .4))
  plot_title <- ggdraw() + draw_label(paste0(id, ' hash densities'))
  pfinal <- plot_grid(plot_title, p1, rel_heights = c(0.1, 1), ncol = 1)
  print(pfinal)
}

write.table(csv_out_dat, out_threshs, quote = F, sep = ',', row.names = F, col.names = T)
dev.off()

### Incomplete. Some lines of code I removed from an early version of Dehash.R
read_CR_metrics <- function(quant_func){
  ## Combine cellranger count output metrics
  if(quant_func == 'multi'){
    runs_subdir <- '/outs/multi/count/raw_feature_bc_matrix/'
    metrics <- lapply(unique(dat$CR_ID), function(x) as.data.frame(read.table(paste0(runs_dir, x, '/outs/per_sample_outs/', x, '/metrics_summary.csv'), sep = ',')[2,]))
    metrics <- rbindlist(metrics)
    colnames(metrics) <- gsub(' ', '_', unlist(read.table(paste0(runs_dir, dat$CR_ID[1], '/outs/per_sample_outs/', dat$CR_ID[1], '/metrics_summary.csv'), sep = ',')[1,]))
  } else if(quant_func == 'count'){
    runs_subdir <- '/outs/raw_feature_bc_matrix/'
    metrics <- lapply(unique(dat$CR_ID), function(x) as.data.frame(read.table(paste0(runs_dir, x, '/outs/metrics_summary.csv'), sep = ',')[2,]))
    metrics <- rbindlist(metrics)
    colnames(metrics) <- gsub(' ', '_', unlist(read.table(paste0(runs_dir, dat$CR_ID[1], '/outs/metrics_summary.csv'), sep = ',')[1,]))
  }
  metrics <- as.data.frame(t(metrics))
  colnames(metrics) <- unique(dat$CR_ID)
  write.table(metrics, out_metrics, sep = '\t', quote = F,  row.names = T, col.names = T)
}

