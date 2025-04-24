
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
  qc_name <- '2024-07-31'
  raw_or_norm <- 'raw'
  probes <- 'S2P_WT_APC,RBD_WT_APC,NTD_WT_APC,gp120_APC,S2P_Beta_PE,RBD_Beta_PE,NTD_Beta_PE,gp120_PE'
  if(raw_or_norm == 'raw'){
    rds_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, 
                       '/results/', qc_name, '/All_data.RDS')
    do_norm <- T
  } else{
    rds_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, 
                       '/results/', qc_name, '/Combined_batches.RDS')
    do_norm <- F
  }
  runs_dir <- '/data/vrc_his/douek_lab/Runs/'
  covs_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name,  
                      '/Sample_sheet.csv')
  out_pdf <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name,
                     '/Libraseq/', raw_or_norm, '/Libraseq_comparisons_threshs.pdf')
  out_threshs <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name,
                        '/Libraseq/', raw_or_norm , '/Libraseq_comparisons.csv')
  demux_methods <- 'MULTIseqDemux,Trough'
  ### Method 1 of input: specifying Trough. But it won't be exactly reproducible due to randomness.
  do_trough <- F
  ### Method 2 of input: a csv of threshold values
  thresh_files <- sapply(strsplit(demux_methods, ',')[[1]], function(x) 
    paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name,
           '/Libraseq/', raw_or_norm, '/thresholds_', x, '.csv'))
  ### Method 3 of input: csv of cell dehashing calls
  calls_files <- sapply(strsplit(demux_methods, ',')[[1]], function(x) 
    paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project,  '/results/', qc_name,
           '/Libraseq/', raw_or_norm, '/calls_', x, '.tsv'))
}else{
  args <- commandArgs(trailingOnly=TRUE)
  runs_dir <- args[1]
  project <- args[2]
  demux_methods <- args[3]
  probes <- args[4]
  do_norm <- as.logical(args[5])
  out_pdf <- args[6]
  out_threshs <- args[7]
  covs_file <- args[8]
  rds_file <- args[9]
  the_rest <- args[10:length(args)]
  calls_files <- the_rest[grepl('calls', the_rest)]
  thresh_files <- the_rest[grepl('thresholds', the_rest)]
  
  do_trough <- F
}

project_dir <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project)

### Read covariates
dat <- read.csv(covs_file, stringsAsFactors = F, header = T,
                check.names = F, 
                colClasses = 'character',
                sep = ',')

row.names(dat) <- dat$Sample_ID

probes <- strsplit(probes, ',')[[1]]
##### Standardize protein names
### space then parentheses -> parentheses
probes <- gsub(' \\(', '\\(', probes)
### spaces, underscores to '.'
probes <- make.names(probes, allow_ = 'F')

downsample_var <- 0.3
if(!interactive()){
  ### Read Seurat data
  sdat <- readRDS(rds_file)
  print('Done reading Seurat object')
  ### Downsample so I can work interactiveley for testing
  assays <- c('RNA', 'prot')
  assays <- assays[which(assays %in% names(sdat@assays))]
  layers <- c('counts', 'data')
  layers <- layers[which(layers %in% names(sdat@assays$RNA@layers))]
  ssub <- DietSeurat(sdat, layers = layers, assays = assays)
  cells <- sample(x = colnames(ssub), size = (length(colnames(ssub)) * downsample_var), replace = F)
  
  ssub <- subset(ssub, cells = cells)
  saveRDS(ssub, file = gsub('.RDS', paste0('_DownSampledTo', downsample_var, '.RDS'), rds_file))
  print('saved downsampled version for testing')
} else{
  sdat <- readRDS(gsub('.RDS', paste0('_DownSampledTo', downsample_var, '.RDS'), rds_file))
}

method_options <- c('MULTIseqDemux', 'Trough', 'Custom')
demux_methods <- strsplit(demux_methods, ',')[[1]]
### Reduce inputs to those in options
demux_methods <- demux_methods[which(toupper(demux_methods) %in% toupper(method_options))]
if(length(demux_methods) < 1){
  warning(paste0('Input demultiplex method(s) are not among the options (', toString(method_options), ')'))
}
### Convert inputs to the exact case of the method_options
demux_methods <- sapply(demux_methods, function(x) method_options[which(toupper(method_options) == toupper(x))])

names(calls_files) <- sapply(calls_files, function(x) gsub('calls_', '', gsub('.tsv', '', basename(x))))
names(thresh_files) <- sapply(thresh_files, function(x) gsub('thresholds_', '', gsub('.csv', '', basename(x))))

read_thresh <- function(thresh_files, demux_method){
  dat <- read.table(thresh_files[demux_method], sep = ',', header = F, quote = '')
  colnames(dat) <- c('CR_ID', 'probe', demux_method)
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
        thresh_dat <- merge(thresh_dat, thresh_list[[i]], by = c('CR_ID', 'probe'))
      }
    }
  }
} else{
  thresh_dat <- NA
}
### Read dehash call files
labels_list <- lapply(calls_files, function(x) read.table(x, sep = '\t', header = T, quote = ''))
demux_methods <- unique(c(demux_methods, names(labels_list)))

line_options <- c('dashed', 'dotted', 'dotdash', 'twodash', 'longdash', 
                  'dashed', 'dotted', 'dotdash', 'twodash', 'longdash')
line_types <- line_options[1:length(demux_methods)]

### Define a color for each dehash
plot_colors <- brewer.pal(n = length(demux_methods), 'Set2')
names(plot_colors) <- demux_methods

### Remove NAs
plot_colors <- plot_colors[which(!is.na(names(plot_colors)))]

csv_out_dat <- as.data.frame(matrix(nrow = 0, ncol = length(demux_methods) + 2))
colnames(csv_out_dat) <- c('CR_ID', 'probe', demux_methods)

pdf(out_pdf)
### For each cellranger count output, read and demultiplex
for(id in unique(dat$CR_ID)){
  print(id)
  ### Subset to those of this CR_ID
  sdatsub <- subset(sdat, subset = CR_ID == id)
  ### Summarize data types in this CR_ID for plot subtitle
  subtitle <- paste(
  paste0('Arms: ', toString(names(table(sdatsub@meta.data$Arm)))),
  paste0('Visits: ', toString(names(table(sdatsub@meta.data$Visit)))),
  toString(names(table(sdatsub@meta.data$Cell_subset))),
  sep = '. ')

  ### Define the data to use
  if(do_norm){
    ### CLR transformation
    print('Doing CLR transformation')
    datNorm <- as.data.frame(apply(sdatsub@assays$prot@data[probes, ], 1, clr))
    xlabel <- 'Tranformed counts'
  } else{
    datNorm <- as.data.frame(t(as.data.frame(sdatsub@assays$prot@data[probes,])))
    xlabel <- 'Normalized counts'
  }

  ### Use density with bandwidth (max - min)/50
  ### Names should be probe IDs, so datNorm should be cells x probes
  densities <- lapply(1:length(colnames(datNorm)), function(i) get_density(datNorm[,i], 75))
  names(densities) <- colnames(datNorm)
  
  ### Subset the labels to the CR_ID
  labels <- lapply(labels_list, function(x) x[which(x$CR_ID == id),])
  
  ### Plot list
  plot_list <- vector(mode = 'list', length = length(probes))
  names(plot_list) <- probes
  maxy <- 0
  miny <- 1
  maxx <- 0
  minx <- 1
  for(i in 1:length(probes)){
    probe <- probes[i]
    print(probe)
    
    ### Determine the threshold used by the input dehash results
    #threshs <- sapply(labels, function(lab) determine_threshold(lab, datNorm, probe))
    threshs <- sapply(demux_methods, function(demux_method) determine_or_fetch_threshold(
      demux_methods, demux_method, thresh_dat, labels, datNorm, id, probe, hash_name = 'probe')) 
    if(length(threshs[1]) > 1){
      warning(paste0('More than one threshold listed for ', id, ', ', probe))
    }
    names(threshs) <- demux_methods
    print(threshs)
    #print(N_thresholds[probe])
    
    thresh_plot_dat <- data.frame(xintercept = threshs, Method_threshold = names(threshs), color = plot_colors, stringsAsFactors = F)
    plot_dat <- data.frame(size = densities[[probe]]$x, density = densities[[probe]]$y)
    ### cells x the 1 HTO (counts)
    cr_dat <- datNorm[, probe, drop = F]

    miny <- min(miny, min(plot_dat$density))
    maxy <- max(maxy, max(plot_dat$density))
    minx <- min(minx, min(plot_dat$size))
    maxx <- max(maxx, max(plot_dat$size))
    ### Plot without thresholds yet
    p <- ggplot(plot_dat) + 
      ggtitle(probe) + theme(plot.title = element_text(size = 10, face = "bold"), 
                               panel.background = element_rect(fill = "lightgrey", colour = "darkgrey",
                                                               size = 0.5, linetype = "solid")) + 
      aes(x = size, y = density) + 
      geom_line(size = 1) +
      #scale_y_continuous(trans = 'log2', labels = scientific) + 
      new_scale_color() +
      geom_vline(aes(xintercept = xintercept, color = Method_threshold), thresh_plot_dat, size = 0.5, linetype = line_types) +
      theme(axis.title.y=element_blank()) + labs(color='Dehash method\nMin threshold')
    ### For those on the right-hand side (even), turn off the y axis tick labels
    if((i %% 2) == 0){
      p <- p + theme(axis.text.y=element_blank())
    }
    ### For those not on the bottom row, turn off the x axis label
    if(!(i %in% c(length(probes), (length(probes) - 1)))){
      p <- p + theme(axis.title.x=element_blank(), axis.text.x = element_blank())
    } else{
      p <- p + xlab(xlabel)
    }
    plot_list[[probe]] <- p
    
    ### Will contain only threshold info
    CR_thresh <- as.data.frame(matrix(nrow = 0, ncol = length(threshs) + 2))
    CR_thresh[1, ] <- c(id, probe, threshs)
    colnames(CR_thresh) <- c('CR_ID', 'probe',  names(threshs))
    csv_out_dat <- rbind(csv_out_dat, CR_thresh)
  }
  ### Cowplot
  legend <- get_legend(plot_list[[1]] + theme(legend.box.margin = margin(0, 0, 0, 12)))
  ## Turn off legends
  plot_list <- lapply(plot_list, function(x)  x + theme(legend.position = "none"))
  ### Set maxy upper bound
  maxy <- min(2.5, maxy) 
  ### Set miny lower bound
  if(miny == 0){
    miny = -0.0000001
  }
  ### Set axis ranges uniformly
  plot_list <- lapply(plot_list, function(x)  x + scale_y_continuous(trans = 'log2', 
                                                              labels = scientific, limits = c(miny, maxy)))
  plot_list <- lapply(plot_list, function(x)  x + xlim(minx, maxx))
  
  rel_heights <- c(rep(1, length(plot_list)/2-1), 1.3)
  p0 <- plot_grid(plotlist = plot_list, ncol = 2, rel_heights = rel_heights)
  p1 <- plot_grid(p0, legend, rel_widths = c(3, .4))
  plot_title <- ggdraw() + draw_label(paste0(id, ' hash densities\n', subtitle))
  #subtitle <- ggdraw() +
   # draw_label(subtitle, x = 0.05, hjust = 0, vjust = 1)
  pfinal <- plot_grid(plot_title, p1, rel_heights = c(0.1, 1), ncol = 1)
  print(pfinal)
}

write.table(csv_out_dat, out_threshs, quote = F, sep = ',', row.names = F, col.names = T)
dev.off()


