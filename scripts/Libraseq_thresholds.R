library('sys')
library('Seurat')
library('ggplot2')
library('dplyr')
library('viridis')
library('data.table')
library('patchwork')
#library('PKI')
library('tinytex')
library('reticulate')
library('gridExtra')
library('cowplot')
library('readxl')
library('Matrix')
library('pastecs')
library('RColorBrewer')
library('ggnewscale')

source('/data/vrc_his/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/sample_sheet_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/sc_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/dehash_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/Utility_functions.R')

if(interactive()){
  project <- '2021614_21-002'
  qc_name <- '2024-07-31'
  runs_dir <- '/data/vrc_his/douek_lab/Runs/'
  demux_method <- 'MULTIseqDemux'
  demux_method <- 'Trough'
  raw_or_norm <- 'norm'
  if(raw_or_norm == 'raw'){
    rds_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, 
                     '/results/', qc_name, '/All_data.RDS')
    do_norm <- T
  } else{
    rds_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, 
                     '/results/', qc_name, '/Combined_batches.RDS')
    do_norm <- F
  }
  
  covs_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, 
                      '/results/', qc_name, '/Sample_sheet.csv')
  thresh_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name,
                        '/Libraseq/', raw_or_norm, '/thresholds_', demux_method, '.csv')
  pdf_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, 
                     '/Libraseq/', raw_or_norm, '/thresholds_', demux_method, '.pdf')
  gtf_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/data/gtf.RDS')
  probes <- 'S2P_WT_APC,RBD_WT_APC,NTD_WT_APC,gp120_APC,S2P_Beta_PE,RBD_Beta_PE,NTD_Beta_PE,gp120_PE'
}else{
  args <- commandArgs(trailingOnly=TRUE)
  
  covs_file <- args[1]
  gtf_file <- args[2]
  rds_file <- args[3]
  runs_dir <- args[4]
  project <- args[5]
  do_norm <- as.logical(args[6])
  demux_method <- as.character(args[7])
  thresh_file <- args[8]
  pdf_file <- args[9]
}
probes <- 'S2P_WT_APC,RBD_WT_APC,NTD_WT_APC,gp120_APC,S2P_Beta_PE,RBD_Beta_PE,NTD_Beta_PE,gp120_PE'
Libraseq_CR_ID <- function(dat, id, probes, project_dir, runs_dir, rds_file = NA, sdat = NA, demux_method = 'Trough', 
                           do_norm = T, autoThresh = T, gtf = NA, thresh_file = NA, thresh_dat = NA, pdf_file = NA){
  print(paste0('Probes: ', toString(probes)))
  probe_name <- 'probe_id'
  ### Sample names rather than CellRanger IDs
  specimens <- as.data.frame(dat[which(dat$CR_ID == id),])[, 'Sample_ID']
  names(specimens) <- probes
  specimens <- c(specimens, 'Negative', 'Doublet', 'Ambiguous')
  names(specimens) <- c(probes, 'Negative', 'Doublet', 'Ambiguous')
  
  sdat_input <- suppressWarnings(is.na(sdat))
  if(is.na(rds_file)){
    print(paste0('No rds file. Reading data from', runs_dir))
    ### Get desired directory of data
    data_dir <- my_read_10x(id, dat, runs_dir, project_dir, cellranger = c('multi', 'count'), 
                            filtered = 'raw')
    print(data_dir)
    CR_raw <- Read10X(data.dir = data_dir)
    
    ### Get some RNA metadata
    rna <- CR_raw$`Gene Expression`
    ### Account for weird "multi_" prefacing cell_id if processed by cellranger multi
    colnames(rna) <- sapply(colnames(rna), function(cn) strsplit(cn, '_')[[1]][2])
    ### Remove the trailing "-1"
    colnames(rna) <- gsub('-1$', '', colnames(rna))
    
    prot <- CR_raw[[names(CR_raw)[which(names(CR_raw) != 'Gene Expression')]]]
    ### Account for weird "multi_" prefacing cell_id if processed by cellranger multi
    colnames(prot) <- sapply(colnames(prot), function(cn) strsplit(cn, '_')[[1]][2])
    ### Remove the trailing "-1"
    colnames(prot) <- gsub('-1$', '', colnames(prot))
    
    ##### protein name standardization
    ### Account for the '.1' appended to protein names if they are also a gene name.....
    if(any(grepl("\\.1$", row.names(prot)))){
      sub <- row.names(prot)[grepl("\\.1$", row.names(prot))]
      key <- sapply(sub, function(x) gsub('\\.1$', '', x))
      names(key) <- sub
      ### Check that the ".1" is there because the version without is in RNA
      key <- key[(key %in% row.names(rna)) & !(key %in% row.names(prot))]
      ### Replace values in row.names that match names(key) with key
      row.names(prot)[which(row.names(prot) %in% names(key))] <- key[ row.names(prot)[which(row.names(prot) %in% names(key))]]
    }
    rm(rna)
    ### space then parentheses -> parentheses
    row.names(prot) <- gsub(' \\(', '\\(', row.names(prot))
    ### spaces, underscores to '.'
    row.names(prot) <- make.names(row.names(prot), allow_ = 'F')
    ### Save space by removing it now that we saved the items separately
    rm(CR_raw)
  } else{
    if(sdat_input){
      print(paste0('Reading data from ', rds_file))
      sdat <- readRDS(rds_file)
    } else{
      print(paste0('Using data already read from  ', rds_file, ' input to function as sdat'))
    }
    ### Subset seurat object, protein data and thresholds to the CR_ID
    sd <- subset(sdat, subset = CR_ID == id)
    prot <- sd@assays$prot@data
  }
  if(!all(probes %in% row.names(prot))){
    warning('Some input probes are not in the protein information')
  }
  
  #prot_size <- log10(Matrix::colSums(prot) + 1)
  print(paste0('Prot size before method selection: ', toString(dim(prot))))
  
  ### Main dehashing wrapper. Output should have columns, CR_ID, hash, and thresh
  thresholds <- Libraseq_method_branch(id, prot, probes, project_dir, demux_method, do_norm, autoThresh, 
                                       thresh_file, thresh_dat, pdf_file)
  print('Thresholds done')
  ### If input method is not in the thresholds (if read from file)
  if(!demux_method %in% colnames(thresholds)){
    demux_method <- colnames(thresholds)[3]
  }
  ##### Density plots with thresholds
  demux_methods <- demux_method
  ### Define line type for each dehash type
  line_options <- c('dashed', 'dotted', 'dotdash', 'twodash', 'longdash', 
                    'dashed', 'dotted', 'dotdash', 'twodash', 'longdash')
  line_types <- line_options[1:length(demux_methods)]
  ### Define a color for each dehash type
  plot_colors <- brewer.pal(n = max(3, length(demux_methods)), 'Set2')
  names(plot_colors) <- demux_methods
  ### Remove NAs
  plot_colors <- plot_colors[which(!is.na(names(plot_colors)))]
  
  protsub <- t(as.data.frame(prot[probes, ]))
  ### CLR transformation
  xlabel <- 'Normalized counts'
  if(do_norm == T){
    print('Doing CLR transformation')
    ### apply over columns (2)
    protsub <- as.data.frame(apply(protsub, 2, clr))
    xlabel <- 'Tranformed counted'
    if(sum(is.na(protsub)) > 0){
      warning('There are NAs in normalized data')
    }
  }
  ### Use kernal density estimation with bandwidth (max - min)/50
  densities <- lapply(1:length(colnames(protsub)), function(i) get_density(protsub[,i], 75))
  names(densities) <- colnames(protsub)
  
  ### Plot list
  probe_ids <- probes
  plot_list <- vector(mode = 'list', length = length(probe_ids))
  names(plot_list) <- probe_ids
  maxy <- 0
  miny <- 1
  maxx <- 0
  minx <- 1

  print('Beginning plot loop')
  for(i in 1:length(probe_ids)){
    probe_id <- probe_ids[i]
    print(probe_id)
    ### Determine the threshold used by the input dehash results
    #threshs <- sapply(labels, function(lab) determine_threshold(lab, protsub, probe_id))
    threshs <- sapply(demux_methods, function(demux_method) 
      thresholds[which(thresholds$CR_ID == id & thresholds[, probe_name] == probe_id), demux_method]) 
    if(length(threshs[1]) > 1){
      warning(paste0('More than one threshold listed for ', id, ', ', probe_id))
    }
    names(threshs) <- demux_methods
    
    thresh_plot_dat <- data.frame(xintercept = threshs, Method_threshold = names(threshs), color = plot_colors, stringsAsFactors = F)
    plot_dat <- data.frame(size = densities[[probe_id]]$x, density = densities[[probe_id]]$y)

    #which(plot_dat[, 'cumulative_CR'] > 2^14)
    miny <- min(miny, min(plot_dat$density))
    maxy <- max(maxy, max(plot_dat$density))
    minx <- min(minx, min(plot_dat$size))
    maxx <- max(maxx, max(plot_dat$size))
    ### Plot without thresholds yet
    p <- ggplot(plot_dat) + 
         ggtitle(probe_id) + 
         theme(plot.title = element_text(size = 10, face = "bold"), 
               panel.background = element_rect(fill = "lightgrey", colour = "darkgrey",
                                  size = 0.5, linetype = "solid")) + 
      aes(x = size, y = density) +
      geom_line(size = 1) +
      #scale_y_continuous(trans = 'log2') + 
      geom_vline(aes(xintercept = xintercept, color = Method_threshold), thresh_plot_dat, size = 0.5, linetype = line_types) +
      labs(color = 'Dehash method\nMin threshold')
    ### For those on the right-hand side (even), turn off the y axis tick labels
    if((i %% 2) == 0){
      p <- p + theme(axis.text.y=element_blank()) + theme(axis.title.y=element_blank())
    }
    ### For those not on the bottom row, turn off the x axis label
    if(!(i %in% c(length(probe_ids), (length(probe_ids) - 1)))){
      p <- p + theme(axis.title.x=element_blank(), axis.text.x = element_blank())
    } else{
      p <- p + xlab(xlabel)
    }
    plot_list[[probe_id]] <- p

    subdat <- as.data.frame(matrix(nrow = 0, ncol = length(threshs) + 2))
    subdat[1, ] <- c(id, probe_id, threshs)
    colnames(subdat) <- c('CR_ID', 'probe_id',  names(threshs))
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
  #miny <- max(0.0000015)
  ### Set axis ranges uniformly
  plot_list <- lapply(plot_list, function(x)  x + scale_y_continuous(trans = 'log2', limits = c(miny, maxy)))
  plot_list <- lapply(plot_list, function(x)  x + xlim(minx, maxx))
  
  rel_heights <- c(rep(1, length(plot_list)/2-1), 1.3)
  p0 <- plot_grid(plotlist = plot_list, ncol = 2, rel_heights = rel_heights)
  p1 <- plot_grid(p0, legend, rel_widths = c(3, .4))
  plot_title <- ggdraw() + draw_label(paste0(id, ' hash densities'))
  pfinal <- plot_grid(plot_title, p1, rel_heights = c(0.1, 1), ncol = 1)
  #print(pfinal)
  out = list(thresholds = thresholds, p = pfinal)
  return(out)
}

### Main demultiplexing function 
### demux_method Should be 'seuratDemux' or 'myDemux', returns dataframe with demultiplexing threshs
# Note that Seurat will give warnings if protein names contain underscores, and will convert them to dashes. 
# This code is meant to handle that conversion, so input underscores are okay. 
# In fact, dashes are a special character in R, so protein names with dashes may cause issues in this code.
Libraseq_method_branch <- function(id, prot, probes, out_dir = NA, demux_method = 'Trough', do_norm,
                                   autoThresh = T, thresh_file = NA, thresh_dat = NA, pdf_file = NA){
  ### Either method will return a 'threshs' dataframe with a row per cell and a column 'Assignment'.
  if(demux_method == 'MULTIseqDemux'){ ### Assignment may be Negative, Doublet, or one of the HTOs
    ### Seurats MULTISeqDemux labeling
    ### Adds protein info and hashtags into seurat object as separate assays with only the 
    ### common cells, normalizes, scales and uses MULTIseqDemux
    ### do_intersect is F because we've already done it
    print('Starting Seurats MULTIseqDemux function')
    #### pdf plots  show the thresholding done by the other (my) method.
    ### Seurat's MULTIseqDemux function
    thresholds <- seurat_demux(id, prot, probes, out = 'thresholds', do_intersect = F, autoThresh = autoThresh, 
                               do_norm, pdf_file, thresh_file = thresh_file, hash_name = 'probe_id')
    print('Completed MULTIseqDemux')
  } else if(demux_method == 'Trough'){ ### Assignment may be Negative, Ambiguous, or one of the HTOs
    ### My demultiplexing and labeling
    print('Determing threshods from trough method')
    thresholds <- trough_dehash(id, prot, hashtags = probes, thresh_file,
                                thresh_N = 2, thresh_n = 2, do_norm, pdf_file, return_labels = F, verbose = F)
    ### Rename columns
    colnames(thresholds) <- c('CR_ID', 'probe_id', 'Trough')
  } else if(demux_method == 'InputThresholdFile'){
    print('Using thresholds from file')
    thresholds <- thresh_dat[which(thresh_dat$CR_ID == id), ]
    ### This function does the calling, not determining of thresholds
    #thresholds <- dehash_by_threshold_input(id, prot, probes, thresh_dat, 
    #                                        hash_colname = 'probe_id', do_norm = do_norm, only1 = F)
  }
  return(thresholds)
}

project_dir <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project)

### Determine whether method is simply reading thresholds from a file and applying them.
demux_method0 <- demux_method
### If thresh_file exists, read it
if(file.exists(thresh_file)){
  print('Threshold file already exists, so will be read and used.')
  thresh_dat <- read.table(thresh_file, sep = ',', header = T, quote = '')
  ### If it seems that there was no header, re read it and assume header names
  if(!any(c('CR_ID', 'probe_id', 'Trough', "MULTIseqDemux") %in% colnames(thresh_dat))){
    thresh_dat <- read.table(thresh_file, sep = ',', header = F, quote = '')
    colnames(thresh_dat) <- c('CR_ID', 'probe_id', demux_method)
  } else{
    ### Allow for the input method to refer to values from the other methods in the file
    demux_methods <- c('Trough', 'MULTIseqDemux', 'Custom')
    demux_methods <- demux_methods[which(demux_methods %in% colnames(thresh_dat))]
    if(!demux_method %in% colnames(thresh_dat)){
      warning('Input demux method is not within the input threshold file columns')
    }
    ### Replace references to other columns with their value in that column
    thresh_dat[, demux_method] <- sapply(1:length(row.names(thresh_dat)), function(i)
      fetch_threshold(demux_methods, demux_method, thresh_dat, 
                      thresh_dat[i, 'CR_ID'], thresh_dat[i, 'probe_id']))
    
    thresh_dat <- thresh_dat[c('CR_ID', 'probe_id', demux_method)]
  }
  demux_method <- 'InputThresholdFile'
} else{
  thresh_dat <- NA
  print('No threshold csv file exists.')
  ### Should not exist if thresh_dat == NA but just in case it does, remove it
  if(!is.na(thresh_file)){
    if(file.exists(thresh_file)){
      file.remove(thresh_file, showWarnings = F)
    }
  }
}
method_options <- c('MULTIseqDemux', 'Trough', 'Custom', 'InputThresholdFile')
if(!demux_method %in% method_options){
  warning(paste0('Input demultiplex method(s) are not among the options (', toString(method_options), ')'))
}
### Convert inputs to the exact case of the method_options
demux_method <- method_options[which(toupper(method_options) == toupper(demux_method))]

print(paste0('Using method ', demux_method))
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

### Do direct matching from gtf
if(grepl('\\.gtf', gtf_file)){
  gtf <- read_gtf(gtf_file, feature_type = 'gene', atts_of_interest = c('gene_id', 'gene_name',  'gene_biotype'))
} else if(grepl('\\.RDS', gtf_file)){
  gtf <- readRDS(gtf_file)
} else{
  warning('No gtf file')
}
row.names(gtf) <- gtf$gene_id

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

print(dim(sdat@assays$prot))
thresh_list <- as.list(rep(NA, length(unique(dat$CR_ID))))
names(thresh_list) <- unique(dat$CR_ID)
### For each cellranger count output, read and demultiplex
print(unique(dat$CR_ID))
pdf(pdf_file)
for(id in unique(dat$CR_ID)){
  print(id)
  pdf_file <- NA
  ### Dehash Wrapper
  out <- Libraseq_CR_ID(dat, id, probes, project_dir, runs_dir, rds_file, sdat, demux_method = demux_method, 
                                     do_norm = do_norm, autoThresh = F, gtf = gtf, thresh_file = thresh_file, 
                                     thresh_dat = thresh_dat, pdf_file = pdf_file)
  thresh_list[[id]] <- out[['thresholds']]
  p <- out[['p']]
  print(p)
}
dev.off()
threshs <- rbindlist(thresh_list)
threshs <- as.data.frame(threshs)
#threshs <- read.table('/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/libraseq.csv', sep = ',')
#colnames(threshs) <- c('CRID', 'x', 'thresh')


