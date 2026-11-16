library('sys')
library('readr')
library('WriteXLS')
library('data.table')

source('/data/vrc_his/douek_lab/snakemakes/Utility_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/DE_functions.R')


### Return the members of list dats that contain gene/geneset gs
get_strats <- function(gs, dats, n = F){
  strats <- c()
  ### For each stratification
  for(i in 1:length(dats)){
    if(gs %in% row.names(dats[[i]])){
      strats <- c(strats, names(dats)[i])
    }
  }
  if(n){
    res <- length(strats)
  } else{
    res <- toString(strats)
  }
  return(res)
}

if(interactive()){
  # args <- c('results/Run2023-05-14/DE/DE.xls', 'results/Run2023-05-14/DE/DE_sig.xls','0.05',
  #           'results/Run2023-05-14/DE/Celltype/All/Celltype-NA-NA_DE_Strat1-All-All-All_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/Time/All/Time-30-31_DE_Strat1-All-All-All_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/Time/Celltype/Time-30-31_DE_Strat1-Celltype-CD4-CD4_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/Time/Celltype/Time-30-31_DE_Strat1-Celltype-CD8-CD8_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/Stim/Celltype/Stim-S205-S42_DE_Strat1-Celltype-CD4-CD4_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/Stim/Celltype/Stim-M37-M45_DE_Strat1-Celltype-CD4-CD4_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/Stim/Celltype/Stim-S302-S42_DE_Strat1-Celltype-CD8-CD8_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/Stim/Celltype/Stim-S302-N27_DE_Strat1-Celltype-CD8-CD8_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/RNA_clusters/All/RNA_clusters-0-1_DE_Strat1-All-All-All_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/RNA_clusters/All/RNA_clusters-2-3_DE_Strat1-All-All-All_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/RNA_clusters/Stim/RNA_clusters-0-3_DE_Strat1-Stim-S42-S42_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-SEB-NS_DE_Strat1-RNA_clusters-0-1_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-SEB-NS_DE_Strat1-RNA_clusters-3-2_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-S42-NS_DE_Strat1-RNA_clusters-0-1_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-S42-NS_DE_Strat1-RNA_clusters-3-2_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-S205-NS_DE_Strat1-RNA_clusters-0-1_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-N88-NS_DE_Strat1-RNA_clusters-0-1_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-M37-NS_DE_Strat1-RNA_clusters-0-1_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-M45-NS_DE_Strat1-RNA_clusters-0-1_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-N33-NS_DE_Strat1-RNA_clusters-0-1_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-S302-NS_DE_Strat1-RNA_clusters-3-2_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-N27-NS_DE_Strat1-RNA_clusters-3-2_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/time/RNA_clusters/time-1-0_DE_Strat1-RNA_clusters-0-0_Strat2-StimCD4-1-1.tsv',  'results/Run2023-05-14/DE/time/RNA_clusters/time-1-0_DE_Strat1-RNA_clusters-3-3_Strat2-StimCD8-1-1.tsv',  'results/Run2023-05-14/DE/time/RNA_clusters/time-2-0_DE_Strat1-RNA_clusters-0-0_Strat2-StimCD4-1-1.tsv',  'results/Run2023-05-14/DE/time/RNA_clusters/time-2-0_DE_Strat1-RNA_clusters-3-3_Strat2-StimCD8-1-1.tsv',  'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-SEB-S42+N27+S302_DE_Strat1-RNA_clusters-3-3_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-SEB-S42+S205+N88+M37+M45+N33_DE_Strat1-RNA_clusters-0-1_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-SEB-S42+S205+N88+M37+M45+N33_DE_Strat1-RNA_clusters-0-0_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-SEB-NS_DE_Strat1-RNA_clusters-0-0_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-SEB-NS_DE_Strat1-RNA_clusters-3-3_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-S42-NS_DE_Strat1-RNA_clusters-0-0_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-S42-NS_DE_Strat1-RNA_clusters-3-3_Strat2-All-All-All.tsv', 'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-S205-NS_DE_Strat1-RNA_clusters-0-0_Strat2-All-All-All.tsv',  'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-N88-NS_DE_Strat1-RNA_clusters-0-0_Strat2-All-All-All.tsv', 'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-M37-NS_DE_Strat1-RNA_clusters-0-0_Strat2-All-All-All.tsv', 'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-M45-NS_DE_Strat1-RNA_clusters-0-0_Strat2-All-All-All.tsv') 
  #                    #'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-N33-NS_DE_Strat1-RNA_clusters-0-0_Strat2-All-All-All.tsv', 'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-S302-NS_DE_Strat1-RNA_clusters-3-3_Strat2-All-All-All.tsv', 'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-N27-NS_DE_Strat1-RNA_clusters-3-3_Strat2-All-All-All.tsv')
  ### CD8
  args <- c('results/Run2023-05-14/DE/DE.xls', 'results/Run2023-05-14/DE/DE_sig.xls','0.05', 'False',
            'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-SEB-NS_DE_Strat1-RNA_clusters-3-2_Strat2-All-All-All.tsv',  
            'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-S42-NS_DE_Strat1-RNA_clusters-3-2_Strat2-All-All-All.tsv',  
            'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-S302-NS_DE_Strat1-RNA_clusters-3-2_Strat2-All-All-All.tsv',  
            'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-N27-NS_DE_Strat1-RNA_clusters-3-2_Strat2-All-All-All.tsv',  
            'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-SEB-S42+N27+S302_DE_Strat1-RNA_clusters-3-3_Strat2-All-All-All.tsv')

  ### CD4
    args <- c('results/Run2023-05-14/DE/DE.xls', 'results/Run2023-05-14/DE/DE_sig.xls','0.05', 'False',
            'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-SEB-NS_DE_Strat1-RNA_clusters-0-1_Strat2-All-All-All.tsv',  
             'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-S42-NS_DE_Strat1-RNA_clusters-0-1_Strat2-All-All-All.tsv',  
            'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-S205-NS_DE_Strat1-RNA_clusters-0-1_Strat2-All-All-All.tsv',  
            'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-N88-NS_DE_Strat1-RNA_clusters-0-1_Strat2-All-All-All.tsv',  
            'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-M37-NS_DE_Strat1-RNA_clusters-0-1_Strat2-All-All-All.tsv',  
            'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-M45-NS_DE_Strat1-RNA_clusters-0-1_Strat2-All-All-All.tsv',  
            'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-N33-NS_DE_Strat1-RNA_clusters-0-1_Strat2-All-All-All.tsv',  
            'results/Run2023-05-14/DE/Stim/RNA_clusters/Stim-SEB-S42+S205+N88+M37+M45+N33_DE_Strat1-RNA_clusters-0-0_Strat2-All-All-All.tsv'
    )

    # args <- c('/data/vrc_his/douek_lab/projects/RNASeq/2021612_finch/results/Treatment_DESeq2.xls',
  #           '/data/vrc_his/douek_lab/projects/RNASeq/2021612_finch/results/Treatment_DESeq2_sig.xls',
  #           '0.05',
  #           '/data/vrc_his/douek_lab/projects/RNASeq/2021612_finch/results/Treatment/All/DESeq2_results.txt',
  #           '/data/vrc_his/douek_lab/projects/RNASeq/2021612_finch/results/Treatment/Euth_Age-D6/DESeq2_results.txt',
  #           '/data/vrc_his/douek_lab/projects/RNASeq/2021612_finch/results/Treatment/Euth_Age-D18/DESeq2_results.txt',
  #           '/data/vrc_his/douek_lab/projects/RNASeq/2021612_finch/results/Treatment/Euth_Age-D35/DESeq2_results.txt',
  #           '/data/vrc_his/douek_lab/projects/RNASeq/2021612_finch/results/Treatment/Euth_Age-D90/DESeq2_results.txt')
  
  # args <- c('/data/vrc_his/douek_lab/projects/RNASeq/2021612_finch/results/Treatment_GSEA.xls',
  #           '/data/vrc_his/douek_lab/projects/RNASeq/2021612_finch/results/Treatment_GSEA_sig.xls',
  #           '0.05',
  #           '/data/vrc_his/douek_lab/projects/RNASeq/2021612_finch/results/Treatment/All/fgsea_results.txt',
  #           '/data/vrc_his/douek_lab/projects/RNASeq/2021612_finch/results/Treatment/Euth_Age-D6/fgsea_results.txt',
  #           '/data/vrc_his/douek_lab/projects/RNASeq/2021612_finch/results/Treatment/Euth_Age-D18/fgsea_results.txt',
  #           '/data/vrc_his/douek_lab/projects/RNASeq/2021612_finch/results/Treatment/Euth_Age-D35/fgsea_results.txt',
  #           '/data/vrc_his/douek_lab/projects/RNASeq/2021612_finch/results/Treatment/Euth_Age-D90/fgsea_results.txt')

    # args <- c('/data/vrc_his/douek_lab/projects/RNASeq/2022612_Petrovas/results/Strain_GSEA.xls',
    #         '/data/vrc_his/douek_lab/projects/RNASeq/2022612_Petrovas/results/Strain_GSEA_sig.xls',
    #         '0.05',
    #         '/data/vrc_his/douek_lab/projects/RNASeq/2022612_Petrovas/results/Strain/All/fgsea_results.txt',
    #         '/data/vrc_his/douek_lab/projects/RNASeq/2022612_Petrovas/results/Strain/Cell_type-pre_Tfh/fgsea_results.txt',
    #         '/data/vrc_his/douek_lab/projects/RNASeq/2022612_Petrovas/results/Strain/Cell_type-Tfh/fgsea_results.txt')
}else{
  args = commandArgs(trailingOnly=TRUE)
}
xls_file <- args[1]
sig_file <- args[2]
pthresh <- args[3]
full_tab_name <- args[4] ### Boolean
files <- args[5:length(args)]

### Return tab name from FindMarkers output file path
FM_path_2_name <- function(basename){
  dat_type <- strsplit(basename, '_')[[1]][1]
  clusts <- strsplit(basename, '-')[[1]][c(2,3)]
  clust_name <- gsub(paste0(dat_type, '_'), '', basename)
  clust_name <- gsub(paste0('-', clusts[1], '-', clusts[2]), '', clust_name)
  dat_type <- gsub('Data', ' data', dat_type)
  tab_name <- paste0(clust_name, ' ', clusts[1], ' v. ', clusts[2])
  return(tab_name)
}

pseudobulk_path_2_name0 <- function(basename){
  ### Path format from new pseudobulk DE rule
  split_piece <- ifelse(grepl('_GSEA_within_', basename), '_GSEA_within_', '_DE_within_')
  
  test_piece <- strsplit(basename, split_piece)[[1]][1]
  strat_piece <- strsplit(basename, split_piece)[[1]][2]
  test_name <- strsplit(test_piece, '-')[[1]][1]
  value1 <- strsplit(test_piece, '-')[[1]][2]
  value2 <- strsplit(test_piece, '-')[[1]][3]
  strat_name <- strsplit(strat_piece, '-')[[1]][1]
  strat_values <-  strsplit(strat_piece, '-')[[1]][2]
  if(value1 == 'NA' | value2 == 'NA'){
    tab_name <- test_name
  } else{
    tab_name <- paste0(value1, ' v ', value2)
  }
  if(!strat_name %in% c("NA", 'All')){
    tab_name <- paste0(tab_name, ' in ', strat_values)
  }
  return(tab_name)
}

### Return tab name from pseduobulk result file path
pseudobulk_path_2_name <- function(basename, full_tab_name){
  split_piece <- ifelse(grepl('_GSEA_', basename), '_GSEA_', '_DE_')
  
  ### Split file name into its two main portions
  test_piece <- strsplit(basename, split_piece)[[1]][1]
  strat_piece <- strsplit(basename, split_piece)[[1]][2]
  test_name <- strsplit(test_piece, '-')[[1]][1]
  value1 <- strsplit(test_piece, '-')[[1]][2]
  value2 <- strsplit(test_piece, '-')[[1]][3]
  if(value1 == 'NA'){
    value1 <- paste0(test_name, '1')
  }
  if(value2 == 'NA'){
    value2 <- paste0(test_name, '2')
  }
  ### Separate and reduce the strat piece
  strat1 <- gsub('Strat1-', '', strsplit(strat_piece, '_Strat2-')[[1]][1])
  strat2 <- strsplit(strat_piece, '_Strat2-')[[1]][2]
  
  strat1_name <- strsplit(strat1, '-')[[1]][1]
  strat1_value1 <-  strsplit(strat1, '-')[[1]][2]
  strat1_value2 <-  strsplit(strat1, '-')[[1]][3]
  ### If there was no strat1
  if(strat1_name %in% c('NA', 'All')){
    strat1 <- ''
  } else{ ### If strat1 exists and has only 1 value
    if(strat1_value1 == strat1_value2){ 
      strat1_value2 <- ''
      strat1 <- paste0('_in_', strat1_name, '_', strat1_value1, '_', strat1_value2)
    } else{ ### If strat1 exists and has a separate value for each test value, so put them there
      value1 <- paste0(value1, '_in_', strat1_name, '_', strat1_value1)
      value2 <- paste0(value2, '_in_', strat1_name, '_', strat1_value2)
      strat1 <- ''
    }
  }
  ### Repeat for strat 2
  strat2_name <- strsplit(strat2, '-')[[1]][1]
  strat2_value1 <-  strsplit(strat2, '-')[[1]][2]
  strat2_value2 <-  strsplit(strat2, '-')[[1]][3]
  if(strat2_name %in% c('NA', 'All')){
    strat2 <- ''
  } else{ ### One strat for both test values
    if(strat2_value1 == strat2_value2){ 
      strat2_value2 <- ''
      strat2 <- paste0('+in_', strat2_name, '_', strat2_value1, '_', strat2_value2)
    } else{ ### Separate strat for each test value, so put them there
      value1 <- paste0(value1, '_in_', strat2_name, '_', strat2_value1)
      value2 <- paste0(value2, '_in_', strat2_name, '_', strat2_value2)
      strat2 <- ''
    }
  }
  
  ### Unless full_tab_name == F, include the test_name. If its false, that info is hopefully somewhere else, like the file name
  if(full_tab_name){
    tab_name <- paste0(test_name, '_', value1, 'v', value2, '_', strat1, '_', strat2)
  } else{
    tab_name <- paste0(value1, 'v', value2, '_', strat1, '_', strat2)
  }
  tab_name <- paste0(test_name, '_', value1, 'v', value2, '_', strat1, '_', strat2)
  tab_name <- trimws(tab_name, whitespace = '_')
  tab_name <- gsub('__', '_', tab_name)
  if(nchar(tab_name) > 31){
    tab_name <- pseudobulk_path_2_name_shorter(basename, full_tab_name)
  }
  return(tab_name)
}
### Used with pesudobulk_path_2_name if abbrreviation will be done
pseudobulk_path_2_name_shorter <- function(basename, full_tab_name){
  split_piece <- ifelse(grepl('_GSEA_', basename), '_GSEA_', '_DE_')
  
  ### Split file name into its two main portions
  test_piece <- strsplit(basename, split_piece)[[1]][1]
  strat_piece <- strsplit(basename, split_piece)[[1]][2]
  test_name <- strsplit(test_piece, '-')[[1]][1]
  value1 <- strsplit(test_piece, '-')[[1]][2]
  value2 <- strsplit(test_piece, '-')[[1]][3]
  if(value1 == 'NA'){
    value1 <- paste0(test_name, '1')
  }
  if(value2 == 'NA'){
    value2 <- paste0(test_name, '2')
  }
  ### Separate and reduce the strat piece
  strat1 <- gsub('Strat1-', '', strsplit(strat_piece, '_Strat2-')[[1]][1])
  strat2 <- strsplit(strat_piece, '_Strat2-')[[1]][2]
  
  strat1_name <- strsplit(strat1, '-')[[1]][1]
  strat1_value1 <-  strsplit(strat1, '-')[[1]][2]
  strat1_value2 <-  strsplit(strat1, '-')[[1]][3]
  ### If there was no strat1
  if(strat1_name %in% c('NA', 'All')){
    strat1 <- ''
  } else{ ### If strat1 exists and has only 1 value
    if(strat1_value1 == strat1_value2){ 
      strat1_value2 <- ''
      strat1 <- paste0('(', strat1_value1, ')')
    } else{ ### If strat1 exists and has a separate value for each test value, so put them there
      test_name <- paste0(test_name, '(', strat1_name, ')')
      value1 <- paste0(value1, '(', strat1_value1, ')')
      value2 <- paste0(value2, '(', strat1_value2, ')')
      strat1 <- ''
    }
  }
  ### Repeat for strat 2
  strat2_name <- strsplit(strat2, '-')[[1]][1]
  strat2_value1 <-  strsplit(strat2, '-')[[1]][2]
  strat2_value2 <-  strsplit(strat2, '-')[[1]][3]
  if(strat2_name %in% c('NA', 'All')){
    strat2 <- ''
  } else{ ### One strat for both test values
    if(strat2_value1 == strat2_value2){ 
      strat2_value2 <- ''
      strat2 <- paste0('+in_', strat2_name, '_', strat2_value1, '_', strat2_value2)
    } else{ ### Separate strat for each test value, so put them there
      test_name <- paste0(test_name, '(', strat2_name, ')')
      value1 <- paste0(value1, '(', strat2_value1, ')')
      value2 <- paste0(value2, '(', strat2_value2, ')')
      strat2 <- ''
    }
  }
  
  if(full_tab_name){
    tab_name <- paste0(test_name, '-', value1, 'v', value2, strat1, '_', strat2)
    
  }else{
    tab_name <- paste0(value1, 'v', value2, strat1, '_', strat2)
  }
  tab_name <- gsub('RNA_clusters', 'clust', tab_name)
  tab_name <- trimws(tab_name, whitespace = '_')
  tab_name <- gsub('__', '_', tab_name)
  ### If its still too big, do as if full_tab_name == F
  if(nchar(tab_name) > 31){
    tab_name <- paste0(value1, 'v', value2,  strat1, '_', strat2)
    tab_name <- gsub('RNA_clusters', 'clust', tab_name)
    tab_name <- trimws(tab_name, whitespace = '_')
    tab_name <- gsub('__', '_', tab_name)
  }
  return(tab_name)
}

### Get tab names for each file using the above functions
tab_names <- sapply(files, function(file){
  ### File name without extension or path
  basename <- tail(strsplit(basename(file), '/')[[1]], 2)
  basename <- gsub('.RDS', '', basename)
  basename <- gsub('.tsv', '', basename)
  ### Path format from cluster FindMarkers rule
  if(grepl('Cluster_DE/One-one/', file)){
    tab_name <- FM_path_2_name(basename)
  } else if(grepl('_within_', file) ){ ### Older version of file names with only 1 strat possible
     tab_name <- pseudobulk_path_2_name0(basename)
  } else{ ### Newest naming convention with two strats possible
     tab_name <- pseudobulk_path_2_name(basename, full_tab_name)
  }
  return(tab_name)
}
)
if(any(table(tab_names) > 1)){
  warning(paste0("Duplicate test entries: ", toString(names(table(tab_names))[which(table(tab_names) > 1)])
))
}
if(any(nchar(tab_names) > 31)){
  warning('Tab names are too long for excel output: ', toString(tab_names[which(nchar(tab_names) > 31)]))
  print('Attempting to shorten while remaining unique')
}

### Read tsv files into a list of data frames
dat_list <- vector(mode = "list", length = length(files))
for(i in 1:length(files)){
  file <- files[i]
  dat_list[[i]] <- read.table(file, sep = '\t', quote = "", header = T, row.names = 1)
}
names(dat_list) <- tab_names

### Subset each by input p-value threshold
p_val <- c('padj', 'p_val_adj')
p_val <- p_val[which(p_val %in% colnames(dat_list[[1]]))]
dats <- lapply(dat_list, function(dat) dat[which(dat[, p_val] < as.numeric(pthresh)),])

### Create summary data
sig_res <- unique(unlist(lapply(dats, function(dat) row.names(dat))))
out <- as.data.frame(sig_res)
row.names(out) <- sig_res
### Add size, gtf_name columns to out if they are in the dat_list column names
cols <- c('size', 'gene_name_gtf')
cols <- cols[cols %in% colnames(dat_list[[1]])]
for(col in cols){
  out[, col] <- sapply(row.names(out), function(x) dat_list[[1]][x, col])
}
out$n_sig <- sapply(row.names(out), function(gs) get_strats(gs, dats, T))
out$sig_strats <- sapply(row.names(out), function(gs) get_strats(gs, dats))
### Specific to finch data
# out$sig_strats <- gsub('Treatment_Euth_Age-', '', out$sig_strats)
# out$sig_strats <- gsub('Treatment_', '', out$sig_strats)

out <- out[order(out$n_sig, decreasing = T),]
out <- out[, colnames(out)[colnames(out) != 'sig_res']]

dat_list[['Summary']] <- out
dats[['Summary']] <- out
### Write all results
WriteXLS(x = dat_list, ExcelFileName = xls_file, SheetNames = names(dat_list), row.names = T)
### Write significant results
WriteXLS(x = dats, ExcelFileName = sig_file, SheetNames = names(dats), row.names = T)
###

