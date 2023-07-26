library('sys')
library('readr')
library('WriteXLS')
library('data.table')

source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/DE_functions.R')

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
  args <- c('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021612_finch/results/Treatment_DESeq2.xls',
            '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021612_finch/results/Treatment_DESeq2_sig.xls',
            '0.05',
            '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021612_finch/results/Treatment/All/DESeq2_results.txt',
            '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021612_finch/results/Treatment/Euth_Age-D6/DESeq2_results.txt',
            '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021612_finch/results/Treatment/Euth_Age-D18/DESeq2_results.txt',
            '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021612_finch/results/Treatment/Euth_Age-D35/DESeq2_results.txt',
            '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021612_finch/results/Treatment/Euth_Age-D90/DESeq2_results.txt')
  
  args <- c('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021612_finch/results/Treatment_GSEA.xls',
            '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021612_finch/results/Treatment_GSEA_sig.xls',
            '0.05',
            '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021612_finch/results/Treatment/All/fgsea_results.txt',
            '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021612_finch/results/Treatment/Euth_Age-D6/fgsea_results.txt',
            '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021612_finch/results/Treatment/Euth_Age-D18/fgsea_results.txt',
            '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021612_finch/results/Treatment/Euth_Age-D35/fgsea_results.txt',
            '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021612_finch/results/Treatment/Euth_Age-D90/fgsea_results.txt')

    args <- c('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022612_Petrovas/results/Strain_GSEA.xls',
            '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022612_Petrovas/results/Strain_GSEA_sig.xls',
            '0.05',
            '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022612_Petrovas/results/Strain/All/fgsea_results.txt',
            '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022612_Petrovas/results/Strain/Cell_type-pre_Tfh/fgsea_results.txt',
            '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022612_Petrovas/results/Strain/Cell_type-Tfh/fgsea_results.txt')
}else{
  args = commandArgs(trailingOnly=TRUE)
}
xls_file <- args[1]
sig_file <- args[2]
pthresh <- args[3]
files <- args[4:length(args)]

tab_names <- sapply(files, function(file){
  basename <- tail(strsplit(basename(file), '/')[[1]], 2)
  basename <- gsub('.RDS', '', basename)
  basename <- gsub('.tsv', '', basename)
  ### Path format from cluster FindMarkers rule
  if(grepl('Cluster_DE/One-one/', file)){
    dat_type <- strsplit(basename, '_')[[1]][1]
    clusts <- strsplit(basename, '-')[[1]][c(2,3)]
    clust_name <- gsub(paste0(dat_type, '_'), '', basename)
    clust_name <- gsub(paste0('-', clusts[1], '-', clusts[2]), '', clust_name)
    dat_type <- gsub('Data', ' data', dat_type)
    ### Compile human-readable name for excel tabs
    #return(paste0(clust_name, ' ', clusts[1], ' v. ', clusts[2], ' with ', dat_type))
    return(paste0(clust_name, ' ', clusts[1], ' v. ', clusts[2]))
  } else{ 
    ### Path format from original pseudobulk DE rule
    # basename <- tail(strsplit(dirname(file), '/')[[1]], 2)
    # return(paste0(basename[1], '_', basename[2]))
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
}
)
print(tab_names)
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
### Add size, gtf_name
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

