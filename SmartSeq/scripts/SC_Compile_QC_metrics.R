### This script assumes you've already completed bcl2fastq, Trimmomatic, STAR and maybe BALDR
library('sys')
#library('Seurat')
library('ggplot2')
library('dplyr')
library('viridis')
library('data.table')
#library('PKI')
library('tinytex')
#library('biomaRt')
#library('harmony')
#library('readxl')
library('stringr')
library('stringi')
library('limma')
library('vegan')

source('/data/vrc_his/douek_lab/snakemakes/sc_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/Utility_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/sample_sheet_functions.R')

sheet <- NA
if(interactive()){
  runs_dir <- '/data/vrc_his/douek_lab/Runs/'
  runs_subdir <- 'demultiplexed_2025'
  project <- '2021600_kristin'
  covs_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/Cell_sheet_2025_sansMissing108.csv')
  # runs_subdir <- 'demultiplexed'
  # project <- '2024615_Boswell'
  # covs_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/Cell_sheet.csv')
  
  mtx_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/data/', 'raw_mtx.RDS')
  #gtf_r_file <- '/data/vrc_his/douek_lab/wakecg/data/gtf.RDS'
  out_csv <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/data/', 'Covariates_QC_metrics.csv')
  #out_txt <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/data/', 'Cell_filters.csv')
  out_pdf <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/data/', 'Histograms.pdf')
  #out_txt <- NA
  #out_pdf <- NA
} else{
  args = commandArgs(trailingOnly=TRUE)
  ### snakemake input
  mtx_file <- args[1]
  ### snakemake params
  runs_dir <- args[2]
  runs_subdir <- args[3]
  project <- args[4]
  covs_file <- args[5]
  #gtf_r_file <- args[2]
  ### snakemake output
  #out_txt <- args[6] ### snakemake pipeline currently (5/14/2023) has these as NA
  out_pdf <- args[6]
  out_csv <- args[7]
}

base_dir  <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/')
### Load gtf info
#gtf <- readRDS(gtf_r_file)

### Read in covariates
sample_info <- read.csv(covs_file, stringsAsFactors = F, header = T,
                 check.names = F, 
                 colClasses = 'character',
                 sep = ',')
### Remove those without FC information
FCID <- c('FC_Name', 'flowcell_full', 'FCID', 'Flowcell_ID')
FCID <- FCID[FCID %in% colnames(sample_info)][1]
sample_info <- sample_info[which(!is.na(sample_info[, FCID])),]
sample_info <- sample_info[which(sample_info[, FCID] != ''),]
### Make sure that cell_id is unique
length(unique(sample_info$Cell_ID)) == length(sample_info$Cell_ID)
row.names(sample_info) <- sample_info$Cell_ID

### Sample_ID from Sample_Name if it doesn't exist
if((!'Sample_ID' %in% colnames(sample_info)) & 'Sample_Name' %in% colnames(sample_info)){
  ### Create abbreviated ID for plotting
  id_match <- paste0('S', 1:length(unique(sample_info$Sample_Name)))
  names(id_match) <- unique(sample_info$Sample_Name)
  sample_info$Sample_ID <- id_match[sample_info$Sample_Name]
}

reps <- split_lanes(sample_info, col = 'Lane', cols = colnames(sample_info))
row.names(reps) <- reps$SeqRep_ID

### Dictionary(ish) to get SeqRep_ID from Cell_ID
key <- lapply(row.names(sample_info), function(x) strsplit(sample_info[x, 'SeqRep_ID'], ' & ')[[1]])
names(key) <- row.names(sample_info)

sample_info$data_dir <- base_dir
### Define more useful columns from parts of existing columns
sample_info$bam_dir <- file.path(base_dir, 'data','bam', sample_info$Cell_ID)
sample_info$count_dir <- file.path(base_dir, 'data','count', sample_info$Cell_ID)

##### Recreate sample_dirs, tab_files, STAR_files, baldr_files, from sample_sheet
#### count csv from bcl2fastq
reps$count_files <- file.path(runs_dir, reps[, FCID], runs_subdir, project, reps$Cell_ID, paste(reps$Cell_Name, reps$S_N, paste0('L00', reps$Lane), 'counts.csv', sep = '_'))
count_files_exist <- file.exists(reps$count_files)
sum(count_files_exist)/dim(reps)[1]
do_count_files_exist <- all(count_files_exist)

### Trimmomatic log files
reps$fastq_dir <- file.path(base_dir, 'data','fastq', 'trimmed', reps[, FCID], reps$SeqRep_ID)
reps$trim_files <- file.path(reps$fastq_dir, 
                   paste0('trim_stderr_', reps$Cell_Name, '_', reps$S_N, '_L00', reps$Lane,'.txt'))
sum(file.exists(reps$trim_files))/dim(reps)[1]
do_trim_files_exist <- all(file.exists(reps$trim_files))
### fastqc summary files
reps$fastqc_summaries_R1 <- file.path(reps$fastq_dir, 
                                          paste0(reps$Cell_Name, '_', reps$S_N, '_L00', reps$Lane, '_R1_paired_fastqc'), 'summary.txt')
reps$fastqc_summaries_R2 <- file.path(reps$fastq_dir, 
                                   paste0(reps$Cell_Name, '_', reps$S_N, '_L00', reps$Lane, '_R2_paired_fastqc'), 'summary.txt')

do_fastqc_files_exist <- all(file.exists(c(reps$fastqc_summaries_R1, reps$fastqc_summaries_R2)))

### Bam files
sample_info$bam_file <- paste0(sample_info$bam_dir, '/', sample_info$Cell_ID, '_Aligned.sortedByCoord.out.bam')
sum(file.exists(sample_info$bam_file))/dim(sample_info)[1]
### tab files
sample_info$tab_files <- file.path(sample_info$count_dir, 'gene_abundances.tab')
sum(file.exists(sample_info$tab_files))/dim(sample_info)[1]
### STAR log files
#sample_info$star_files <- file.path(sample_info$bam_dir,paste0('flowcell_', sample_info$flowcell_full, '_lane_', sample_info$Lane, '_index_', sample_info$UDP_ID, '_star_align_summary.txt'))
sample_info$star_files <- file.path(sample_info$bam_dir, paste0(sample_info$Cell_ID, '_Log.final.out'))
sum(file.exists(sample_info$star_files))/dim(sample_info)[1]

### Make bam_file list per sample
sample_info <- as.data.frame(sample_info)
### Create txt files with lists of bam files, per sample ID
for(sample in unique(sample_info$Sample_ID)){
  bam_file <- paste0(base_dir, '/data/', sample, '.list')
  bams <- sample_info[which(sample_info$Sample_ID == sample), 'bam_file']
  write.table(bams, file = bam_file, quote = F, col.names = F, row.names = F)
}

print(paste0(toString(sample_info[which(!file.exists(sample_info$tab_files)), 'Cell_ID']), ' does not have tab output.'))
### Can't do anything without the count files, so remove cells (from the sample sheet) if their tab file is missing.
sample_info <- sample_info[which(file.exists(sample_info$tab_files)),]

mtx <- readRDS(mtx_file)
print(paste0('Reading mtx file with ', dim(mtx)[1], ' features and ', dim(mtx)[2], ' cells'))
### Diversity
sample_info$Shannon_diversity <- diversity(x = mtx, index = 'shannon', MARGIN = 2)

###### Compile extra QC information
if(do_count_files_exist){
  ### Read bcl2fastq count info (may not exist)
  bfcounts <- as.data.frame(rbindlist(lapply(unique(sample_info$count_files), function(cfile) bcl2fastq_count(cfile, sample_info, FCID = 'Flowcell_ID'))))
  row.names(bfcounts) <- bfcounts$Cell_ID
  ### Match order
  bfcounts <- bfcounts[row.names(sample_info),]
  ### Add bfcounts to sample_info
  sample_info$bcl2fastq_counts <- bfcounts$bcl2fastq_count
}

### Read trimmomatic files
if(do_trim_files_exist){
  trims <- bind_rows(lapply(reps$trim_files, function(x) read_trim_summary(x)))
  trims <- as.data.frame(trims)
  colnames(trims) <- paste0('Trim_', colnames(trims))
  row.names(trims) <- row.names(reps)
  ### Those with empty inputs
  a <- trims[which(trims$Trim_Input_Read_Pairs %in% c('Empty', 'Error', 'Absent')),]
  table(a$Trim_Input_Read_Pairs)
  ### Those actually trimmed
  is <- which(!(trims$Trim_Input_Read_Pairs %in% c('Empty', 'Error', 'Absent')))
  ### Create status column
  trims$Trim_status <- trims$Trim_Input_Read_Pairs
  trims[which(!trims$Trim_status %in%c('Empty', 'Error', 'Absent')), 'Trim_status'] <- 'Pass'
  ### Set rest of columns to NA
  trims[which(trims$Trim_Input_Read_Pairs %in% c('Empty', 'Error', 'Absent')), c("Trim_Input_Read_Pairs","Trim_Both_Surviving","Trim_Forward_Only_Surviving","Trim_Reverse_Only_Surviving","Trim_Dropped")] <- NA
  
  trims$Trim_Input_Read_Pairs <- as.numeric(trims$Trim_Input_Read_Pairs)
  trims$Trim_Both_Surviving <- as.numeric(trims$Trim_Both_Surviving)
  trims$Trim_Forward_Only_Surviving <- as.numeric(trims$Trim_Forward_Only_Surviving)
  trims$Trim_Reverse_Only_Surviving <- as.numeric(trims$Trim_Reverse_Only_Surviving)
  trims$Trim_Dropped <- as.numeric(trims$Trim_Dropped)
  trims[is, 'Trim_loss_freq'] <- (trims[is,]$Trim_Input_Read_Pairs - trims[is, ]$Trim_Both_Surviving)/trims[is, ]$Trim_Input_Read_Pairs
  
  reps <- cbind(reps, trims)  
  ### Remove those with Error or Empty trim logs
  reps <- reps[which(!reps$Trim_status %in% c('Empty', 'Error')),]
  
  trim_cols <- c("Trim_Input_Read_Pairs", "Trim_Both_Surviving", "Trim_Forward_Only_Surviving", "Trim_Reverse_Only_Surviving", "Trim_Dropped")
  for(cid in row.names(sample_info)){
    sids <- key[cid][[1]]
    for(col in trim_cols){
      sample_info[cid, col] <- sum(reps[sids, col]) 
    }
  }
  sample_info[, 'Trim_loss_freq'] <- (sample_info$Trim_Input_Read_Pairs - sample_info$Trim_Both_Surviving)/sample_info$Trim_Input_Read_Pairs
}

### Read fastqc summary files
if(do_fastqc_files_exist){
  ### Fail if either R1 or R2 fails. Warn if either R1 or R2 warns.
  fastqc <- bind_rows(lapply(row.names(reps), function(rn) read_fastqc_summary(reps[rn, 'fastqc_summaries_R1'], reps[rn, 'fastqc_summaries_R2'])))
  fastqc <- as.data.frame(fastqc)
  colnames(fastqc) <- paste0('Fastqc_', colnames(fastqc))
  row.names(fastqc) <- row.names(reps)
  ### Fail if eithre replicate fails. Warn if either replicates warns.
  reps <- cbind(reps, fastqc)  
  for(cid in row.names(sample_info)){
    sids <- key[cid][[1]]
    for(col in colnames(fastqc)){
      sample_info[cid, col] <- ifelse(any(reps[sids, col] == 'FAIL'), 'FAIL', ifelse(any(reps[sids, col] == 'WARN'), 'WARN', 
                                  ifelse(any(reps[sids, col] == 'PASS'), 'PASS', unique(reps[sids,col]))))
    }
  }
}

### Read STAR logs
stars <- star_files_to_df(sample_info$star_files)
stars$Mapped_freq <- stars$N_mapped/stars$N_reads
sample_info <- cbind(sample_info, stars)

if(!is.na(out_pdf) & out_pdf != ''){
  pdf(out_pdf)
  hist(sample_info$Trim_Input_Read_Pairs)
  hist(sample_info$Mapped_freq)
  hist(sample_info$Shannon_diversity)
  dev.off()
}

### Read fastqc summaries


### Cell filtering
# if(!is.na(out_txt) & out_txt != ''){
#   sample_info_full <- sample_info
#   ### To record the number of cells left at each step
#   filteredN <- as.data.frame(data = NA, matrix(nrow = 6, ncol = 3))
#   colnames(filteredN) <- c('Order', 'Filter', 'N_Remaining')
#   filteredN$Order <- 0:(length(row.names(filteredN))-1)
#   filteredN[1, 'N_Remaining'] <- length(row.names(sample_info))
#   ### First, filter cells without too few or too many reads, based on Trim_Input_Read_Pairs
#   sample_info_complete <- sample_info[which(sample_info$Trim_Input_Read_Pairs < (mean(sample_info_full$Trim_Input_Read_Pairs) + 2*sd(sample_info_full$Trim_Input_Read_Pairs))
#                                             & sample_info$Trim_Input_Read_Pairs > (mean(sample_info_full$Trim_Input_Read_Pairs) - 2*sd(sample_info_full$Trim_Input_Read_Pairs))),]
#   filteredN[2, 'N_Remaining'] <- length(row.names(sample_info_complete))
#   ##### All further filters are based on mean and standard deviations of the data just after the first filter
#   ### Filter on Mapped frequency
#   sample_info <- sample_info_complete[which(sample_info_complete$Mapped_freq > (mean(sample_info_complete$Mapped_freq) - 2*sd(sample_info_complete$Mapped_freq))),]
#   filteredN[3, 'N_Remaining'] <- length(row.names(sample_info))
#   ### Filter on trim frequency
#   sample_info <- sample_info[which(sample_info$Trim_loss_freq < (mean(sample_info_complete$Trim_loss_freq) + 2*sd(sample_info_complete$Trim_loss_freq))),]
#   filteredN[4, 'N_Remaining'] <- length(row.names(sample_info))
#   ### Filter on N_Remaining
#   #sample_info <- sample_info[which(sample_info$Trim_Both_Surviving < (mean(sample_info$Trim_Both_Surviving) + 2*sd(sample_info$Trim_Both_Surviving))),]
#   ### Filter on Shannon diversity
#   sample_info <- sample_info[which(sample_info$Shannon_diversity < (mean(sample_info_complete$Shannon_diversity) + 1*sd(sample_info_complete$Shannon_diversity))),]
#   filteredN[5, 'N_Remaining'] <- length(row.names(sample_info))
#   sample_info <- sample_info[which(sample_info$Shannon_diversity > (mean(sample_info_complete$Shannon_diversity) - 1.8*sd(sample_info_complete$Shannon_diversity))),]
#   filteredN[6, 'N_Remaining'] <- length(row.names(sample_info))
#   filteredN$Filter <- c('Original', 'Trim_Input_Read_Pairs', 'Mapped_freq', 'Trim_loss_freq', 'diversity_high', 'diversity_low')
#   
#   write.table(filteredN, out_txt, row.names=F, col.names = T, quote = F, sep = ',')
# }

write.table(sample_info, out_csv, row.names = F, col.names = T, quote  = F, sep = ',')
