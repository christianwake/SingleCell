### This script assumes you've already completed bcl2fastq, Trimmomatic, STAR and maybe BALDR
library('sys')
library('Seurat')
library('ggplot2')
library('dplyr')
library('viridis')
library('data.table')
library('PKI')
#library('tinytex')
#library('biomaRt')
#library('harmony')
library('readxl')
library('stringr')
library('stringi')
#library('limma')
library('vegan')

source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R')

sheet <- NA
if(interactive()){
  runs_dir <- '/hpcdata/vrc/vrc1_data/douek_lab/Runs/'
  project <- '2021600_kristin'
  mtx_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/', 'raw_mtx.rds')
  covs_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/All_covariates.csv')
  gtf_file <- '/hpcdata/vrc/vrc1_data/douek_lab/reference_sets/tenX/Homo_sapiens.GRCh38.93/GRCh38_protein_coding_only/genes/genes.gtf'
  gtf_r_file <- '/hpcdata/vrc/vrc1_data/douek_lab/wakecg/data/GRCh38_reference.RDS'
  out_txt <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/', 'Cell_filters.csv')
  out_pdf <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/', 'Cell_filters.pdf')
  out_csv <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/', 'Covariates_QC_metrics.csv')
} else{
  args = commandArgs(trailingOnly=TRUE)
  ### snakemake input
  mtx_file <- args[1]
  ### snakemake params
  runs_dir <- args[2]
  project <- args[3]
  covs_file <- args[4]
  gtf_file <- args[5]
  ### snakemake output
  gtf_r_file <- args[6]
  out_txt <- args[7]
  out_pdf <- args[8]
  out_csv <- args[9]
}

base_dir  <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/')
if(file.exists(gtf_r_file)){
  gtf <- readRDS(gtf_r_file)
} else {
  gtf <- read_gtf(gtf_file, atts_of_interest = c('gene_id', 'gene_biotype'))
  saveRDS(gtf, gtf_r_file)
}

### Read in covariates
sample_info <- read.csv(covs_file, stringsAsFactors = F, header = T,
                 check.names = F, 
                 colClasses = 'character',
                 sep = ',')
sample_info <- sample_info[which(!is.na(sample_info$FCID)),]
sample_info <- sample_info[which(sample_info$FCID != ''),]

sample_info$data_dir <- base_dir
### Define more useful columns from parts of existing columns
sample_info$cell_dir <- file.path(base_dir, 'data', sample_info$flowcell_full, sample_info$Lane, sample_info$UDP_ID)

### Make bam_file list per sample
sample_info$bam_file <- paste0(sample_info$cell_dir, '/STAR/Aligned.sortedByCoord.out.bam')
any(file.exists(sample_info$bam_file))
all(file.exists(sample_info$bam_file))
sample_info <- as.data.frame(sample_info)
### Create txt files with lists of bam files, per sample ID
for(sample in unique(sample_info$Sample_ID)){
  bam_file <- paste0(base_dir, '/data/', sample, '.list')
  bams <- sample_info[which(sample_info$Sample_ID == sample), 'bam_file']
  write.table(bams, file = bam_file, quote = F, col.names = F, row.names = F)
}

### Make sure that cell_id is unique
length(unique(sample_info$Cell_ID)) == length(sample_info$Cell_ID)

### Recreate sample_dirs, tab_files, STAR_files, baldr_files, from sample_sheet
sample_info$tab_files <- file.path(sample_info$cell_dir, 'STAR', 'gene_abundances.tab')
any(file.exists(sample_info$tab_files))
all(file.exists(sample_info$tab_files))

### Create abbreviated flowcell ID for plotting
fc_match <- paste0('F', 1:length(unique(sample_info$FCID)))
names(fc_match) <- unique(sample_info$FCID)
sample_info$flowcell_plot <- fc_match[sample_info$FCID]

row.names(sample_info) <- sample_info$Cell_ID

print(paste0(toString(sample_info[which(!file.exists(sample_info$tab_files)), 'Cell_ID']), ' does not have tab output.'))
### Can't do anything without the count files, so remove cells (from the sample sheet) if their tab file is missing.
sample_info <- sample_info[which(file.exists(sample_info$tab_files)),]

mtx <- readRDS(out_mtx)
print(paste0('Reading mtx file with ', dim(mtx)[1], ' features and ', dim(mtx)[2], ' cells'))
### Diversity
sample_info$diversity <- diversity(x = mtx, index = 'shannon', MARGIN = 2)

###### Compile extra QC information
#### count csv from bcl2fastq
sample_info$count_files <- file.path(runs_dir, sample_info$flowcell_full, 'demultiplexed', sample_info$Lane, paste(sample_info$flowcell_full, sample_info$Lane, 'counts.csv', sep = '_'))
any(file.exists(unique(sample_info$count_files)))
all(file.exists(unique(sample_info$count_files)))

bfcounts <- as.data.frame(rbindlist(lapply(unique(sample_info$count_files), function(cfile) bcl2fastq_count(cfile, sample_info))))
row.names(bfcounts) <- bfcounts$Cell_ID
### Match order
bfcounts <- bfcounts[row.names(sample_info),]
### Add bfcounts to sample_info
sample_info$bcl2fastq_counts <- bfcounts$bcl2fastq_count

### Get list of STAR log files
sample_info$star_files <- file.path(sample_info$cell_dir, 'STAR', 
                                    paste0('flowcell_', sample_info$flowcell_full, '_lane_', sample_info$Lane, '_index_', sample_info$UDP_ID, '_star_align_summary.txt'))
any(file.exists(sample_info$star_files))
all(file.exists(sample_info$star_files))
sum(file.exists(sample_info$star_files))/dim(sample_info)[1]

### Same for Trimmomatic log files
sample_info$trim_files <- paste0(runs_dir, sample_info$flowcell_full, '/demultiplexed/', sample_info$Lane, '/', sample_info$UDP_ID, '/trim_stderr.txt')
any(file.exists(sample_info$trim_files))
all(file.exists(sample_info$trim_files))
sum(file.exists(sample_info$trim_files))/dim(sample_info)[1]

trims <- bind_rows(lapply(sample_info$trim_files, function(x) read_trim_summary(x)))
colnames(trims) <- paste0('Trim_', colnames(trims))
a <- trims[which(trims$Trim_Input_Read_Pairs %in% c('Empty', 'Error', 'Absent')),]
table(a$Trim_Input_Read_Pairs)

is <- which(!(trims$Trim_Input_Read_Pairs %in% c('Empty', 'Error', 'Absent')))
trims$Trim_status <- trims$Trim_Input_Read_Pairs
trims[which(!trims$Trim_status %in%c('Empty', 'Error', 'Absent')), 'Trim_status'] <- 'Pass'
trims[which(trims$Trim_Input_Read_Pairs %in% c('Empty', 'Error', 'Absent')), c("Trim_Input_Read_Pairs","Trim_Both_Surviving","Trim_Forward_Only_Surviving","Trim_Reverse_Only_Surviving","Trim_Dropped")] <- NA

trims$Trim_Input_Read_Pairs <- as.numeric(trims$Trim_Input_Read_Pairs)
trims$Trim_Both_Surviving <- as.numeric(trims$Trim_Both_Surviving)
trims$Trim_Input_Read_Pairs <- as.numeric(trims$Trim_Input_Read_Pairs)
trims[is, 'Trim_loss_freq'] <- (trims[is,]$Trim_Input_Read_Pairs - trims[is, ]$Trim_Both_Surviving)/trims[is, ]$Trim_Input_Read_Pairs

sample_info <- cbind(sample_info, trims)
### Remove those with Error or Empty trim logs
sample_info <- sample_info[which(!sample_info$Trim_status %in% c('Empty', 'Error')),]

stars <- star_files_to_df(sample_info$star_files)
stars$Mapped_freq <- stars$N_mapped/stars$N_reads
sample_info <- cbind(sample_info, stars)

#match <- sample_info[which(sample_info$bcl2fastq_counts == sample_info$Trim_Input_Read_Pairs)]
#mismatch <- sample_info[which(sample_info$bcl2fastq_counts != sample_info$Trim_Input_Read_Pairs)]

pdf(out_pdf)
hist(sample_info$Trim_Input_Read_Pairs)
hist(sample_info$Mapped_freq)
hist(sample_info$diversity)
dev.off()
sample_info_full <- sample_info
### To record the number of cells left at each step
filteredN <- as.data.frame(data = NA, matrix(nrow = 6, ncol = 3))
colnames(filteredN) <- c('Order', 'Filter', 'N_Remaining')
filteredN$Order <- 0:(length(row.names(filteredN))-1)
filteredN[1, 'N_Remaining'] <- length(row.names(sample_info))
### First filter cells without too few or too many reads, based on Trim_Input_Read_Pairs
sample_info_complete <- sample_info[which(sample_info$Trim_Input_Read_Pairs < (mean(sample_info_full$Trim_Input_Read_Pairs) + 2*sd(sample_info_full$Trim_Input_Read_Pairs))
                                 & sample_info$Trim_Input_Read_Pairs > (mean(sample_info_full$Trim_Input_Read_Pairs) - 2*sd(sample_info_full$Trim_Input_Read_Pairs))),]
filteredN[2, 'N_Remaining'] <- length(row.names(sample_info_complete))
##### All further filters are based on mean and standard deviations of the data just after the first filter
### Filter on Mapped frequency
sample_info <- sample_info_complete[which(sample_info_complete$Mapped_freq > (mean(sample_info_complete$Mapped_freq) - 2*sd(sample_info_complete$Mapped_freq))),]
filteredN[3, 'N_Remaining'] <- length(row.names(sample_info))
### Filter on trim frequency
sample_info <- sample_info[which(sample_info$Trim_loss_freq < (mean(sample_info_complete$Trim_loss_freq) + 2*sd(sample_info_complete$Trim_loss_freq))),]
filteredN[4, 'N_Remaining'] <- length(row.names(sample_info))
### Filter on N_Remaining
#sample_info <- sample_info[which(sample_info$Trim_Both_Surviving < (mean(sample_info$Trim_Both_Surviving) + 2*sd(sample_info$Trim_Both_Surviving))),]
### Filter on Shannon diversity
sample_info <- sample_info[which(sample_info$diversity < (mean(sample_info_complete$diversity) + 1*sd(sample_info_complete$diversity))),]
filteredN[5, 'N_Remaining'] <- length(row.names(sample_info))
sample_info <- sample_info[which(sample_info$diversity > (mean(sample_info_complete$diversity) - 1.8*sd(sample_info_complete$diversity))),]
filteredN[6, 'N_Remaining'] <- length(row.names(sample_info))
filteredN$Filter <- c('Original', 'Trim_Input_Read_Pairs', 'Mapped_freq', 'Trim_loss_freq', 'diversity_high', 'diversity_low')

write.table(filteredN, out_txt, row.names=F, col.names = T, quote = F, sep = ',')
write.table(sample_info, out_csv, row.names = F, col.names = T, quote  = F, sep = ',')

# a<-read.table('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021600_kristin/temp7.txt', sep = ',')
# sample_info$key <- paste(sample_info$flowcell_full, sample_info$Lane, sample_info$UDP_ID, sep = '_')
# a$key <- paste(a$V1, a$V2, a$V3, sep = '_')
# all(a$key %in% sample_info$key)
# 
# b<-sample_info[which(sample_info$key %in% a$key),]
# 
# write.table(b$trim_files, '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021600_kristin/temp.sh', row.names = F, col.names = T, quote  = F, sep = ',')
