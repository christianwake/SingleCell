library('sys')
library('dplyr')
library('data.table')

source('/data/vrc_his/douek_lab/snakemakes/sc_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/Utility_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/sample_sheet_functions.R')

if(interactive()){
  runs_dir <- '/data/vrc_his/douek_lab/Runs/'

  runs_subdir <- 'demultiplexed_2025'
  project <- '2021600_kristin'
  in_csv <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/Cell_sheet_2025.csv')

  # runs_subdir <- 'demultiplexed'
  # project <- '2024615_Boswell'
  # in_csv <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/Cell_sheet.csv')
  
  mtx_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/data/', 'raw_mtx.rds')
  out_csv <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/data/', 'Covariates_QC_metrics.csv')
  out_pdf <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/data/', 'Histograms.pdf')
} else{
  args = commandArgs(trailingOnly=TRUE)
  in_csv <- args[1]
  out_csv <- args[2]
}

### Read in covariates
cells <- read.csv(in_csv, stringsAsFactors = F, header = T,
                        check.names = F, 
                        colClasses = 'character',
                        sep = ',')

### Make sure that cell_id is unique
length(unique(cells$Cell_ID)) == length(cells$Cell_ID)
row.names(cells) <- cells$Cell_ID

evals <- c('Sample_Project', 'Operator', 'Sample_ID', 'FC_ID', 'FC_Name', 'Description', 'Recipe', 'Date_Sequenced', 'Cell_ID')
samples <- replicates2samples(cells, name_col = 'Sample_Name', evals = evals)
row.names(samples) <- samples$Sample_ID
### Confirm that Sample_Name is now unique
all(table(samples$Sample_Name) == 1)

ns <- table(cells$Sample_ID)
ns <- ns[row.names(samples)]
samples$Ncells <- ns
cols <- c('Sample_ID', 'Sample_Name', 'Sample_Project', 'Operator', 'Description', 'Recipe', 'FC_ID', 'FC_Name', 'Date_Sequenced', 'Ncells')
cols <- cols[which(cols %in% colnames(samples))]
samples <- samples[, ]

write.table(samples, out_csv, row.names = F, col.names = T, quote  = F, sep = ',')
