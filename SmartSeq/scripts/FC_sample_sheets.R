
library('sys')
library('dplyr')
library('viridis')
library('readxl')
library('stringr')
library('WriteXLS')
library('data.table')
library('tools')
library('readxl')
library('tools')

source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/sample_sheet_functions.R')
#source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R')

args = commandArgs(trailingOnly=TRUE)
project <- args[1]
investigator <- args[2]
reference <- args[3]
run_path <- args[4]
csv_file <- args[5]

project <- '2021600_kristin'
investigator <- 'Boswell'
reference <- 'hsapiens'
run_path <- '/hpcdata/vrc/vrc1_data/douek_lab/Runs/'
csv_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/All_covariates.csv')

project <- '2021618_galt'
investigator <- 'Mexico'
reference <- 'hsapiens'
run_path <- '/hpcdata/vrc/vrc1_data/douek_lab/Runs/'
csv_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/All_covariates.csv')

# ## bcl2fastq sample sheets from excel docs
# excel_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/', project, '/RD_TLbulkrnaseq_052721.xlsx')
# excel_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021612_finch/RD_TLbulkrnaseq_052721.xlsx'
# csv_file <- paste0(file_path_sans_ext(excel_file), '.csv')

### Read csv
dat <- read.table(csv_file, sep = ',', quote = "", header = T)
### Conforming column names and other formatting
dat <- cell_file_to_csv(dat, project, reference, bcl2fastq2 = F)
### Full flowcell paths using ls and the abbreviated name
if(!'full_flowcell' %in% colnames(dat)){
  full_flowcell <- sapply(unique(dat$FCID), function(x) list.files(path = paste0(run_path), pattern = paste0('*', x), include.dirs = T)[1])
  dat$flowcell_full <- full_flowcell[dat$FCID]
}

##### Now per flowcell bcl2fastq format csv
### Get unique flowcells from dat
a <- unique(dat[, c('FCID', 'flowcell_full')])
flowcells <- a$flowcell_full
names(flowcells) <- a$FCID

### bcl2fastq format header
header <- paste0('[Header]\nDate,', Sys.Date(), '\nWorkflow,bcl2fastq2\nInvestigator,', investigator, '\n[Data]')
#### For each flowcell, write an individual csv file
for(i in 1:length(flowcells)){
  ### Subset rows
  fc_dat <- dat[which(dat$FCID == names(full_flowcell)[i]), ]
  ### Subset columns and change format
  fc_dat <- bcl2fastq2_format(fc_dat)
  file_path <- paste0(run_path, full_flowcell[i], '/', 'SampleSheet.csv')
  write(header, file = file_path, append = F)
  suppressWarnings(write.table(fc_dat, file_path, quote = F, sep = ',', col.names = T, row.names = F, append = T))
}


