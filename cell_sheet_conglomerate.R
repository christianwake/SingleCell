
### For when the input sheets have a row per cell rather than per sample or replicate
library('sys')
library('dplyr')
library('viridis')
library('readxl')
library('stringr')
library('stringi')
library('WriteXLS')
library('data.table')
library('tools')
library('readxl')
library('tools')
library('Biostrings')
library('VennDiagram')

if(getwd() == "C:/Users/wakecg/Documents"){
  source('C:/Users/wakecg/Documents/Code/sample_sheet_functions.R')
  source('C:/Users/wakecg/Documents/Code/Utility_functions.R')
} else{
  source('/data/vrc_his/douek_lab/snakemakes/sample_sheet_functions.R')
  source('/data/vrc_his/douek_lab/snakemakes/Utility_functions.R')
}

if(interactive()){
  project <- '2021600_kristin'
  investigator <-'Kristin Boswell'
  reference_rna <- '/data/vrc_his/douek_lab/reference_sets/tenX/Homo_sapiens.GRCh38.93/GRCh38_protein_coding_only_with_AD8/'
  reference_prot <- '/data/vrc_his/douek_lab/Runs/211123_A00243_0139_BHLKCFDSX2/feature_refs_plus_cmos.csv'
  reference_vdj <- '/data/vrc_his/douek_lab/reference_sets/tenX/refdata-cellranger-vdj-GRCh38-alts-ensembl-3.1.0/'
  datatype <- 'RNASeq'
  cellranger <- NA
  bcl2fastq2 <- F
  lib_rep_definition = c('Sample_Name', 'plate', 'well')
  cell_name_definition <- lib_rep_definition
  SeqRep_definition = c('Cell_Name', 'FC_ID', 'Lane')

  reference <- 'hsapiens'
  
  #sample_files <- c('789_31 785_33 768_31 Index.xlsx')
  sample_files <- c('Boswell_SMARTseq_March2020.xlsx', 'Boswell_SMARTseq_Dec2020.xlsx', 'Boswell_SMARTseq_May2021.xlsx', 'Boswell_SMARTseq_June2021.xlsx')

  # project <- '2024615_Boswell'
  # investigator <- 'Kristin Boswell'
  # # reference_rna <- '/data/vrc_his/douek_lab/reference_sets/tenX/Homo_sapiens.GRCh38.93/GRCh38_protein_coding_only_with_AD8/'
  # # reference_prot <- '/data/vrc_his/douek_lab/Runs/211123_A00243_0139_BHLKCFDSX2/feature_refs_plus_cmos.csv'
  # # reference_vdj <- '/data/vrc_his/douek_lab/reference_sets/tenX/refdata-cellranger-vdj-GRCh38-alts-ensembl-3.1.0/'
  # cellranger <- NA
  # bcl2fastq2 <- T
  # lib_rep_definition = c('Sample_Name', 'plate', 'well')
  # cell_name_definition <- lib_rep_definition
  # SeqRep_definition = c('Cell_Name', 'FC_ID', 'Lane')
  # reference <- 'hsapiens'
  # sample_files <- c('Boswell_SMARTseq_July2024.xlsx', 'Boswell_SMARTseq_Nov2024.xlsx')
 }else{
  args = commandArgs(trailingOnly=TRUE)
  
  project <- args[1]
  investigator <- args[2]
  cellranger <- args[3]
  bcl2fastq2 <- args[4]
  lib_rep_definition <- args[5]
  cell_name_definition <- args[6]
  SeqRep_definition <- args[7]
  reference <- args[8]
  sample_files <- args[9:length(args)]
  
  lib_rep_definition <- strsplit(lib_rep_definition, ',')[[1]]
  cell_name_definition <- strsplit(cell_name_definition, ',')[[1]]
  SeqRep_definition <- strsplit(SeqRep_definition, ',')[[1]]
}

datatype <- 'RNASeq'
drop_samples = c()
drop_subjects = c()
sheet <- NA

### bcl2fastq sample sheets from excel docs
if(getwd() == "C:/Users/wakecg/Documents"){
  project_dir <- paste0('C:/Users/wakecg/Documents/SCSeq/', project, '/') ### Laptop
  run_path <- file.path(project_dir, 'Runs')
} else{
  project_dir <- paste0('/data/vrc_his/douek_lab/projects/', datatype, '/', project, '/') ### Skyline
  run_path <- '/data/vrc_his/douek_lab/Runs/'
  #run_path <- '/data/vrc_seq/Illumina-Sequencer-Data/'
}

rep_file <- file.path(project_dir, 'Rep_sheet.csv')
cell_file <- file.path(project_dir, 'Cell_sheet.csv')
sub_file <- file.path(project_dir, 'sub_sheet_', Sys.Date(), '.csv')
base_dir <- file.path(project_dir, 'SampleSheets')

### Read a single file, either csv or xlsx
read_single_file <- function(sample_file, base_dir, bcl2fastq2, sheet = NA){
  print(sample_file)
  ### If a csv file
  if(grepl('\\.csv', sample_file)){
    df <- read.table(file.path(base_dir, sample_file), sep = ',', quote = '', stringsAsFactors = F, header = T)
  } else if(grepl('\\.xls', sample_file)){ ### If an excel file
    #### If sheet is entered read only that sheet. If not, assume each sheet is a separate sample entry, so rbind them all.
    excel_file <- file.path(base_dir, sample_file)
    if(is.na(sheet)){
      df <- as.data.frame(rbindlist(lapply(1:length(excel_sheets(excel_file)), function(i) read_excel(excel_file, sheet = i))) )
    } else {
      df <- read_excel(excel_file, sheet = sheet)
    }
  }
  df <- cell_file_to_csv(df, cellranger, bcl2fastq2, drop_samples, drop_subjects, cell_name_definition)
  return(df)
}

print('Beginning combination of input sheets.')
df <- as.data.frame(rbindlist(lapply(sample_files, function(sample_file) read_single_file(sample_file, base_dir, NA)), fill = T))

### Info may be by replicate or by cell at this point. 
any(duplicated(df$Cell_Name))
### Make it by replicate
### Flowcell &s
if(any(grepl('&', df$FC_Name) | grepl(' and ', df$FC_Name))){
  ### Not tested 
  df$FC_Name <- gsub(' & ', '&', df$FC_Name)
  if(grepl(' and' , df$FC_Name)){
    df$FC_Name <- gsub(' and ', '&', df$FC_Name)
  }
  df <- split_lanes(df, col = 'FC_Name')
}
### For lanes, first convert any dashes to series of &s
### lane values that have a dash
lns <- names(table(df$Lane))[grepl('-', names(table(df$Lane)))]
### List of integer ranges
rgs <- lapply(lns, function(ln) range(as.numeric(strsplit(ln, '-')[[1]])))
### From range to sequence
seqs <- lapply(rgs, function(rg) seq(rg[1], rg[2]))
### From sequence to &-delimited string
key <- lapply(seqs, function(s) gsub(', ', '&', toString(s)))
names(key) <- lns
### Use key to replace Lane values with dashes to their respective & string
is <- which(grepl('-', df$Lane))
for(i in is){
  df[i, 'Lane'] <- key[df[i, 'Lane']]
}
### 
print('Beginning conversion of data frame to per-replicate.')
reps <- split_lanes(df, col = 'Lane')
print('Complete')

### Get full flowcell name from ls path and the abbreviated name in FC_ID (data must be by replicate in case of &s in FC_Name)
if(('FC_ID' %in% colnames(reps)) & (!'FC_Name' %in% colnames(reps))){
  reps <- reps[which(!is.na(reps$FC_ID)),]
  if(run_path == 'C://Users/wakecg/Documents'){
    reps$FC_Name <- reps$FC_ID
  } else{
    FC_Name <- sapply(unique(reps$FC_ID), function(x) list.files(path = paste0(run_path), pattern = paste0('*', x), include.dirs = T)[1])
    reps$FC_Name <- FC_Name[reps$FC_ID]
  }
}

if(project == "2021600_kristin"){
  #df$Cell_Name <- paste(df$FC_Name,  df$Lane,  df$UDP_ID, sep = '_')
  
  ### Specific to Kristins T
  df$Sample_Name <- gsub('-d', '-', df$Sample_Name)
  df$Subject <- sapply(df$Sample_Name, function(x) strsplit(x, '-')[[1]][1])
  df$Time <- gsub('d', '', sapply(df$Sample_Name, function(x) strsplit(x, '-')[[1]][2]))
  # df$Description <- gsub('_ ', '_', df$Description)
  # df$Description <- gsub('SEB activated CD4_CD8', 'SEB activated CD4/CD8', df$Description)
  # df$Description <- gsub('No stim CD4mem_CD8mem', 'No stim CD4mem/CD8mem', df$Description)
  # df$Description <- gsub('No Stim CD4 mem_CD8 mem', 'No stim CD4mem/CD8mem', df$Description)
  
  df$Description <- NULL
  
  ### Kristin's early specificity data
  conf <- read_excel('/data/vrc_his/douek_lab/projects/RNASeq/2021600_kristin/specificity/220921_Confirmed_and_non-specific_TCRs.xlsx', sheet = 'Confirmed TCRs')
  ns <- read_excel('/data/vrc_his/douek_lab/projects/RNASeq/2021600_kristin/specificity/220921_Confirmed_and_non-specific_TCRs.xlsx', sheet = 'Non-specific TCRs')
  tcr1 <- rbind(conf,ns)
  tcr1 <- tcr1[, c('sampleID', 'TCR status', 'cell_prediction', 'specificity')]
  colnames(tcr1) <- c('Cell_Name', 'TCR_status', 'cell_prediction', 'specificity')

  ### Final TCR
  tcr2 <- read_excel('/data/vrc_his/douek_lab/projects/RNASeq/2021600_kristin/SampleSheets/Final_TCR.xlsx',
                     sheet = 'Sheet1')
  colnames(tcr2) <- gsub('sampleID', 'Cell_Name', colnames(tcr2))
  tcr2[, c('geneCounts', 'cloneID', 'TRA', 'TRB', 'UDP')] <- NULL

  tcr <- merge(tcr2, tcr1[, c('Cell_Name', 'TCR_status')], by = 'Cell_Name', all = T)
  colnames(tcr) <- gsub('Well', 'well', colnames(tcr))
  colnames(tcr) <- gsub('Plate', 'plate', colnames(tcr))
  colnames(tcr) <- gsub('Flowcell_Name', 'FC_Name', colnames(tcr))
  colnames(tcr) <- gsub('Sample', 'Sample_Name', colnames(tcr))
  
  ### Index files (combined)
  sams <- read.table('/data/vrc_his/douek_lab/projects/RNASeq/2021600_kristin/specificity/Index_seq.csv', header = T, sep = ',')
  sams$Cell_Name <- paste0(sams$FC_Name, '_', sams$Lane, '_', sams$UDP_ID)
  colnames(sams) <- gsub('Sample_ID', 'Sample_Name', colnames(sams))
  sams$Stim <- gsub('Nuc 27', 'N27', sams$Stim)
  sams$Stim <- gsub('Spike 205', 'S205', sams$Stim)
  sams$Stim <- gsub('Spike 42', 'S42', sams$Stim)
  sams$Stim <- gsub('No Stim', 'NS', sams$Stim)
  sams$Stim <- gsub('Mem 45', 'M45', sams$Stim)
  sams$specificity <- paste0(sams$CD4_or_CD8, '_', sams$Stim)
  sams[, c('Stim', 'CD4_or_CD8')]<- NULL
  
  tcr[, c('well')] <- NULL ### For those confirmed to be equal
  adds <- merge(sams, tcr, by = 'Cell_Name', all = T) # plate and TNF not the same
  adds$Sample_Name <- sapply(1:length(adds$Cell_Name), function(i) ifelse(!is.na(adds[i, 'Sample_Name.x']), adds[i, 'Sample_Name.x'], adds[i, 'Sample_Name.y']))
  adds$TNF <- sapply(1:length(adds$Cell_Name), function(i) ifelse(!is.na(adds[i, 'TNF.x']), adds[i, 'TNF.x'], adds[i, 'TNF.y']))
  adds$specificity <- sapply(1:length(adds$Cell_Name), function(i) ifelse(!is.na(adds[i, 'specificity.x']), adds[i, 'specificity.x'], adds[i, 'specificity.y']))
  adds[, c('Sample_Name.x', 'Sample_Name.y', 'TNF.x', 'TNF.y', 'specificity.x', 'specificity.y', 'Sample_Project')] <- NULL

  colnames(adds) <- gsub('plate.x', 'plate_index', colnames(adds))
  colnames(adds) <- gsub('plate.y', 'plate_tcr', colnames(adds))
  row.names(adds) <- adds$Cell_Name
  #my_venn(c('Seq', 'TCR1', 'TCR2', 'Index'), list(df$Cell_Name, tcr1$Cell_Name, tcr2$Cell_Name, sams$Cell_Name))
  my_venn(c('TCR', 'Sams', 'All'), list(tcr$Cell_Name, sams$Cell_Name, adds$Cell_Name))
  my_venn(c('seq', 'adds'), list(df$Cell_Name, adds$Cell_Name))
  
  ### Updated TCR_status
  excel_file <- '/data/vrc_his/douek_lab/projects/RNASeq/2021600_kristin/specificity/Confirmed_and_High_Confidence_Tcells.xlsx'
  sheets <- excel_sheets(excel_file)
  excel_list <- lapply(sheets, function(sheet) read_excel(excel_file, sheet = sheet, col_names = T))
  names(excel_list) <- sheets
  tstatus <- as.data.frame(rbindlist(excel_list, idcol = 'TCR_status'))
  tstatus[which(tstatus$TCR_status == 'nonspecific'), 'TCR_status'] <- 'non-specific'
  table(tstatus$TCR_status)
  tstatus$TCR_status <- ifelse(tstatus$TCR_status %in% c('CD4_SEB', 'CD8_SEB', 'CD4_NS', 'CD8_NS', 'nonspecific'), tstatus$TCR_status, 'confirmed')
  table(tstatus$TCR_status)
  row.names(tstatus) <- tstatus$Cell_Name
  tstatus$Subject <- sapply(tstatus$Sample_Name, function(x) strsplit(x, '-')[[1]][1])
  tstatus$Time <- gsub('d', '', sapply(tstatus$Sample_Name, function(x) strsplit(x, '-')[[1]][2]))
  
  ### Compare to adds and to df
  my_venn(c('T-status', 'Seq', 'Sams+TCR'), list(tstatus$Cell_Name, df$Cell_Name, adds$Cell_Name))
  ### Time, Subject, cell_prediction, specificity, Lane agree and have no NAs, and plate agrees with plate_tcr not plate_index
  tstatus[, c('Time', 'Subject', 'Sample_Name', 'cell_prediction', 'specificity', 'Lane', 'plate')] <- NULL
  ### well, TNF agree if NAs are removed
  ### TCR_status 
        ### All new (tstatus) TCR_status that are controls are not in the old data (adds)
        ### The two sources agree for all that have data
  
  adds[which(adds$TCR_status =='high confidence'), 'TCR_status'] <- 'confirmed'
  tstatus[, c('geneCounts', 'cloneID', 'TRA', 'TRB', 'UDP_ID')] <- NULL
  adds <- merge(adds, tstatus, by = 'Cell_Name', all = T) # plate and TNF not the same
  adds$TNF <- sapply(1:length(adds$Cell_Name), function(i) ifelse(!is.na(adds[i, 'TNF.x']), adds[i, 'TNF.x'], adds[i, 'TNF.y']))
  adds$TCR_status <- sapply(1:length(adds$Cell_Name), function(i) ifelse(!is.na(adds[i, 'TCR_status.x']), adds[i, 'TCR_status.x'], adds[i, 'TCR_status.y']))
  adds$well <- sapply(1:length(adds$Cell_Name), function(i) ifelse(!is.na(adds[i, 'well.x']), adds[i, 'well.x'], adds[i, 'well.y']))
  
  adds[, c('well.x', 'well.y', 'TNF.x', 'TNF.y', 'TCR_status.x', 'TCR_status.y')] <- NULL
  
  
  colnames(df) <- gsub('plate', 'plate_seq', colnames(df))
  
  dat <- merge(df, adds, by = 'Cell_Name', all = T)
  ### Two which the well's disagree
  #dat <- dat[-which(!is.na(dat$well.x) & !is.na(dat$well.y)  & dat$well.x != dat$well.y),]
 
  ### Do they agree (for those that are not NA for either)?
  # val <- 'plate'
  # all(dat[which(!is.na(dat[, paste0(val, '.x')]) & !is.na(dat[, paste0(val, '.y')])), paste0(val, '.x')] == 
  #       dat[which(!is.na(dat[, paste0(val, '.x')]) & !is.na(dat[, paste0(val, '.y')])), paste0(val, '.y')])
  
  
  dat$FC_Name <- sapply(1:length(dat$Cell_Name), function(i) ifelse(!is.na(dat[i, 'FC_Name.x']), dat[i, 'FC_Name.x'], dat[i, 'FC_Name.y']))
  dat$Sample_Name <- sapply(1:length(dat$Cell_Name), function(i) ifelse(!is.na(dat[i, 'Sample_Name.x']), dat[i, 'Sample_Name.x'], dat[i, 'Sample_Name.y']))
  dat$well <- sapply(1:length(dat$Cell_Name), function(i) ifelse(!is.na(dat[i, 'well.x']), dat[i, 'well.x'], dat[i, 'well.y']))
  dat$Subject <- sapply(1:length(dat$Cell_Name), function(i) ifelse(!is.na(dat[i, 'Subject.x']), dat[i, 'Subject.x'], dat[i, 'Subject.y']))
  dat$Time <- sapply(1:length(dat$Cell_Name), function(i) ifelse(!is.na(dat[i, 'Time.x']), dat[i, 'Time.x'], dat[i, 'Time.y']))
  dat$UDP_ID <- sapply(1:length(dat$Cell_Name), function(i) ifelse(!is.na(dat[i, 'UDP_ID.x']), dat[i, 'UDP_ID.x'], dat[i, 'UDP_ID.y']))
  dat$Notes <- sapply(1:length(dat$Cell_Name), function(i) ifelse(!is.na(dat[i, 'Notes.x']), dat[i, 'Notes.x'], dat[i, 'Notes.y']))
  dat$Lane <- sapply(1:length(dat$Cell_Name), function(i) ifelse(!is.na(dat[i, 'Lane.x']), dat[i, 'Lane.x'], dat[i, 'Lane.y']))
  ### Seq, Index and second TCR sources agree on plate, but first TCR does not always
  dat$plate <- sapply(1:length(dat$Cell_Name), function(i) ifelse(!is.na(dat[i, 'plate_seq']), dat[i, 'plate_seq'], dat[i, 'plate_index']))
  ### 222 rows which are not NA and differ
  sum(dat[which(!is.na(dat[, 'plate']) & !is.na(dat[, 'plate_tcr'])), 'plate'] != 
        dat[which(!is.na(dat[, 'plate']) & !is.na(dat[, 'plate_tcr'])), 'plate_tcr'])
  
  d <- dat[which(!is.na(dat[, 'plate']) & !is.na(dat[, 'plate_tcr'])), c('Cell_Name', 'Subject', 'Time', 'plate', 'plate_tcr')]
  d <- d[which(d$plate != d$plate_tcr),]
  
  dat[, c('Sample_ID', 'Cell_ID', 'FC_Name.x', 'FC_Name.y', 'Sample_Name.x', 'Sample_Name.y', 'well.x', 'well.y', 'Subject.x', 'Subject.y',
          'Time.x', 'Time.y', 'UDP_ID.x', 'UDP_ID.y', 'Notes.x', 'Notes.y', 'Lane.x', 'Lane.y', 'plate_seq', 'plate_index')] <- NULL
  colnames(dat)
  
  ### Columns that are parts of larger columns
  dat$UDP_ID <- sapply(dat$Cell_Name, function(x) strsplit(x, '_')[[1]][length(strsplit(x, '_')[[1]])])
  dat$Lane <- sapply(dat$Cell_Name, function(x) strsplit(x, '_')[[1]][length(strsplit(x, '_')[[1]])- 1])
  dat$FC_Name <- sapply(1:length(dat$Cell_Name), function(i) gsub(paste0('_', dat[i, 'Lane'], '_', dat[i, 'UDP_ID']), '', dat[i, 'Cell_Name']))
  dat$Celltype <- sapply(dat$specificity, function(x) ifelse(is.na(x), NA, strsplit(x, '_')[[1]][1]))
  dat$Stim <- sapply(dat$specificity, function(x) ifelse(is.na(x), NA, strsplit(x, '_')[[1]][2]))
  
  my_venn(c('T-status', 'Seq', 'Sams', 'TCR'), list(tstatus$Cell_Name, df$Cell_Name, 
                                                    sams$Cell_Name, tcr$Cell_Name))
  
  ### Post filtering
  my_venn(c('Sams', 'Seq'), list(sams$Cell_Name, df$Cell_Name))
  qc_name <- 'Run2023-05-14'
  sdat <- readRDS(paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, 
                         '/results/', qc_name, '/PostQC3.RDS'))
  sdat$Cell_Name <- paste(sdat$FC_Name, sdat$Lane, sdat$UDP_ID, sep = '_')
  dat$in_filtered <- ifelse(dat$Cell_Name %in% sdat$Cell_Name, T, F)
  dat$in_sams <- ifelse(dat$Cell_Name %in% sams$Cell_Name, T, F)
  table(dat$in_filtered, dat$in_sams)
  #df <- dat
  #df <- df[which(!is.na(df$Index_Name)),]
}

### Add project info if absent
if(!'Sample_Project' %in% colnames(reps)){
  reps$Sample_Project <- project
}

### FC_ID from, FC_Name mapping (Needs to be done to replicate data, not cell data)
### Remove those with NA FC_Name
reps <- reps[which(!is.na(reps$FC_Name)),]
names_map <- paste0('FC', 1:length(unique(reps[, 'FC_Name'])))
names(names_map) <- unique(reps[, 'FC_Name'])
reps[, 'FC_ID'] <- names_map[reps[, 'FC_Name']]

### Sample_ID from Sample_Name if it doesn't exist
if((!'Sample_ID' %in% colnames(reps)) & 'Sample_Name' %in% colnames(reps)){
  ### Create abbreviated ID for plotting
  id_match <- paste0('S', 1:length(unique(reps$Sample_Name)))
  names(id_match) <- unique(reps$Sample_Name)
  reps$Sample_ID <- id_match[reps$Sample_Name]
}

### Cell_ID from Cell_Name if it doesn't exist
if((!'Cell_ID' %in% colnames(reps)) & 'Cell_Name' %in% colnames(reps)){
  ### Create abbreviated ID for plotting
  id_match <- paste0('C', 1:length(unique(reps$Cell_Name)))
  names(id_match) <- unique(reps$Cell_Name)
  reps$Cell_ID <- id_match[reps$Cell_Name]
}

### SeqRep_Name and ID
SeqRep_definition <- SeqRep_definition[which(SeqRep_definition %in% colnames(reps))]
reps[, 'SeqRep_Name'] <- sapply(row.names(reps), function(rn) gsub(', ', '_', toString(reps[rn, SeqRep_definition])))
print(table(table(reps$SeqRep_Name)))
### Create abbreviated ID for plotting
id_match <- paste0('SR', 1:length(unique(reps$SeqRep_Name)))
names(id_match) <- unique(reps$SeqRep_Name)
reps$SeqRep_ID <- id_match[reps$SeqRep_Name]

### Write Sample Number, only used by bcl2fastq2 and quite annoying.
reps$S_N <- NA
flowcells <- unique(reps$FC_Name)
for(i in 1:length(flowcells)){
  flowcell <- flowcells[i]
  subdat <- reps[which(reps$FC_Name == flowcell), ]
  reps[row.names(subdat), 'S_N'] <- paste0('S', 1:length(row.names(subdat)))
}

hto <- c("Hashtag", "HTO")
hto <- hto[which(hto %in% colnames(reps))]
cols <- c('Sample_Project', 'Operator', 'Lane', 'Lane_string', 'Sample_Name', colnames(reps)[grepl('_ID', colnames(reps))], hto, 'FC_Name', 'FC_ID', 
          'FC_Name_prot', 'FC_ID_prot', 'FC_Name_VDJ', 'FC_ID_VDJ')
cols <- unique(cols[which(cols %in% colnames(reps))])
reps <- reps[, c(cols, colnames(reps)[which(!colnames(reps) %in% cols)])]

### Choose the name of the gene expression index column from a few options
if(!exists('GEx_index_name')){
  GEx_index_cols <- c('GEx_index_seq', 'index', 'index2')
  GEx_index_name <- GEx_index_cols[which(GEx_index_cols %in% colnames(reps))]
}

### For checking if pipeline files already exist
# ssss_dirs <- file.path(run_path, reps$FC_Name, 'demultiplexed', reps$Lane, reps$UDP_ID)
# b2q_dirs <- file.path(run_path, reps$FC_Name, 'bcl2fastq2', project, reps$Cell_ID)
# proj_dirs <- file.path(project_dir, 'data', 'fastq', 'trimmed', reps$Cell_ID)
# 
# reps$ssss_fq <- paste0(ssss_dirs, '/read1_index_', reps$UDP_ID, '.fq.gz')
# reps$b2f_fq <- paste0(b2q_dirs, '/', reps$Cell_Name, '_', reps$S_N, '_L00', reps$Lane, '_R1_001.fastq.gz')
# 
# sum(file.exists(reps$ssss_fq))
# sum(file.exists(reps$b2f_fq))
# sum(file.exists(proj_dirs))
# 
# hist(file.size(reps$ssss_fq))
# hist(file.size(reps$b2f_fq))
# 
# reps$trim_log <- paste0(proj_dirs, '/trim_stderr.txt')
# reps$trim_log_exists <- file.exists(reps$trim_log)
# 
# ci <- reps[which(reps$trim_log_exists == T), ]
# logs <- reps[which(reps$trim_log_exists == T), 'trim_log']
# 
# #err <- sapply(row.names(ci), function(rn) grep('Invalid FASTQ',readLines(ci[rn, 'trim_log']))) != 'integer(0)'
# err <- sapply(row.names(ci), function(rn) grep('Completed successfully', readLines(ci[rn, 'trim_log']))) == 'integer(0)'
# reps[, 'trim_error'] <- NA
# reps[names(err), 'trim_error'] <- err
# table(reps$trim_error, useNA = 'always')


### Create sample-level info from cell-level info
# sample_info <- reps
# cols <- c('Operator', 'Sample_ID', 'Sample_Name', 'Subject', 'Time', 'Description', 'plate', 'Lane', 'FC_ID', 'FC_Name', 
#           'Library_Preparation_Method', 'Recipe', 'Date_Sequenced')
# sample_info <- sample_info[, cols[which(cols %in% colnames(sample_info))]]
# sample_info <- sample_info[!duplicated(sample_info),]

print(paste0('Writing replicate sheet to: ', rep_file))
write.table(reps, rep_file, row.names = F, col.names = T, quote  = F, sep = ',')

flowcells <- unique(reps$FC_Name)
### Loop over flowcells to write per flowcell SampleSheet used by cellranger mkfastq
submit_file <-  file.path(project_dir, 'bcl2fastq_submit.sh')
file.remove(submit_file)

if(!is.na(bcl2fastq2)){
  ### First print if they exist already
  file.exists(sapply(flowcells, function(flowcell) file.path(run_path, flowcell, 'CellSheet.csv')))
  for(i in 1:length(flowcells)){
    flowcell <- flowcells[i]
    ### Declare flowcell file path
    fc_csv <- file.path(run_path, flowcell, 'CellSheet.csv')
    ### Subset to rows with this flowcell
    fcsub <- reps[which(reps[, 'FC_Name'] == flowcell), ]
  
    if(all(fcsub$Lane == '*')){
      Lane <- NA
    } else{
      Lane <- 'Lane'
    }
    
    ### Get bcl2fastq2 format
    # samplesheet <- bcl2fastq2_format(fcsub, Lane = Lane, sample_name = 'Cell_Name',
    #                                  sample_id = 'Cell_ID', sample_project = 'Sample_Project',
    #                                  indeces = GEx_index_name, RC = c(F, T))
    samplesheet <- bcl2fastq2_format(fcsub, Lane = Lane, sample_name = 'Cell_Name',
                                     sample_id = 'SeqRep_ID', sample_project = 'Sample_Project',
                                     indeces = GEx_index_name, RC = c(F, T))
    
    header <- paste0('[Header]\nDate,', Sys.Date(), '\nWorkflow,bcl2fastq2\nInvestigator,', investigator, '\n[Data]')
    cat(header, sep = '\n', file = fc_csv)
    suppressWarnings(write.table(samplesheet, fc_csv, row.names = F, col.names = T, quote  = F, sep = ',', append = T))
    ### Add submission line to sh file
    cat(paste0('sbatch bcl2fastq.sbatch ', flowcell), sep = '\n', file = submit_file, append = T)
  }
}

### Switch to cell-level information
print('Beginning conversion of sheet from per-replicate to per-cell.')
cells <- replicates2samples(reps, name_col = 'Cell_Name', evals = colnames(reps))

print(paste0('Writing cell sheet to: ', cell_file))
cells$`Sample#` <- NULL
write.table(cells, cell_file, row.names = F, col.names = T, quote  = F, sep = ',')
df <- read.table(cell_file, header = T, sep = ',')

### scratch
#fcsub <- df[which(df[, 'FC_Name'] == flowcells[1]), ]
#path1 <- file.path(run_path, flowcells[1], 'bcl2fastq2', project)
#files1 <- list.files(path1, 'R1', recursive = T, full.names = T)
#files1 <- files1[which(file.size(files1) > 0)]
#length(files1)

#path2 <- file.path(run_path, flowcells[1], 'demultiplexed', '2')
#files2 <- list.files(path2, 'read1', recursive = T, full.names = T)
#files2 <- files2[which(file.size(files2) > 10000)]
#length(files2)

#hist(file.size(files1))
#hist(file.size(files2))

### Write a shell script with the following: sh NextSeq-HiSeq3K4K_NexteraUDPSubmitter.sh Flowcell_Nam 'Lane'
#submitter_script <- '/data/vrc_his/douek_lab/wakecg/preprocessing/NextSeq-HiSeq3K4K_NexteraUDPSubmitter.sh'
#shell_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/Submitter.sh')
#subdat <- df[, c('FC_Name', 'Lane')]
#subdat <- subdat[!duplicated(subdat), ]

#subdat$completion_status <- sapply(row.names(subdat), function(rn) 
#  all(df[which(df$Lane == subdat[rn, 'Lane'] & df$FC_Name == subdat[rn, 'FC_Name']), 'ssss_fq'])
#)
#subdat <- subdat[which(subdat$completion_status != T),]

#subdat$submitter <- paste0("sh ", submitter_script, " ", subdat$FC_Name, " \'", subdat$Lane, "\'")

#writeLines(subdat$submitter,shell_file)
