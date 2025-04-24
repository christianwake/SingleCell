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
  #demux_method <- 'Trough'
  demux_method <- 'MULTIseqDemux'
  raw_or_norm <- 'raw'
  if(raw_or_norm == 'raw'){
    rds_in <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, 
                     '/results/', qc_name, '/PostQC1.RDS')
    do_norm <- T
  } else{
    rds_in <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, 
                     '/results/', qc_name, '/Combined_batches.RDS')
    do_norm <- F
  }
  thresh_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name,
                        '/Libraseq/', raw_or_norm, '/thresholds_', demux_method, '.csv')
  rds_out <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, 
                   '/results/', qc_name, '/Libraseq/', raw_or_norm, '/calls_', demux_method, '.RDS')
  pdf_out <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, 
                    '/results/', qc_name, '/Libraseq/', raw_or_norm, '/calls_', demux_method, '.pdf')
  probes <- 'S2P_WT_APC,RBD_WT_APC,NTD_WT_APC,gp120_APC,S2P_Beta_PE,RBD_Beta_PE,NTD_Beta_PE,gp120_PE'
 }else{
  args <- commandArgs(trailingOnly=TRUE)
  rds_in <- args[1]
  thresh_file <- args[2]
  rds_out <- args[3]
  tsv_out <- args[4]
  pdf_out <- args[5]
  probes <- args[6]
  do_norm <- as.logical(args[7])
 }

probes <- strsplit(probes, ',')[[1]]
##### Standardize protein names
### space then parentheses -> parentheses
probes <- gsub(' \\(', '\\(', probes)
### spaces, underscores to '.'
probes <- make.names(probes, allow_ = 'F')

downsample_var <- 0.3
if(!interactive()){
  ### Read Seurat data
  sdat <- readRDS(rds_in)
  print('Done reading Seurat object')
  ### Downsample so I can work interactiveley for testing
  assays <- c('RNA', 'prot')
  assays <- assays[which(assays %in% names(sdat@assays))]
  layers <- c('counts', 'data')
  layers <- layers[which(layers %in% names(sdat@assays$RNA@layers))]
  ssub <- DietSeurat(sdat, layers = layers, assays = assays)
  cells <- sample(x = colnames(ssub), size = (length(colnames(ssub)) * downsample_var), replace = F)
  
  ssub <- subset(ssub, cells = cells)
  saveRDS(ssub, file = gsub('.RDS', paste0('_DownSampledTo', downsample_var, '.RDS'), rds_in))
  print('saved downsampled version for testing')
} else{
  sdat <- readRDS(gsub('.RDS', paste0('_DownSampledTo', downsample_var, '.RDS'), rds_in))
}
thresholds <- read.table(thresh_file, sep = ',')
colnames(thresholds) <- c('CR_ID', 'probes', 'thresh')

### Use thresholds to add metadata columns 
# Probe S2P Beta_plus-minus : S2P.Beta.PE
# Probe S2P WA1_plus-minus : S2P
# Probe S2P BA1_plus-minus
# Probe RBD Beta_plus-minus
# Probe RBD WA1_plus-minus
# Probe RBD BA1_plus-minus
# Probe NTD Beta_plus-minus
# Probe NTD WA1_plus-minus
# Probe NTD BA1_plus-minus

label_list <- as.list(rep(NA, length(unique(sdat$CR_ID))))
names(label_list) <- unique(sdat$CR_ID)
for(id in unique(sdat$CR_ID)){
  ### Subset seurat object, protein data and thresholds to the CR_ID
  sd <- subset(sdat, subset = CR_ID == id)
  prot <- sd@assays$prot@data[probes,]
  threshs <- thresholds[which(thresholds$CR_ID == id), ]
  ### Label by threshold
  label_list[[id]] <- dehash_by_threshold_input(id, prot, probes, threshs, hash_colname = 'probes', do_norm = do_norm, only1 = F)
}
print('Combining labels across CRs')
labels <- rbindlist(label_list)
labels <- as.data.frame(labels)
row.names(labels) <- unlist(lapply(label_list, function(x) row.names(x)))
### Reorder rows to match seurat object
labels <- labels[row.names(sdat@meta.data),]

#"S2P.WT.APC"  "RBD.WT.APC"  "NTD.WT.APC"  "gp120.APC"   "S2P.Beta.PE" "RBD.Beta.PE" "NTD.Beta.PE" "gp120.PE"
### Rename columns
colnames(labels) <- sapply(colnames(labels), function(x) gsub('\\.', '_', x))
colnames(labels) <- sapply(colnames(labels), function(x) gsub('WT', 'WA1', x))
colnames(labels) <- sapply(colnames(labels), function(x) gsub('WA1_APC', 'WA1', x))
colnames(labels) <- sapply(colnames(labels), function(x) gsub('RBD_APC', 'RBD', x))
colnames(labels) <- sapply(colnames(labels), function(x) gsub('Beta_APC', 'Beta', x))
colnames(labels) <- sapply(colnames(labels), function(x) gsub('WA1_PE', 'WA1', x))
colnames(labels) <- sapply(colnames(labels), function(x) gsub('RBD_PE', 'RBD', x))
colnames(labels) <- sapply(colnames(labels), function(x) gsub('Beta_PE', 'Beta', x))
colnames(labels) <- sapply(colnames(labels), function(x) paste0(x, '_plusminus'))

print('Adding metadata to seurat object')
sdat <- AddMetaData(sdat, labels)

strains <- c('WA1', 'Beta')
### Create list of strains whose values are vectors of labels column names
slist <- as.list(rep(NA, length(strains)))
names(slist) <- strains
for(strain in strains){
  slist[[strain]] <- colnames(labels)[grepl(strain, colnames(labels))]
}

epitopes <- c('S2P', 'RBD', 'NTD')
### Create list of strains whose values are vectors of labels column names
elist <- as.list(rep(NA, length(epitopes)))
names(elist) <- epitopes
for(epitope in epitopes){
  elist[[epitope]] <- colnames(labels)[grepl(epitope, colnames(labels))]
}
##### Creating new columns that contain epitope and variant binding specificity
### Variants. If any value for that strains is positive, the value for that cell becomes True
print('Adding strain info to seurat object')
for(strain in strains){
  sdat@meta.data[, strain] <- rowSums(sdat@meta.data[, slist[[strain]]] == '+') > 0 
}
### Epitopes. If any column for that epitope is positive, the value for that cell becomes True
print('Adding epitope info to seurat object')
for(epitope in epitopes){
  sdat@meta.data[, epitope] <- rowSums(sdat@meta.data[, elist[[epitope]]] == '+') > 0 
}

#Creating variant and epitope binding columns
# sdat@meta.data$variant_binding <- ifelse(sdat@meta.data$Beta == TRUE,
#                                          ### if Beta is True
#                                          ifelse(sdat@meta.data$WA1 == TRUE | sdat@meta.data$BA1 == TRUE, "cross", "Beta"), 
#                                          ### If Beta is False
#                                          ifelse(sdat@meta.data$WA1 == TRUE,
#                                                 ifelse(sdat@meta.data$BA1 == TRUE, "cross", "WA1"),
#                                                 ifelse(sdat@meta.data$BA1 == TRUE, "BA1", "negative")))

# sdat@meta.data$variant_binding <- ifelse(sdat@meta.data$Beta == TRUE,
#                                          ### if Beta is True
#                                          ifelse(sdat@meta.data$WA1 == TRUE, "cross", "Beta"), 
#                                          ### If Beta is False
#                                          ifelse(sdat@meta.data$WA1 == TRUE, "WA1", 'negative'))

### For each cell (row name), take the first column name that has a True value. 
### This obviously fails for all Falses, or more than one True, but we will override those in the next 2 lines
bool <- sdat@meta.data[, names(slist)]
print('Adding variant binding info to seurat object')
sdat@meta.data$variant_binding <- sapply(row.names(bool), function(rn) colnames(bool)[which(bool[rn,] == T)[1]])
sdat@meta.data[which(rowSums(sdat@meta.data[, names(slist)]) > 1), 'variant_binding'] <- 'cross'
sdat@meta.data[which(rowSums(sdat@meta.data[, names(slist)]) == 0), 'variant_binding'] <- 'negative'

### Negative control filter
neg <- length(which(sdat@meta.data$variant_binding == 'negative')) / length(row.names(sdat@meta.data))

### maybe filter out RBD+ S2P- if there are many
print(table(sdat@meta.data$RBD, sdat@meta.data$S2P, useNA = 'ifany'))

rs_neg <- length(which(sdat@meta.data$RBD == TRUE & sdat@meta.data$S2P == FALSE))/length(row.names(sdat@meta.data))
rs_pos <- length(which(sdat@meta.data$RBD == TRUE & sdat@meta.data$S2P == TRUE))/length(row.names(sdat@meta.data))

### Complex should be minimum
print('Adding epitope binding info to seurat object')
sdat@meta.data$epitope_binding <- ifelse(sdat@meta.data$RBD == TRUE,
                                         ### If RBD is True
                                         ifelse(sdat@meta.data$NTD == TRUE, "complex", "RBD"), 
                                         ### If RBD is False
                                         ifelse(sdat@meta.data$NTD == TRUE,
                                                ### If RBD and NTD are False
                                                "NTD",
                                                ifelse(sdat@meta.data$S2P == TRUE, "S2P", "negative")))

plot_dat <- sdat@meta.data
plot_dat$positive_or_negative <- ifelse(plot_dat$epitope_binding == 'negative', 'Negative', 'Positive')
plot_dat[which(plot_dat$gp120_APC_plusminus == '+' | plot_dat$gp120_PE_plusminus == '+'), 'positive_or_negative'] <-
  'gp120+'

### First plot (just pos/neg)
df <- plot_dat %>% 
  group_by(positive_or_negative) %>% # Variable to be transformed
  count() %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

p1 <- ggplot(df, aes(x = "", y = perc, fill = positive_or_negative)) +
  geom_col() +
  geom_text(aes(label = labels),
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  theme_void() + 
  ggtitle('Libraseq calls', subtitle = paste0('Total N cells: ', dim(plot_dat)[1]))

### 2nd plot (Epitope binding of just postives)
plot_dat <- plot_dat[which(plot_dat$positive_or_negative == 'Positive'),]
### Remove gp120 positives
plot_dat <- plot_dat[which(plot_dat$gp120_APC_plusminus == '-'), ]
plot_dat <- plot_dat[which(plot_dat$gp120_PE_plusminus == '-'), ]

df <- plot_dat %>% 
  group_by(epitope_binding) %>% # Variable to be transformed
  count() %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

p2 <- ggplot(df, aes(x = "", y = perc, fill = epitope_binding)) +
  geom_col() +
  geom_text(aes(label = labels),
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  theme_void() + 
  ggtitle('Libraseq epitope calls (+ only)', subtitle = paste0('Total N cells: ', dim(plot_dat)[1]))

### 3rd plot (S2P + or - within RBD +)
plot_dat <- plot_dat[which(plot_dat$epitope_binding == 'RBD'),]
plot_dat$S2P <- ifelse(plot_dat$S2P == T, '+', '-')
df <- plot_dat %>% 
  group_by(S2P) %>% # Variable to be transformed
  count() %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

p3 <- ggplot(df, aes(x = "", y = perc, fill = S2P)) +
  geom_col() +
  geom_text(aes(label = labels),
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  theme_void() + 
  ggtitle('Libraseq S2P calls (RBD+ only)', subtitle = paste0('Total N cells: ', dim(plot_dat)[1]))

table(sdat@meta.data$epitope_binding, useNA = 'always')
table(sdat@meta.data$variant, useNA = 'always')

### Expect mostly RBD
print(table(sdat@meta.data$epitope_binding))
print(table(sdat@meta.data$variant_binding))
print(table(sdat@meta.data$variant_binding, sdat@meta.data$epitope_binding))

length(which(sdat@meta.data$epitope_binding == 'RBD' & sdat@meta.data$S2P == T))
### These should be NTD not 'complex'
length(which(sdat@meta.data$epitope_binding == 'NTD' & sdat@meta.data$S2P == T))
### So its RBD > NTD > S2P 
saveRDS(sdat, file = rds_out)

dat <- sdat@meta.data[, c('CR_ID', 'bc', 'rna_size', 'ngene', 'prot_size', 'variant_binding', 'epitope_binding', epitopes)]
write.table(dat, tsv_out, quote = F, sep = '\t', col.names = T, row.names = F)

pdf(pdf_out)
print(p1)
print(p2)
print(p3)

dev.off()
