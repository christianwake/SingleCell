library('sys')
library('Seurat')
library('ggplot2')
library('dplyr')
library('viridis')
library('data.table')
library('patchwork')
library('PKI')
library('tinytex')
library('dsb')
#library('tidyverse')
#library('pastecs')
library('reticulate')
library('umap')
library('gridExtra')
library('cowplot')
library('logspline')

source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

args = commandArgs(trailingOnly=TRUE)

barcodes <- args[1]
features <- args[2]
matrix <- args[3]
labels_file <- args[4]
RDS_file <- args[5]
hto_string <- args[6]

### The Read10X R function uses input directory path, but for snakemake it is better to have files directly as the input/output.
if(dirname(barcodes) == dirname(features) & dirname(features) == dirname(matrix)){
  data_dir <- dirname(barcodes)
} else{
  warning('Input data files are not in the same directory. This will cause an error downstream.')
}
# data_dir <- '/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/Flu/Data/'
# data_dir <- '/hpcdata/vrc/vrc1_data/douek_lab/Runs/210511_A00243_0102_AHC3YJDSX2/wakecg/CR6/3/outs//filtered_feature_bc_matrix/'
# data_dir <- '/hpcdata/vrc/vrc1_data/douek_lab/Runs/211123_A00243_0139_BHLKCFDSX2/count_unnormalized/outs/count/filtered_feature_bc_matrix/'
# data_dir <- '/hpcdata/vrc/vrc1_data/douek_lab/Runs/211123_A00243_0139_BHLKCFDSX2/CR1/outs/filtered_feature_bc_matrix/'
# hto_string <- 'C0251,C0252,C0253,C0254,C0255,C0256,C0257,C0258,C0259,C0260'
# ## CSV file with features as row.names and a column for each HTO group (plate). Values are 1 if the feature is present for that plate, 0 if not.
# feature_file <- '/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/Flu/features.csv'

print(data_dir)
### Should be 'seuratDemux' (default) or 'myDemux'
demux_method <- 'seuratDemux'

out <- read_dat(data_dir, hto_string)

prot <- out[['prot']]
htos <- out[['htos']]

### Demultiplex
labels <- solve_hashes_main(prot, htos, out_dir, demux_method)
table(labels)
### Save outputs
saveRDS(out, RDS_file)
write.table(labels, labels_file, sep = ',', quote = F, row.names = T, col.names = T)
