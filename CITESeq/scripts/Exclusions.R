library('sys')
library('Seurat')
library('stringr')
library('pheatmap')
library('ggplot2')
library('reticulate'); use_virtualenv("r-reticulate")
library('textshape')
library('dplyr')
library('biomaRt')
library('grid')
library('scales')

source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

args = commandArgs(trailingOnly=TRUE)

cseq_file <- args[1]
key_words <- args[2]
out_file <- args[3]

cseq <- readRDS(cseq_file)

all_names <- row.names(cseq@assays$RNA@counts)
key_words <- strsplit(key_words, ',')[[1]]
exclude_gene_names <- all_names[sapply(all_names, function(x) any(sapply(key_words, function(y) grepl(y, x))))]
exclude_gene_names <- sort(exclude_gene_names)
write.table(exclude_gene_names, out_file, quote = F, sep = ',', row.names = F, col.names = F)
