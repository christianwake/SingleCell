library('sys')
library('Seurat')
library('stringr')
library('pheatmap')
library('ggplot2')
#library('textshape')
library('dplyr')
library('biomaRt')
library('grid')
library('scales')
library('rsconnect')

source('/data/vrc_his/douek_lab/wakecg/sc_functions.R')
source('/data/vrc_his/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

args = commandArgs(trailingOnly=TRUE)
out_rdata <- args[1]
rds_files <- args[2:length(args)]
