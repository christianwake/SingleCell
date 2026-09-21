library('sys')
library('Seurat')
library('stringr')
library('pheatmap')
library('ggplot2')
library('umap')
library('textshape')
library('dplyr')
library('biomaRt')
library('grid')
library('scales')
library('rsconnect')

source('/home/cwake/snakemakes/sc_functions.R')
source('/home/cwake/snakemakes/CITESeq_functions.R')

args = commandArgs(trailingOnly=TRUE)
sdat_file <- args[1]
de_file <- args[2]
out_rdata <- args[3]
project <- args[4]
username <- args[5]
server <- args[6]

# sdat_file <- '/home/cwake/projects/2021618_galt/results/Go2/Filtered.RDS'
# de_file <- '/home/cwake/projects/2021618_galt/results/Go2/Cluster_DE.RDS'
# out_rdata <-'/home/cwake/projects/2021618_galt/results/Go2/App/Data.RData'
# project <- '2021618_galt'
# username <- 'wakecg'
# server <- 'rstudio-connect.niaid.nih.gov'

print(sdat_file)
print(de_file)

sdat <- readRDS(sdat_file)
print(names(sdat@reductions))
### removing some assays and reductions
DefaultAssay(sdat) <- 'RNA'
umaps <- names(sdat@reductions)[grepl('umap', names(sdat@reductions))]
print(umaps)
sdat <- DietSeurat(sdat, layers = c('data'), dimreducs = umaps, assays = c('RNA'))

rna_by_rna <- readRDS(de_file)

save.image(out_rdata)

if(project != '' & username != '' & server != ''){
  app_dir <- paste0(dirname(sdat_file), '/App/')
  deployApp(appDir = app_dir, appName = project, account = username, server = server)
}

