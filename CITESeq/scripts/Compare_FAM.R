library('sys')
library('Seurat')
library('stringr')
library('pheatmap')
library('ggplot2')
#library('umap')
#library('textshape')
library('dplyr')
library('biomaRt')
library('grid')
library('scales')
library('rsconnect')
library('data.table')
library('lisi')

source('/data/vrc_his/douek_lab/wakecg/sc_functions.R')
source('/data/vrc_his/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

get_topN_genes_fam <- function(dat, topN, celltype = NA){
  if(!is.na(celltype)){
    subdat <- dat[which(dat$cluster == celltype), ]
  } else{
    subdat <- dat
  }
  subdat <- subdat[order(abs(subdat$avg_log2FC), decreasing = T),]
  tops <- subdat[1:topN, 'gene']
  return(tops)
}

get_topN_genes_fm <- function(dat, topN){
  dat <- dat[order(abs(dat$avg_log2FC), decreasing = T),]
  tops <- row.names(dat[1:topN, ])
  return(tops)
}

res_sub <- 'DS_by_batch'

out_pdf <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/', res_sub, 'Cluster_DE_Comparison.pdf')
rds_files <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/', res_sub, '/Downsample_',
                    c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), '/Cluster_DE/One-all/RNA-predicted.celltype.l1.RDS')
rds_files <- c(rds_files, '/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/Go2/Cluster_DE/One-all/RNA-predicted.celltype.l1.RDS')
rds_files <- rds_files[file.exists(rds_files)]

rna_dats <- lapply(rds_files, function(file) readRDS(file))
names(rna_dats) <- sapply(rds_files, function(file) gsub('/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/', '', gsub('/Cluster_DE/One-all/RNA-predicted.celltype.l1.RDS', '', file)))
### Go2 should've been named Downsample_1.0
names(rna_dats)[which(names(rna_dats) == 'Go2')] <- 'Downsample_1.0'
### Remove "Downsample_"
names(rna_dats) <- sapply(names(rna_dats), function(x) gsub('Downsample_', '', x))

rds_files <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/', res_sub, '/Downsample_',
                    c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), '/Cluster_DE/One-all/prot-predicted.celltype.l1.RDS')
rds_files <- c(rds_files, '/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/Go2/Cluster_DE/One-all/prot-predicted.celltype.l1.RDS')
rds_files <- rds_files[file.exists(rds_files)]
prot_dats <- lapply(rds_files, function(file) readRDS(file))
names(prot_dats) <- sapply(rds_files, function(file) gsub('/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/', '', gsub('/Cluster_DE/One-all/prot-predicted.celltype.l1.RDS', '', file)))
### Go2 should've been named Downsample_1.0
names(prot_dats)[which(names(prot_dats) == 'Go2')] <- 'Downsample_1.0'
### Remove "Downsample_"
names(prot_dats) <- sapply(names(prot_dats), function(x) gsub('Downsample_', '', x))

pdf(out_pdf)

topN <- 50
celltypes <- unique(rna_dats[[1]]$cluster)
celltype <- celltypes[1]
for(celltype in celltypes){
  ### Subset to cell type
  tops <- lapply(rna_dats, function(dat) get_topN_genes_fam(dat, topN, celltype))
  names(tops) <- names(rna_dats)
  pdat <- as.data.frame(sapply(tops, function(top) length(top[which(top %in% tops[['1.0']])])))
  pdat$DS <- row.names(pdat)
  colnames(pdat) <- c('TopN_genes', 'Downsample')
  p <- ggplot(data=pdat, aes(x=Downsample, y=TopN_genes)) +
    geom_bar(stat="identity") + ggtitle(paste0('RNA - ', celltype, ' v. all others'))
  print(p)
}

topN <- 20
for(celltype in celltypes){
  ### Subset to cell type
  tops <- lapply(prot_dats, function(dat) get_topN_genes_fam(dat, topN, celltype))
  names(tops) <- names(prot_dats)
  pdat <- as.data.frame(sapply(tops, function(top) length(top[which(top %in% tops[['1.0']])])))
  pdat$DS <- row.names(pdat)
  colnames(pdat) <- c('TopN_proteins', 'Downsample')
  p <- ggplot(data=pdat, aes(x=Downsample, y=TopN_proteins)) +
    geom_bar(stat="identity") + ggtitle(paste0('protein -', celltype, ' v. all others'))
  print(p)
}

celltypes <- c('B-Mono', 'B-T')
assays <- c('prot', 'RNA')
for(celltype in celltypes){
  for(assay in assays){
    later_path <- paste0('/Cluster_DE/One-one/', assay, '_predicted.celltype.l1-', celltype, '.RDS')
    rds_files <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/', res_sub, '/Downsample_',
                        c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), later_path)
    rds_files <- c(rds_files, paste0('/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/Go2', later_path))
    rds_files <- rds_files[file.exists(rds_files)]
    de_dats <- lapply(rds_files, function(file) readRDS(file))
    names(de_dats) <- sapply(rds_files, function(file) gsub(paste0('/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/', res_sub, '/'), 
                                                              '', gsub(later_path, '', file)))
    ### Go2 should've been named Downsample_1.0
    names(de_dats)[which(names(de_dats) == 'Go2')] <- 'Downsample_1.0'
    ### Remove "Downsample_"
    names(de_dats) <- sapply(names(de_dats), function(x) gsub('Downsample_', '', x))
    
    topN <- 50
    ### Subset to cell type
    tops <- lapply(de_dats, function(dat) get_topN_genes_fm(dat, topN))
    names(tops) <- names(de_dats)
    pdat <- as.data.frame(sapply(tops, function(top) length(top[which(top %in% tops[['1.0']])])))
    pdat$DS <- row.names(pdat)
    colnames(pdat) <- c('TopN_proteins', 'Downsample')
    p <- ggplot(data=pdat, aes(x=Downsample, y=TopN_proteins)) +
      geom_bar(stat="identity") + ggtitle(paste0(assay, '- ', celltype))
    print(p)
  }
}

dev.off()
