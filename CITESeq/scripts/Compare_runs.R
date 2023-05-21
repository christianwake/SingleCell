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
library('data.table')
library('lisi')
library('WriteXLS')

source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

# args = commandArgs(trailingOnly=TRUE)
# out_pdf <- args[1]
# rds_files <- args[2:length(args)]
res_sub <- 'DS_by_batch'
#res_sub <- ''
cor_pdf <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/', res_sub, '/DSB_correlations.pdf')
cor_hist__pdf <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/', res_sub, '/DSB_correlations_hist.pdf')
cor_xls <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/', res_sub, '/DSB_correlations.xls')
cor_dir <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/', res_sub, '/DSB_correlations/')
out_pdf <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/', res_sub, '/Comparison.pdf')
rds_files <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/', res_sub, '/Downsample_', 
       c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), '/Mapped.RDS')
rds_files <- c(rds_files, '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/Go2/Mapped.RDS')
rds_files <- rds_files[file.exists(rds_files)]
#rds_files <- rds_files[1:2]
names(rds_files) <- sapply(rds_files, function(file) gsub(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/', res_sub, '/'), '', gsub('/Mapped.RDS', '', file)))
names(rds_files) <- sapply(names(rds_files), function(file) gsub('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/', '', gsub('/Mapped.RDS', '', file)))

### Go2 should've been named Downsample_1.0
names(rds_files)[which(names(rds_files) == 'Go2')] <- 'Downsample_1.0'
### Remove "Downsample_"
names(rds_files) <- sapply(names(rds_files), function(x) gsub('Downsample_', '', x))

names_key_file <- '/hpcdata/vrc/vrc1_data/douek_lab/wakecg/Biolegend_panel.csv'
names_key_dat <- read.table(names_key_file, sep = ',', quote = "", header = T)
#names_key_dat$Name_brev <- gsub('', names_key_dat$Name
names_key <- names_key_dat$Name
names(names_key) <- names_key_dat$ID

get_topN_genes <- function(sdat, topN, datatype = 'RNA', celltype = NA){
  if(!is.na(celltype)){
    subdat <- subset(sdat, cells = names(which(sdat$predicted.celltype.l1 == celltype)))
  } else{
    subdat <- sdat
  }
  sums <- rowSums(subdat@assays[[datatype]])
  sums <- sums[order(sums, decreasing = T)]
  tops <- names(sums[1:topN])
  return(tops)
}

read_diet_seurat <- function(rds_file, umaps, assays = c('RNA', 'prot')){
  sdat <- readRDS(rds_file)
  DefaultAssay(sdat) <- assays[1]
  sdat <- DietSeurat(sdat, counts = F, dimreducs = umaps, assays = assays)
  return(sdat)
}

evaluating_downsample_fraction <- function(ncounts, ds){
  ### Get the cells from this downsample, but for the last downsample set (or just '1.0')
  orig <- ncounts[['1.0']]
  ### Remove cells in DSX but not in DS1.0?
  n <- length(which(!(names(ncounts[[ds]]) %in% names(orig))))
  if(n > 0){
    print(paste0(n, ' cells in the downsampled set but not in the original....'))
  }
  ncounts[[ds]] <- ncounts[[ds]][which(names(ncounts[[ds]]) %in% names(orig))]
  orig <- orig[names(ncounts[[ds]])]
  ### fractions
  frac <- as.data.frame(ncounts[[ds]] / orig)
  frac$DS <- ds
  colnames(frac) <- c('actual_fraction', 'downsample_category')
  nna <- sum(is.na(frac$actual_fraction))
  if(nna > 0){
    print(paste0(nna, ' NAs for some reason')) 
  }
  return(frac)
}

ilisi_dat <- function(dat, sdat_name){
  dat$downsample <- sdat_name
  colnames(dat) <- c('ilisi', 'downsample')
  return(dat)
}

umaps <- c("RNA_umap", "prot_umap", "dsb_wnn_umap")
### Correlation plots for each protein
# pdf(cor_pdf)
# sdats <- lapply(rds_files[c('0.5', '1.0')],function(file) read_diet_seurat(file, umaps, c('prot')))
# b_summary <- as.data.frame(matrix(nrow = length(row.names(sdats[['1.0']]@assays$prot)), ncol = 3))
# row.names(b_summary) <- row.names(sdats[['1.0']]@assays$prot)
# colnames(b_summary) <- c('name', 'cor', 'rsq')
# all_summary <- b_summary
# for(i in 1:length(row.names(sdats[['1.0']]@assays$prot))){
#   adt <- row.names(sdats[[1]]@assays$prot@data)[i]
#   x <- sdats[['1.0']]@assays$prot@data[adt, ]
#   y <- sdats[['0.5']]@assays$prot@data[adt, ]
#   cells <- intersect(names(x), names(y))
#   print(paste0(length(cells), ' cells in intersection of 1.0 (', length(names(x)), ') and 0.5 (', length(names(y)), ').'))
#   x <- x[cells]
#   y <- y[cells]
#   correlation <- signif(cor(x,y), digits = 3)
#   rsq <- signif(summary(lm(x ~y))$adj.r.squared, digits = 3)
#   ### Update the summary data frame
#   all_summary[adt, 'name'] <- names_key[adt]
#   all_summary[adt, 'cor'] <- correlation
#   all_summary[adt, 'rsq'] <- rsq
#   plotdat <- as.data.frame(matrix(nrow = length(x), ncol = 3))
#   plotdat[, 1] <- x
#   plotdat[, 2] <- y
#   plotdat[, 3] <- sdats[['1.0']]@meta.data[names(x), 'predicted.celltype.l1']
#   colnames(plotdat) <- c('DS1.0', 'DS0.5', 'celltype')
#   p1 <- ggplot(plotdat, aes(DS1.0, DS0.5, colour = celltype)) +
#     geom_point() + geom_abline(slope=1, intercept=0) + ggtitle(names_key[adt], subtitle = paste0('corr: ', correlation, ', rsquared: ', rsq))
#   print(p1)
#   
#   cells <- colnames(sdats[['1.0']])[which(sdats[['1.0']]$predicted.celltype.l1 == 'B')]
#   cells <- cells[cells %in% names(x)]
#   x <- x[cells]
#   y <- y[cells]
#   correlation <- signif(cor(x,y), digits = 3)
#   rsq <- signif(summary(lm(x ~y))$adj.r.squared, digits = 3)
#   ### Update the summary data frame
#   b_summary[adt, 'name'] <- names_key[adt]
#   b_summary[adt, 'cor'] <- correlation
#   b_summary[adt, 'rsq'] <- rsq
#   
#   plotdat <- as.data.frame(matrix(nrow = length(x), ncol = 3))
#   plotdat[, 1] <- x
#   plotdat[, 2] <- y
#   plotdat[, 3] <- sdats[['1.0']]@meta.data[names(x), 'predicted.celltype.l1']
#   colnames(plotdat) <- c('DS1.0', 'DS0.5', 'celltype')
#   p2 <- ggplot(plotdat, aes(DS1.0, DS0.5, colour = celltype)) +
#     geom_point() + geom_abline(slope=1, intercept=0) + ggtitle(names_key[adt], subtitle = paste0('corr: ', correlation, ', rsquared: ', rsq))
#   print(p2)
#   
#   name_short <- names_key_dat[which(names_key_dat$ID == adt), 'Name_short']
#   png(paste0(cor_dir, '/', name_short, '.png'))
#   print(p2)
#   dev.off()
# }
# dev.off()
# 
# WriteXLS(c('b_summary', 'all_summary'), cor_xls, c('B_cells', 'All_celltypes'), row.names = T, col.names = T)
# pdf(cor_hist__pdf)
# ggplot(b_summary, aes(x = cor)) + 
#   geom_histogram() + ggtitle('B cell correlation values')
# ggplot(b_summary, aes(x = rsq)) + 
#   geom_histogram() + ggtitle('B cell r-squared values')
# ggplot(all_summary, aes(x = cor)) + 
#   geom_histogram() + ggtitle('All cell types, correlation values')
# ggplot(all_summary, aes(x = rsq)) + 
#   geom_histogram() + ggtitle('All cell types, r-squared values')
# dev.off()

# out_rdata <-  '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/DS_by_cell/Go2/App/Data.RData'
sdats <- lapply(rds_files, function(file) read_diet_seurat(file, umaps))
names(sdats) <- names(rds_files)
for(i in 1:length(sdats)){
  print(i)
  sdat <- sdats[[i]]
  print(mean(colSums(sdat@assays$prot@data)))
  print(mean(sdat$nCount_prot))
  print(mean(sdat$nCount_dsb))
}

#subdats <- lapply(sdats, function(sdat) subset(sdat, downsample = 100))
subdats <- sdats
subdats[[1]] <- subset(subdats[[1]], downsample = 100)
subdats <- lapply(subdats, function(sdat) subset(sdat, cells = colnames(subdats[[1]]@assays$RNA)))
saveRDS(subdats, paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/', res_sub, '/Comparison.RDS'))

### Check that per cell fraction of total is the same per downsample
ncounts <- lapply(sdats, function(sdat) sdat$nCount_RNA)
#saveRDS(ncounts, '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/DS_by_cell/downsample_fractions.RDS')
pdf(out_pdf)
dat <- rbindlist(lapply(names(ncounts), function(ds) evaluating_downsample_fraction(ncounts, ds)))
p <- ggplot(dat, aes(x=downsample_category, y=actual_fraction)) +
 geom_boxplot() + ggtitle('Downsample fraction reduced per cell (RNA)')
p

### Check that per cell fraction of total is the same per downsample
ncounts <- lapply(sdats, function(sdat) sdat$nCount_prot)
dat <- rbindlist(lapply(names(ncounts), function(ds) evaluating_downsample_fraction(ncounts, ds)))
p <- ggplot(dat, aes(x=downsample_category, y=actual_fraction)) +
  geom_boxplot() + ggtitle('Downsample fraction reduced per cell (protein)')
p

### Bind all the predicted cell type information together for plotting
freqdat <- rbindlist(lapply(sdats, function(sdat) as.data.frame(table(sdat$predicted.celltype.l1))), idcol = 'DS')
colnames(freqdat) <- c("DS", "Cell_Type", "N_cells")

### Stacked barchart
ggplot(freqdat, aes(fill=Cell_Type, y=N_cells, x=DS)) +
  geom_bar(position="fill", stat="identity")

topN <- 50
celltypes <- unique(sdats[[1]]$predicted.celltype.l1)
celltype <- celltypes[1]
for(celltype in celltypes){
  ### Subset to cell type
  tops <- lapply(sdats, function(sdat) get_topN_genes(sdat, topN, 'RNA', celltype))
  names(tops) <- names(sdats)
  dat <- as.data.frame(sapply(tops, function(top) length(top[which(top %in% tops[['1.0']])])))
  dat$DS <- row.names(dat)
  colnames(dat) <- c('TopN_genes', 'Downsample')
  p <- ggplot(data=dat, aes(x=Downsample, y=TopN_genes)) +
    geom_bar(stat="identity") + ggtitle(celltype)
  print(p)
}

topN <- 20
for(celltype in celltypes){
  ### Subset to cell type
  tops <- lapply(sdats, function(sdat) get_topN_genes(sdat, topN, 'prot', celltype))
  names(tops) <- names(sdats)
  dat <- as.data.frame(sapply(tops, function(top) length(top[which(top %in% tops[['1.0']])])))
  dat$DS <- row.names(dat)
  colnames(dat) <- c('TopN_proteins', 'Downsample')
  p <- ggplot(data=dat, aes(x=Downsample, y=TopN_proteins)) +
    geom_bar(stat="identity") + ggtitle(celltype)
  print(p)
}
ilisi_rna <- lapply(sdats, function(sdat) compute_lisi(sdat@reductions[['RNA_umap']]@cell.embeddings, sdat@meta.data[, 'RNA_clusters', drop = F], c('RNA_clusters')))
ilisi_prot <- lapply(sdats, function(sdat) compute_lisi(sdat@reductions[['prot_umap']]@cell.embeddings, sdat@meta.data[, 'prot_clusters', drop = F], c('prot_clusters')))
names(ilisi_rna) <- names(sdats)
names(ilisi_prot) <- names(sdats)
ilisi_rna <- as.data.frame(rbindlist(lapply(names(sdats), function(sdat_name) ilisi_dat(ilisi_rna[[sdat_name]], sdat_name))))
ilisi_prot <- as.data.frame(rbindlist(lapply(names(sdats), function(sdat_name) ilisi_dat(ilisi_prot[[sdat_name]], sdat_name))))

print('plot ilisi')
p <- ggplot(ilisi_rna, aes(x=downsample, y=ilisi)) +
  geom_boxplot() + ggtitle('RNA cluster ilisi')
p

p <- ggplot(ilisi_prot, aes(x=downsample, y=ilisi)) +
  geom_boxplot() + ggtitle('Protein cluster ilisi')
p
dev.off()

# Idents(sdat) <- sdat[['predicted.celltype.l1']]
# DefaultAssay(sdat) <- 'RNA'
# mono <- FindMarkers(sdat, assay = 'RNA', slot = 'data', logfc.threshold = 0, ident.1 = c('B'), ident.2 = c('Mono'))
# cd <- FindMarkers(sdat, assay = 'RNA', slot = 'data', logfc.threshold = 0, ident.1 = c('B'), ident.1 = c('CD4 T', 'CD8 T', 'other T'))