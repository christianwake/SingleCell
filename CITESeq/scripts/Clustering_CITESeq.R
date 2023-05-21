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
library('lisi')

source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

if(interactive()){
  # project <- '2021614_21-002'
  # qc_name <- 'Go2'
  # negative_markers <- ''
  
  project <- '2022619_857.3b'
  qc_name <- 'QC_first_pass'
  negative_markers <- ''
  RNA_resolution <- 0.08
  prot_resolution <- 0.05
  sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/DSB_normalized_data.RDS')
  exclude_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Excluded_genes.txt')
  out_rds <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Clustered.RDS')
  out_pdf <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/ADT_clusters.pdf')
}else{
  args = commandArgs(trailingOnly=TRUE)
  sdat_file <- args[1]
  exclude_file <- args[2]
  gtf_file <- args[3]
  gene_file <- args[4]
  RNA_resolution <- as.numeric(args[5])
  prot_resolution <- as.numeric(args[6])
  test_var <- args[7]
  integration_file <- args[8]
  out_rds <- args[9]
  out_pdf <- args[10]
}

### Do direct matching from gtf
if(grepl('\\.gtf', gtf_file)){
  gtf <- read_gtf(gtf_file, feature_type = 'gene', atts_of_interest = c('gene_id', 'gene_name',  'gene_biotype'))
} else if(grepl('\\.RDS', gtf_file)){
  gtf <- readRDS(gtf_file)
}

### Defaults
if(is.na(RNA_resolution) | RNA_resolution == ''){
  RNA_resolution <- 0.8
}
print(paste0('Input RNA resolution:', RNA_resolution))
### Defaults
if(is.na(prot_resolution) | prot_resolution == ''){
  prot_resolution <- 0.8
}
print(paste0('Input prot resolution:', prot_resolution))

min_clusters <- 3

sdat <- readRDS(sdat_file)
DefaultAssay(sdat) <- 'RNA'
print(colnames(sdat@meta.data))

if(file.info(exclude_file)[, 'size'] > 0){
  exclude_gene_names <- read.table(exclude_file)[,1]
} else{
  exclude_gene_names <- c()
}

### Seurat's default normalization
sdat <- NormalizeData(sdat)

sdat <- cluster_RNA(sdat, exclude_gene_names, RNA_resolution)

sdat <- cluster_prot(sdat, prot_resolution)

### Remove isotype controls
#dsb <- dsb[which(!row.names(dsb) %in% isotype_controls), ]

### Original method on DSB github
# p_dist <- dist(t(dsb))
# print('Protein UMAP4')
# p_dist <- as.matrix(p_dist)
# # Cluster using Seurat
# print('Protein find neighbors')
# sdat[["p_dist"]] <- Seurat::FindNeighbors(p_dist)$snn
# print('Protein clustering')
# sdat <- Seurat::FindClusters(sdat, resolution = resolution, graph.name = "p_dist")

ilisi_rna <- compute_lisi(sdat@reductions[['RNA_umap']]@cell.embeddings, sdat@meta.data[, 'RNA_clusters', drop = F], c('RNA_clusters'))
ilisi_prot <- compute_lisi(sdat@reductions[['prot_umap']]@cell.embeddings, sdat@meta.data[, 'prot_clusters', drop = F], c('prot_clusters'))


print('Protein cluster sizes:')
print(table(sdat$prot_clusters))
print('Save seurat object')
saveRDS(sdat, file = out_rds)

### Plot
# calculate the average of each protein separately for each cluster 
# prots = rownames(dsb)
# adt_data = cbind(sdat@meta.data, as.data.frame(t(dsb)))
# adt_plot = adt_data %>% 
#   group_by(prot_clusters) %>% 
#   summarize_at(.vars = prots, .funs = mean) %>% 
#   column_to_rownames("prot_clusters")
# adt_plot <- t(adt_plot)
# 
# # plot a heatmap of the average dsb normalized values for each cluster
# p <- pheatmap::pheatmap(adt_plot, color = viridis::viridis(25, option = "B"), fontsize_row = 6, border_color = NA, 
#                         main = 'Protein averages within protein clusters')
# 
# adt_plot <- as.data.frame(adt_plot)
# adt_plot$colors <- ifelse(row.names(adt_plot) %in% negative_markers, 'red', 'black')
# 
# cols = adt_plot[order(match(row.names(adt_plot), p$gtable$grobs[[5]]$label)), ]$colors
# 
# p$gtable$grobs[[5]]$gp=gpar(col = cols, fontsize = 6)
# 
# ### Using ScaledData
# sdat <- ScaleData(sdat, assay = "prot", features = row.names(sdat@assays$prot@data))
# dsb_relative <- as.matrix(sdat@assays$prot@scale.data)
# prots = rownames(dsb_relative)
# adt_data = cbind(sdat@meta.data, as.data.frame(t(dsb_relative)))
# adt_plot2 = adt_data %>% 
#   group_by(prot_clusters) %>% 
#   summarize_at(.vars = prots, .funs = mean) %>% 
#   column_to_rownames("prot_clusters")
# adt_plot2 <- t(adt_plot2)
# # plot a heatmap of the average dsb normalized values for each cluster
# p2 <- pheatmap::pheatmap(adt_plot2, color = viridis::viridis(25, option = "B"), fontsize_row = 6, border_color = NA, 
#                         main = 'Z-scaled protein averages within protein clusters')
# adt_plot2 <- as.data.frame(adt_plot2)
# adt_plot2$colors <- ifelse(row.names(adt_plot2) %in% c('CD3', 'CD4', 'CD14'), 'red', 'black')
# cols = adt_plot2[order(match(row.names(adt_plot2), p2$gtable$grobs[[5]]$label)), ]$colors
# p2$gtable$grobs[[5]]$gp=gpar(col = cols, fontsize = 6)


pdf(out_pdf)
#p
#plot.new()
#p2
#plot.new()
### RNA UMAP
DimPlot(sdat, group.by = 'RNA_clusters', label = T, reduction = 'RNA_umap') + ggtitle('RNA UMAP, RNA clusters')
DimPlot(sdat, group.by = 'prot_clusters', label = T, reduction = 'RNA_umap') + ggtitle('RNA UMAP, protein clusters')
### Protein UMAP
DimPlot(sdat, group.by = 'prot_clusters', label = T, reduction = 'prot_umap') + ggtitle('Protein UMAP, protein clusters')
DimPlot(sdat, group.by = 'RNA_clusters', label = T, reduction = 'prot_umap') + ggtitle('Protein UMAP, RNA clusters')

FeaturePlot(sdat, 'nCount_RNA', reduction = 'RNA_umap')
FeaturePlot(sdat, 'nFeature_RNA', reduction = 'RNA_umap')

FeaturePlot(sdat, 'nCount_RNA', reduction = 'prot_umap')
FeaturePlot(sdat, 'nFeature_RNA', reduction = 'prot_umap')
# FeaturePlot(sdat, 'fraction_MT')
#DimPlot(sdat, group.by = 'Sample_ID', label = T, reduction = 'umap')
#DimPlot(sdat, group.by = 'CR_ID', label = T, reduction = 'umap')

### Specific gene violin plots
if(gene_file != '' & file.exists(gene_file) & file.info(gene_file)$size != 0){
  genes <- as.data.frame(read_excel(gene_file, sheet = 1))
  colnames(genes) <- c('gene_name', 'gene_id', 'category')
  genes <- complete_gene_table(genes, gtf)
  cats <- unique(genes$category)
  
  p <- expression_plots2(sdat, genes, gtf, features, cats, 'seurat_clusters')
  print(p)
  
  # ### For each category, plot genes together (up to 4 each)
  # for(cat in unique(genes$category)){
  #   ### Gene names of this category
  #   gene_names <- genes[which(genes$category == cat), 'gene_name'][[1]]
  #   ### Split into groups if there are more than 4.
  #   if(length(gene_names) > 4){
  #     groups <- split(gene_names, cut(1:length(gene_names), length(gene_names) %/% 4 + 1))
  #   } else{
  #     groups <- list(gene_names)
  #   }
  #   
  #   title <- ggdraw() + draw_label(cat, fontface = 'bold', x = 0, hjust = 0) + theme(plot.margin = margin(0, 0, 0, 250))
  #   ### For each group (may be only 1) of up to four genes, plot violin plots together
  #   for(gnames in groups){
  #     ### For each gene in the group
  #     vplots <- lapply(gnames, function(gene_name) my_vln_plot(sdat, gene_name, de = NA, groupby = 'seurat_clusters'))
  #     pmain <- plot_grid(plotlist = vplots)
  #     print(plot_grid(title, pmain, ncol = 1 , rel_heights = c(0.1, 1)))
  #     #print(FeaturePlot(sdat, features = gnames))
  #   }
  # }
}
dev.off()
