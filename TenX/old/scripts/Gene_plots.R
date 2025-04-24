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
library('readxl')
library('cowplot')

source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

args = commandArgs(trailingOnly=TRUE)

sdat_file <- args[1]
gene_file <- args[2]
out_pdf <- args[3]

# sdat_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021618_galt/results/RNA_cell_filtered.RDS'
# gene_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021618_galt/genes.xlsx'
# out_pdf <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021618_galt/RNA_clusters.pdf'

sdat <- readRDS(sdat_file)
DefaultAssay(sdat) <- 'RNA'

pdf(out_pdf)
### Specific gene violin plots
if(gene_file != '' & file.exists(gene_file) & file.info(gene_file)$size != 0){
  genes <- read_excel(gene_file, sheet = 1)
  colnames(genes) <- c('gene_name', 'category')
  genes <- genes[which(genes$gene_name %in% row.names(sdat)),]
    
  ### For each category, plot genes together (up to 4 each)
  for(cat in unique(genes$category)){
    ### Gene names of this category
    gene_names <- genes[which(genes$category == cat), 'gene_name'][[1]]
    ### Split into groups if there are more than 4.
    if(length(gene_names) > 4){
      groups <- split(gene_names, cut(1:length(gene_names), length(gene_names) %/% 4 + 1))
    } else{
      groups <- list(gene_names)
    }
    
    title <- ggdraw() + draw_label(cat, fontface = 'bold', x = 0, hjust = 0) + theme(plot.margin = margin(0, 0, 0, 250))
    ### For each group (may be only 1) of up to four genes, plot violin plots together
    for(gnames in groups){
      ### For each gene in the group
      vplots <- lapply(gnames, function(gene_name) my_vln_plot(sdat, gene_name, de = NA, genes, groupby = 'seurat_clusters'))
      pmain <- plot_grid(plotlist = vplots)
      print(plot_grid(title, pmain, ncol = 1 , rel_heights = c(0.1, 1)))
      #print(FeaturePlot(sdat, features = gnames))
    }
  }
}

dev.off()
