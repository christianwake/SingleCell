library('sys')
library('Seurat')
library('stringr')
library('ggplot2')
#library('textshape')
library('plyr')
library('dplyr')
library('cowplot')

source('/data/vrc_his/douek_lab/wakecg/sc_functions.R')
source('/data/vrc_his/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

if(interactive()){
  # sdat_file <- '/data/vrc_his/douek_lab/projects/RNASeq/2021617_mis-c/results/Dropout_mitigated/PostQC3.RDS'
  # gene_file <- '/data/vrc_his/douek_lab/projects/RNASeq/2021617_mis-c/genes.xlsx'
  # test_var <- 'Condition'
  # out_pdf <- '/data/vrc_his/douek_lab/projects/RNASeq/2021617_mis-c/results/RNA_clusters.pdf'
  
  sdat_file <- '/data/vrc_his/douek_lab/projects/RNASeq/2021600_kristin/results/Run2022/PostQC3.RDS'
  gene_file <- '/data/vrc_his/douek_lab/projects/RNASeq/2021600_kristin/220921_transcriptome_analysis.xlsx'
  name_match_file <- '/data/vrc_his/douek_lab/wakecg/data/Gene_names.RDS'
  test_var <- 'Time,TCR.status'
  out_pdf <- '/data/vrc_his/douek_lab/projects/RNASeq/2021600_kristin/results/Run2022/Violins.pdf'
} else{
  args = commandArgs(trailingOnly=TRUE)
  
  sdat_file <- args[1]
  gene_file <- args[2]
  name_match_file <- args[3]
  test_var <- args[4]
  out_pdf <- args[5]
}

sdat <- readRDS(sdat_file)
DefaultAssay(sdat) <- 'RNA'

### Load gtf info
dat <- readRDS(name_match_file)
dat <- dat[which(!is.na(dat$plot_name)), ]
key <- sapply(dat$plot_name, function(x) dat[which(dat$plot_name == x), 'gtf_id'])

pdf(out_pdf)
### Specific gene violin plots
if(gene_file != '' & file.exists(gene_file) & file.info(gene_file)$size != 0){
  genes <- read_excel(gene_file, sheet = 1)
  colnames(genes) <- c('gene_name', 'category')
  #genes <- genes[which(genes$gene_name %in% row.names(sdat)),]
  genes <- genes[which(genes$gene_name %in% names(key)),]
  
  ### For each category, plot genes together (up to 4 each)
  for(cat in unique(genes$category)){
    ### Gene names of this category
    gene_names <- key[genes[which(genes$category == cat), 'gene_name'][[1]]]
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
      vplots <- lapply(1:length(gnames), function(i) my_vln_plot(sdat, gnames[i], title = names(gnames)[i], de = NA, groupby = 'seurat_clusters'))
      pmain <- plot_grid(plotlist = vplots)
      print(plot_grid(title, pmain, ncol = 1 , rel_heights = c(0.1, 1)))
      #print(FeaturePlot(sdat, features = gnames))
    }
  }
}

dev.off()
