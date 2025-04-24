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

source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

args = commandArgs(trailingOnly=TRUE)
cseq_file <- args[1]
exclude_file <- args[2]
out_rds <- args[3]
assay <- args[4]
cluster <- args[5]
#cseq_file <- '/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/snakemake/results/50Ab/Processed_data.RDS'
#exclude_file <- '/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/snakemake/results/50Ab/Excluded_genes.txt'
# cseq_file <- '/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/snakemake/results/Biolegend/Processed_data.RDS'
# exclude_file <- '/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/snakemake/results/Biolegend/Excluded_genes.txt'
# assay = 'RNA'
# cluster = 'RNA_clusters'

cseq <- readRDS(cseq_file)
DefaultAssay(cseq) <- 'prot'
exclude_gene_names <- read.table(exclude_file)[,1]

#cseq <- ScaleData(cseq, assay = "prot", features = row.names(cseq))
### Exclude the IG and ribosomal RNAs
if(assay == 'RNA'){
  features <- row.names(cseq@assays$RNA@counts)[which(!row.names(cseq@assays$RNA@counts) %in% exclude_gene_names)]
} else{
  features <- NULL
}
print(table(cseq[[cluster]]))
de_dat <- run_fam(cseq, assay, cluster, 5, features = features, pdf_file = NA)
### Order
de_dat <- de_dat[order(de_dat$p_val),]
saveRDS(de_dat, file = out_rds)

