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
out_rdata <- args[3]
project <- args[4]
username <- args[5]
server <- args[6]

# cseq_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/Go2/Processed_data.RDS'
# exclude_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/Go2/Excluded_genes.txt'

cseq <- readRDS(cseq_file)
DefaultAssay(cseq) <- 'prot'
exclude_gene_names <- read.table(exclude_file)[,1]

cseq <- ScaleData(cseq, assay = "prot", features = row.names(cseq))
### Exclude the IG and ribosomal RNAs
rna_features <- row.names(cseq@assays$RNA@counts)[which(!row.names(cseq@assays$RNA@counts) %in% exclude_gene_names)]
rna_by_rna <- run_fam(cseq, 'RNA', 'RNA_clusters', 5, features = rna_features)
rna_by_prot <- run_fam(cseq, 'RNA', 'prot_clusters', 5, features = rna_features)
rna_by_wnn <- run_fam(cseq, 'RNA', 'wnn_clusters', 5, features = rna_features)
prot_by_rna <- run_fam(cseq, 'prot', 'RNA_clusters', 5)
prot_by_prot <- run_fam(cseq, 'prot', 'prot_clusters', 5)
prot_by_wnn <- run_fam(cseq, 'prot', 'wnn_clusters', 5)
### Order
rna_by_rna <- rna_by_rna[order(rna_by_rna$p_val),]
rna_by_prot <- rna_by_prot[order(rna_by_prot$p_val),]
rna_by_wnn <- rna_by_wnn[order(rna_by_wnn$p_val),]
prot_by_rna <- prot_by_rna[order(prot_by_rna$p_val),]
prot_by_prot <- prot_by_prot[order(prot_by_prot$p_val),]
prot_by_wnn <- prot_by_wnn[order(prot_by_wnn$p_val),]


# to_keep <- c('cseq', 'rna_by_rna', 'rna_by_prot', 'rna_by_wnn','prot_by_rna', 'prot_by_prot', 'prot_by_wnn')
# to_remove  <- ls()[which(!ls() %in% to_keep)]
# rm(list = to_remove)

save.image(out_rdata)

if(project !='' & username != '' & server != ''){
  app_dir <- paste0(dirname(cseq_file), '/App/')
  deployApp(appDir = app_dir, appName = project, account = username, server = server)
}
