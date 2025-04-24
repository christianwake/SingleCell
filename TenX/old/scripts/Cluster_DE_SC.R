library('sys')
library('Seurat')
library('stringr')
library('pheatmap')
library('ggplot2')
#library('reticulate'); use_virtualenv("r-reticulate")
library('umap')
library('textshape')
library('dplyr')
library('biomaRt')
library('grid')
library('scales')

source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

args = commandArgs(trailingOnly=TRUE)
sdat_file <- args[1]
exclude_file <- args[2]
out_rdata <- args[3]
project <- args[4]
username <- args[5]
server <- args[6]

#sdat_file <- '/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/snakemake/results/50Ab/Processed_data.RDS'
#exclude_file <- '/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/snakemake/results/50Ab/Excluded_genes.txt'
# sdat_file <- '/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/snakemake/results/Biolegend/Processed_data.RDS'
# exclude_file <- '/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/snakemake/results/Biolegend/Excluded_genes.txt'

sdat <- readRDS(sdat_file)
### If the file exists and is not empty, read it
if(file.exists(exclude_file) & file.info(exclude_file)$size != 0){
  exclude_gene_names <- read.table(exclude_file)[,1] 
} else{
  exclude_gene_names <- c()
}

sdat <- ScaleData(sdat, assay = "RNA", features = row.names(sdat))

### Exclude the IG and ribosomal RNAs
rna_features <- row.names(sdat@assays$RNA@counts)[which(!row.names(sdat@assays$RNA@counts) %in% exclude_gene_names)]
rna_by_rna <- run_fam(sdat, 'RNA', 'seurat_clusters', 5, features = rna_features)
### Order
rna_by_rna <- rna_by_rna[order(rna_by_rna$p_val),]

# to_keep <- c('sdat', 'rna_by_rna')
# to_remove  <- ls()[which(!ls() %in% to_keep)]
# rm(list = to_remove)

save.image(out_rdata)

if(project != '' & username != '' & server != ''){
  app_dir <- paste0(dirname(sdat_file), '/App/')
  deployApp(appDir = app_dir, appName = project, account = username, server = server)
}
