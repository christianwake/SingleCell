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

source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

if(interactive()){
  project <- '2021614_21-002'
  qc_name <- '2023-Jan'
  out_rdata <-  paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/App/Data.RData')
  username <- 'wakecg'
  server <- 'rstudio-connect.niaid.nih.gov'
  cseq_file <-  paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2', project, '/results/', qc_name, '/Mapped.RDS')
  rds_files <- c('RNA-RNA_clusters', 'RNA-prot_clusters', 'RNA-wnn_clusters', 'prot-RNA_clusters', 'prot-prot_clusters', 'prot-wnn_clusters')
  rds_files <- c('RNA-RNA_clusters', 'RNA-prot_clusters', 'prot-RNA_clusters', 'prot-prot_clusters', 'prot-wnn_clusters')
  rds_files <- c(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Cluster_DE/', rds_files, '.RDS'),
                 paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/SubClusters/predicted.celltype.l1-B.RDS'))
  
}else{
  args = commandArgs(trailingOnly=TRUE)
  out_rdata <- args[1]
  project <- args[2]
  username <- args[3]
  server <- args[4]
  findmarkers <- args[5]
  cseq_file <- args[6]
  rds_files <- args[7:length(args)]
  
}

print('Beginning app data')
### Divide support RDS files into Cluster DE or sub clusters
cluster_de_files <- rds_files[grepl('/Cluster_DE/', rds_files)]
sub_cluster_files <- rds_files[grepl('/SubClusters/', rds_files)]

### Read Seurat data
cseq <- readRDS(cseq_file)

### Save size by removing some unnecessary meta data
cols <- c('MT_sum', 'orig.ident', 'propmt', 'rna_size', 'prot_size', 'ngene', 'bc', 'Assignment_CR', 'Assignment_simple', 'cell_id', 'FCID', 'HTO_index_name', 'Cell_subset', 'Cell_Count',
          'Total_cells', 'cDNA_Concentration', 'Date_of_library_prep', 'Surface_protein_index_name', 'CSO_target_reads', 'GEx_index_name', 'VDJ_index_name', 'VDJ_target_reads', 
          'Date_Sequenced')
cols <- cols[cols %in% colnames(cseq@meta.data)]
for(c in cols){
  cseq[[c]] <- NULL
}

### removing some assays and reductions
DefaultAssay(cseq) <- 'RNA'
umaps <- names(cseq@reductions)[grepl('_umap', names(cseq@reductions))]
cseq <- DietSeurat(cseq, counts = F, dimreducs = umaps, assays = c('RNA', 'prot'))

dat <- lapply(cluster_de_files, function(file) readRDS(file))
names(dat) <- sapply(cluster_de_files, function(file) gsub('.RDS', '', basename(file)))

### Order
rna_by_rna <- dat[['RNA-RNA_clusters']]
rna_by_prot <- dat[['RNA-prot_clusters']]
rna_by_wnn <- dat[['RNA-wnn_clusters']]
prot_by_rna <- dat[['prot-RNA_clusters']]
prot_by_prot <- dat[['prot-prot_clusters']]
prot_by_wnn <- dat[['prot-wnn_clusters']]
rm(dat)

### Subcluster data
dat <- lapply(sub_cluster_files, function(file) readRDS(file))
names(dat) <- sapply(sub_cluster_files, function(file) gsub('.RDS', '', basename(file)))
### Add each dimensionality reduction and meta data column to the real cseq
for(sc in names(dat)){
  ### Add cluster info to the original object
  cseq@meta.data[, colnames(dat[[sc]]@meta.data)] <- dat[[sc]]@meta.data
  ### Add UMAP info to the original object
  for(umap in names(dat[[sc]]@reductions)){
    cseq@reductions[umap] <- dat[[sc]]@reductions[umap] 
  }
}
rm(dat)

# to_keep <- c('cseq', 'rna_by_rna', 'rna_by_prot', 'rna_by_wnn','prot_by_rna', 'prot_by_prot', 'prot_by_wnn')
# to_remove  <- ls()[which(!ls() %in% to_keep)]
# rm(list = to_remove)
total_bytes <- as.numeric(object.size(cseq)[1]) + as.numeric(object.size(rna_by_rna)[1]) + as.numeric(object.size(rna_by_prot)[1]) +
  as.numeric(object.size(rna_by_wnn)[1]) + as.numeric(object.size(prot_by_rna)[1]) + as.numeric(object.size(prot_by_prot)[1]) +
  as.numeric(object.size(prot_by_wnn)[1])
deploy_app <- T
if(total_bytes >= 3145728000){
  deploy_app <- F
  print(paste0('App will not be deployed because it is too large (', total_bytes, ' bytes)'))
}

### Maximum size for deployment is 3145728000 bytes
save.image(out_rdata)

print(project)
print(username)
print(server)
print(deploy_app)
if(project !='' & username != '' & server != '' & deploy_app == T){
  print(paste0('Deploying app titled ', project))
  app_dir <- paste0(dirname(cseq_file), '/App/')
  deployApp(appDir = app_dir, appName = project, account = username, server = server)
}
print('done')
