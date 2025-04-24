library('sys')
library('Seurat')
library('stringr')
library('pheatmap')
library('ggplot2')
#library('umap')
#library('textshape')
library('dplyr')
#library('biomaRt')
library('grid')
library('scales')
library('rsconnect')
library('readxl')
library('data.table')
library('plyr')
library('lisi')

source('/data/vrc_his/douek_lab/snakemakes/sc_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/Utility_functions.R')
source('/data/vrc_his/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

if(interactive()){
  project <- '2021614_21-002'
  qc_name <- 'DSB_by_sample'
  out_rdata <-  paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/App/Data.RData')
  username <- 'wakecg'
  server <- 'rstudio-connect.niaid.nih.gov'
  server <- ''
  cseq_file <-  paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Mapped.RDS')
  rds_files <- c('RNA-RNA_clusters', 'RNA-prot_clusters', 'RNA-wnn_clusters', 'prot-RNA_clusters', 'prot-prot_clusters', 'prot-wnn_clusters')
  rds_files <- c('RNA-RNA_clusters', 'RNA-prot_clusters', 'prot-RNA_clusters', 'prot-prot_clusters', 'prot-wnn_clusters')
  rds_files <- c(paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Cluster_DE/One-all/', 
                        rds_files, '.RDS'),
                 paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/SubClusters/predicted.celltype.l1-B.RDS'))
  
  gene_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/genes.xlsx')
  pdf_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/SomeGenes.pdf')
  gtf_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/data/gtf.RDS')
  
}else{
  args = commandArgs(trailingOnly=TRUE)
  out_rdata <- args[1]
  pdf_file <- args[2]
  gene_file <- args[3]
  project <- args[4]
  username <- args[5]
  server <- args[6]
  gtf_file <- args[7]
  findmarkers <- args[8]
  cseq_file <- args[9]
  rds_files <- args[10:length(args)]
  
}

print('Beginning app data')
### Divide support RDS files into Cluster DE or sub clusters
cluster_de_files <- rds_files[grepl('/Cluster_DE/', rds_files)]
sub_cluster_files <- rds_files[grepl('/SubClusters/', rds_files)]

downsample_var <- 0.3
if(!interactive()){
  ### Read Seurat data
  cseq <- readRDS(cseq_file)
  print('Done reading Seurat object')
  ### Downsample so I can work interactiveley for testing
  csub <- DietSeurat(cseq, layers = c('data'), assays = c('RNA', 'prot'))
  cells <- sample(x = colnames(csub), size = (length(colnames(csub)) * downsample_var), replace = F)
  
  csub <- subset(csub, cells = cells)
  saveRDS(csub, file = gsub('.RDS', paste0('_DownSampledTo', downsample_var, '.RDS'), cseq_file))
  print('saved downsampled version for testing')
} else{
  cseq <- readRDS(gsub('.RDS', paste0('_DownSampledTo', downsample_var, '.RDS'), cseq_file))
}

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

print("Diet Seurat")
cseq <- DietSeurat(cseq, layers = c('data'), dimreducs = umaps, assays = c('RNA', 'prot'))
print('Reading cluster DE files')
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

### Specific gene expression dot plots
if(pdf_file != '' & gene_file != '' & file.exists(gene_file) & file.info(gene_file)$size != 0){
  ### Do direct matching from gtf
  if(grepl('\\.gtf', gtf_file)){
    gtf <- read_gtf(gtf_file, feature_type = 'gene', atts_of_interest = c('gene_id', 'gene_name',  'gene_biotype'))
  } else if(grepl('\\.RDS', gtf_file)){
    gtf <- readRDS(gtf_file)
  }
  
  print(paste0('Violin plots for: ', gene_file))
  genes <- as.data.frame(read_excel(gene_file, sheet = 1))
  colnames(genes) <- c('gene_name', 'gene_id', 'category')
  ### Does not currently allow duplicate gene names
  genes <- complete_gene_table(genes, gtf)
  cats <- unique(genes$category)
  
  ### Get names of column names that contain 'cluster'
  clusters <- colnames(cseq@meta.data)[grepl('cluster', colnames(cseq@meta.data))]
  ### And don't begin with 'Subset'
  clusters <- clusters[!grepl('^Subset', clusters)]
  pdf(pdf_file)
  for(clust in clusters){
    print(clust)
    p <- expression_plots2(cseq, genes, gtf, features, cats, group_by = clust)
    print(p)
  }

  if(length(umaps) > 1){
    print(paste0('doing ilisi to compare ', toString(umaps[1:2])))
    plot_ilisi(cseq, umaps, 'Sample_ID')
  }
  dev.off()
}




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
