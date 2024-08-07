library('sys')
library('Seurat')
library('SeuratDisk')
library('stringr')
library('pheatmap')
library('ggplot2')
library('umap')
library('textshape')
library('dplyr')
library('biomaRt')
library('grid')
library('scales')
library('data.table')
library('SingleR')
library('scuttle')
library('scRNAseq')

source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R')

if(interactive()){
  # project <- '2021617_mis-c'
  # qc_name <- 'Dropout_mitigated'
  # species <- 'hsapiens'
  # cell_type <- 'PBMC'
  # method <- 'Seurat'
  # 
  # project <- '2021600_kristin'
  # qc_name <- 'Run2022'
  # species <- 'hsapiens'
  # cell_type <- 'T'
  # method <- 'SingleR'
  
  project <- '2021614_21-002'
  qc_name <- '2024-01-20'
  method <- 'SingleR'
  #reference <- 'celldex:HumanPrimaryCellAtlasData'
  reference_name <- 'celldex:ImmGenData'
  sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Cluster_filtered.RDS')
  out_rds <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Mapped.RDS')
  out_pdf <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Mapping.pdf')
  gtf_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/gtf.RDS')
  exclude_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Excluded_genes.txt')
}else{
  args = commandArgs(trailingOnly=TRUE)
  sdat_file <- args[1]
  exclude_file <- args[2]
  gtf_file <- args[3]
  reference_name <- args[4]
  method <- args[5]
  species <- args[6]
  cell_type <- args[7]
  out_rds <- args[8]
  #out_pdf <- args[6]
}

### Do direct matching from gtf
if(grepl('\\.gtf', gtf_file)){
  gtf <- read_gtf(gtf_file, feature_type = 'gene', atts_of_interest = c('gene_id', 'gene_name',  'gene_biotype'))
} else if(grepl('\\.RDS', gtf_file)){
  gtf <- readRDS(gtf_file)
}

#sdat <- readRDS(sdat_file)
downsample_var <- 0.2
if(!interactive()){
  ### Read Seurat data
  sdat <- readRDS(sdat_file)
  print('Done reading Seurat object')
  ### Downsample so I can work interactiveley for testing
  ssub <- DietSeurat(sdat, counts = F, assays = c('RNA', 'prot'))
  cells <- sample(x = colnames(ssub), size = (length(colnames(ssub)) * downsample_var), replace = F)
  
  ssub <- subset(ssub, cells = cells)
  saveRDS(ssub, file = gsub('.RDS', paste0('_DownSampledTo', downsample_var, '.RDS'), sdat_file))
  print('saved downsampled version for testing')
} else{
  sdat <- readRDS(gsub('.RDS', paste0('_DownSampledTo', downsample_var, '.RDS'), sdat_file))
}
DefaultAssay(sdat) <- 'RNA'

### the Emily hack to fix the weird named list thing
colnames(sdat@assays$RNA@data) <- rownames(sdat@meta.data)
if(dim(sdat@assays$RNA@counts)[1] > 0){
  colnames(sdat@assays$RNA@counts) <- rownames(sdat@meta.data)
}

### If the file exists and is not empty, read it
if(file.exists(exclude_file) & file.info(exclude_file)$size != 0){
  exclude_gene_names <- read.table(exclude_file)[,1]
} else{
  exclude_gene_names <- c()
}
features <- row.names(sdat)[which(!(row.names(sdat) %in% exclude_gene_names))]

#sdat <- DietSeurat(sdat, counts = T, assays = c('RNA', 'RNA_orig'))
# sdat <- subset(x = sdat, downsample = 1000)
# print(dim(sdat))
# saveRDS(sdat, file = '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021617_mis-c/results/Dropout_mitigated/Test_suberset.RDS')

# if(method == 'Seurat'){
#   reference <- LoadH5Seurat("/hpcdata/vrc/vrc1_data/douek_lab/reference_sets/pbmc/pbmc_multimodal.h5seurat")
#   method <- 'Seurat'
# }
### Single R with celldex HPCA reference
if(method == 'SingleR'){
  ### If celldex, grab reference names. Options are ImmGen, MouseRNAseqData, BlueprintEncodeData, HumanPrimaryCellAtlasData, DatabaseImmuneCellExpressionData, NovershternHematopoieticData, MonacoImmuneData
  if(grepl('CELLDEX', toupper(reference_name))){
    library('celldex')
    reference <- strsplit(reference_name, ':')[[1]][2]
    print(paste0('Using Celldex ', reference))
    ### Get the celldex fetch function but named as fetch_func
    fetch_func <- match.fun(reference)
    reference <- fetch_func()
  }
  row.names(reference) <- toupper(row.names(reference)) 
  genes <- row.names(reference)
}

#res <- gene_matching_v1(row.names(sdat), genes, gtf)
### For converting from gene_name to gene_id
res <- gene_matching_v1(genes, features, gtf, one2multi = 'exclude')
### Update reference to match sdat (for those that do)
row.names(reference) <- res[row.names(reference), 'gene_id']

reference <- reference[which(!is.na(row.names(reference))),]

if(toupper(method) == 'SEURAT'){
  reference <- LoadH5Seurat(reference)
  sdat <-  seurat_annotation(sdat, reference)
} else if(toupper(method) == 'SINGLER'){
  ### Convert from Seurat object to SingleCellExperiment
  dat <- as.SingleCellExperiment(sdat)
  ### SingleR() expects reference datasets to be normalized and log-transformed.
  ### Requires RNA 'counts' assay data
  dat <- logNormCounts(dat) 
  
  sdat_pred <- SingleR(test = dat, ref = reference, assay.type.test = 1, labels = reference$label.main)
  table(sdat_pred$labels)
  sdat@meta.data[row.names(sdat_pred), 'predicted_label_main'] <- sdat_pred$labels
  
  sdat_pred <- SingleR(test = dat, ref = reference, assay.type.test = 1, labels = reference$label.fine)
  table(sdat_pred$labels)
  sdat@meta.data[row.names(sdat_pred), 'predicted_label_fine'] <- sdat_pred$labels
}

# red <- 'umap'
# DimPlot(sdat, group.by = 'predicted_label_main', label = T, reduction = red) + ggtitle(paste0('predicted_label_main', ', ', red))

### predicted.celltype.l1 and 2 contain the predicted cell type from the reference
### predicted.celltype.l1.score and 2.score contain the ... maximum score???
### 
# Idents(sdat) <- 'predicted_label_fine'
# treg_markers <- FindMarkers(sdat, ident.1 = "Treg", only.pos = TRUE, logfc.threshold = 0.1)
# print(head(treg_markers))
saveRDS(sdat, out_rds)
# 
# pdf(out_pdf)
# DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
# ### ref.umap is defined by the reference
# DimPlot(sdat, reduction = "ref.umap", group.by = "predicted_label_main", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
# DimPlot(sdat, reduction = "ref.umap", group.by = "predicted_label_fine", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
# 
# ### Cell type scores, to compare to the cell type declarations
# FeaturePlot(sdat, features = c("pDC", "CD16 Mono", "Treg"),  reduction = "ref.umap", cols = c("lightgrey", "darkred"), ncol = 3) & theme(plot.title = element_text(size = 10))
# 
# VlnPlot(sdat, features = c("CLEC4C", "LILRA4"), sort = TRUE) + NoLegend()
# 
# DefaultAssay(sdat) <- 'predicted_ADT'
# # see a list of proteins: rownames(sdat)
# FeaturePlot(sdat, features = c("CD3-1", "CD45RA", "IgD"), reduction = "ref.umap", cols = c("lightgrey", "darkgreen"), ncol = 3)
# 
# 
# ### "De novo" visualization, to catch if there are cell types in query that are not within reference
# #merge reference and query
# reference$id <- 'reference'
# sdat$id <- 'query'
# refquery <- merge(reference, sdat)
# refquery[["spca"]] <- merge(reference[["spca"]], sdat[["ref.spca"]])
# refquery <- RunUMAP(refquery, reduction = 'spca', dims = 1:50)
# DimPlot(refquery, group.by = 'id', shuffle = TRUE)
# 
# dev.off()

### To create a test set for interactive use
# csub <- subset(x = sdat, downsample = 10000)
# saveRDS(csub, file = '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/Go2/Test_subset.RDS')

### To create a test set for interactive use
# csub <- subset(x = sdat, downsample = 1000)
# print(dim(csub))
# saveRDS(csub, file = '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/Go2/Test_suberset.RDS')
# 
# cols <- c('MT_sum', 'orig.ident', 'propmt', 'rna_size', 'ngene', 'bc')
# for(c in cols){
#   sdat[[c]] <- NULL
# }
