library('sys')
library('Seurat')
library('SeuratWrappers')
library('harmony')
library('stringr')
library('pheatmap')
library('ggplot2')
#library('umap')
library('textshape')
library('dplyr')
library('biomaRt')
library('grid')
library('scales')
library('vegan')
library('data.table')
library('cowplot')

source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R')

run_harmony_integration <- function(sdat, batch_feature, exclude_gene_names){
  ### NormalizeData uses 'count' slot and overrides the 'data' slot
  #sdat <- NormalizeData(sdat, assay = 'RNA')
  
  sdat <- FindVariableFeatures(sdat, assay = 'RNA', selection.method = "vst")
  #VariableFeaturePlot(sdat)
  variable_genes <- sdat@assays$RNA@var.features
  variable_genes <- variable_genes[which(!variable_genes %in% exclude_gene_names)]
  
  ### data slot is unchanged, scaled.data slot is added
  ### By default, it seems that ScaleData wants the raw 'counts', so if it is absent it will use 'data' but throw a warning
  sdat <- ScaleData(sdat, assay = "RNA", vars.to.regress = c('fraction_MT', 'nFeature_RNA'))
  sdat <- RunPCA(sdat, verbose = F, assay = "RNA", features = variable_genes)
  
  ### Results go in sdat@reductions$harmony
  sdat <- RunHarmony(sdat, batch_feature, plot_convergence = TRUE, assay.use = 'RNA')
  return(sdat)
}

run_seurat_integration <- function(sdat, batch_feature, exclude_gene_names, norm_method){
  ### Seurat's integration. This ends with data in new slot "integrated"
  # split the dataset into a list of two seurat objects (stim and CTRL)
  sdat_list <- SplitObject(sdat, split.by = batch_feature)
  
  ### Normalize for each dataset independently
  if(norm_method == 'Seurat'){### Standard Seurat log normalization, should be used on cellranger-processed data like from TenX
    ### NormalizeData uses 'count' slot and overrides the 'data' slot
    sdat_list <- lapply(X = sdat_list, FUN = function(x) {
      x <- NormalizeData(x)
    }) 
  } else if(norm_method == 'TPM'){
    sdat_list <- lapply(X = sdat_list, FUN = function(x) {
      x@assays$RNA@data <- as.matrix(log(x@assays$RNA@counts + 1))
    })
  }
  ### Identify variable features for each dataset independently
  sdat_list <- lapply(X = sdat_list, FUN = function(x) {
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })
  
  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = sdat_list)
  features <- features[which(!features %in% exclude_gene_names)]
  anchors <- FindIntegrationAnchors(object.list = sdat_list, anchor.features = features)
  # this command creates an 'integrated' data assay
  anc <- IntegrateData(anchorset = anchors)
  sdat@assays$integrated <- anc@assays$integrated
  sdat$nCount_integrated <- colSums(sdat@assays$integrated@data)
  sdat$nFeature_integrated <- colSums(sdat@assays$integrated@data > 0 )
  return(sdat)
}

run_integration <- function(sdat, integration_method, batch_feature, exclude_gene_names = c(), norm_method = 'Seurat'){
  print(paste0('Integrating over ', batch_feature, ' with method ', integration_method))
  if(toupper(integration_method) == 'HARMONY'){
    print('Harmony integration')
    sdat <- run_harmony_integration(sdat, batch_feature, exclude_gene_names)
  } else if(toupper(integration_method) == 'SEURAT'){
    print('Seurat integration')
    sdat <- run_seurat_integration(sdat, batch_feature, exclude_gene_names, norm_method)
  }
  return(sdat)
}

if(interactive()){
  #project <- '2021617_misc-c'
  #qc_name <- 'Filter1'
  project <- '2022620_857.1'
  qc_name <- '2022-12-09'
  sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/PostQC3.RDS')
  QC_input_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/QC_steps/step4_integration.csv')
  exclude_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/Excluded_genes.txt')
  norm_method <- 'TPM'
  out_rds <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/PostQC4.RDS')
  out_pdf <-  paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/Integration.pdf')
} else{
  args = commandArgs(trailingOnly=TRUE)
  sdat_file <- args[1]
  QC_input_file <- args[2]
  exclude_file <- args[3]
  norm_method <- args[4]
  out_rds <- args[5]
  pdf_file <- args[6]
}

sdat <- readRDS(sdat_file)
DefaultAssay(sdat) <- 'RNA'


### If the file exists and is not empty, read it
if(file.exists(exclude_file) & file.info(exclude_file)$size != 0){
  exclude_gene_names <- read.table(exclude_file)[,1] 
} else{
  exclude_gene_names <- c()
}

qc_input <- read.table(QC_input_file, sep = ',', header = T)

### For each row, run imputation of that method and value
for(r in row.names(qc_input)){
  sdat <- run_integration(sdat, qc_input[r, 'method'], qc_input[r, 'value'], exclude_gene_names, norm_method)
}

saveRDS(sdat, file = out_rds)

# features <- c('nFeature_RNA', 'nCount_RNA', 'nFeature_integrated', 'nCount_integrated')
# features <- features[features %in% colnames(sdat@meta.data)]
features <- colnames(sdat@meta.data)[grepl('nFeature_', colnames(sdat@meta.data))]
features <- c(features, colnames(sdat@meta.data)[grepl('nCount_', colnames(sdat@meta.data))])
reductions <- names(sdat@reductions)
pdf(pdf_file)
### For each batch_feature entered (probaly only 1)
for(batch_feature in unique(qc_input[,2])){
  ### For each feature that we can plot
  for(feature in features){
    print(VlnPlot(object = sdat, features = feature, group.by = batch_feature, pt.size = 0))
  }
  ### For each reduction that we can plot
  for(red in reductions){
    print(DimPlot(object = sdat, reduction = red, pt.size = 0, group.by = batch_feature))
  } 
}
dev.off()

