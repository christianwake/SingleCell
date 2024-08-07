library('sys')
library('Seurat')
library('stringr')
library('pheatmap')
library('ggplot2')
#library('umap')
library('textshape')
library('plyr')
library('dplyr')
library('biomaRt')
library('grid')
library('scales')
library('readxl')
library('cowplot')
library('lisi')
library('data.table')

source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R')
#source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

if(interactive()){
  # project <- '2023600_21-0012'
  # qc_name <- '2023-12-01'
  # negative_markers <- ''
  # gene_file <- 'genes.xlsx'
  # QC_input_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/Test_comparisons.csv')
  # integration_file <- ''
  # reses <- c('RNA-0.6', 'prot-0.6')
  
  project <- '2021614_21-002'
  qc_name <- '2024-01-20'
  negative_markers <- ''
  # project <- '2022620_857.1'
  # qc_name <- '2023-01-09' ### published
  # qc_name <- '2023-05-02' ### A few more cells salvaged
  # test_var <- 'timepoint'
  # integration_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/QC_steps/step4_integration.csv')
  # gene_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/genes.xlsx')
  
  ### 1st cluster
  sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, 
                      '/results/', qc_name, '/Normalized.RDS')
  sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, 
                      '/results/', qc_name, '/DSB_normalized_data.RDS')
  
  out_rds <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Clustered.RDS')
  out_pdf <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/RNA_clusters.pdf')
  ### 2nd cluster
  # sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Filtered.RDS')
  # out_rds <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Filtered_clustered.RDS')
  # out_pdf <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Clusters_filtered.pdf')
  # reses <- c('RNA-0.7')

  # project <- '2022619_857.3b'
  # qc_name <- 'QC_first_pass'
  # negative_markers <- ''
  # sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/DSB_normalized_data.RDS')
  # gene_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/genes.xlsx')
  # test_var <- 'Timepoint'
  # integration_file <- ''
  # out_rds <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Clustered.RDS')
  # out_pdf <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/ADT_clusters.pdf')
  # reses <- c('RNA-0.08', 'prot-0.05')
  
  ### Should be universal
  QC_input_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, 
  '/Test_comparisons.csv')
  gtf_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/gtf.RDS')
  exclude_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Excluded_genes.txt')
  
} else{
  args = commandArgs(trailingOnly=TRUE)
  
  sdat_file <- args[1]
  exclude_file <- args[2]
  gtf_file <- args[3]
  gene_file <- args[4]
  QC_input_file <- args[5]
  integration_file <- args[6]
  out_rds <- args[7]
  out_pdf <- args[8]
  reses <- args[9:length(args)]
}

qc_input <- read.table(QC_input_file, sep = ',', header = T)
tests <- c(qc_input[, 'test_name'], qc_input[, 'strat_name1'], qc_input[, 'strat_name2'])
tests <- unique(tests[which(!is.na(tests) & tests != '')])

downsample_var <- 0.3
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

### res will contain the assay names
res <- sapply(reses, function(x) as.numeric(strsplit(x, '-')[[1]][2]))
names(res) <- sapply(reses, function(x) strsplit(x, '-')[[1]][1])
res <- res[which(names(res) %in% names(sdat@assays))]
print(paste0('Will do clustering for assay(s): ', toString(names(res))))

RNA_resolution <- as.numeric(res['RNA'])
### Defaults
if(is.na(RNA_resolution) | RNA_resolution == ''){
  RNA_resolution <- 0.8
}
print(paste0('Input RNA resolution:', RNA_resolution))
min_clusters <- 3

print(sdat_file)
print(QC_input_file)
print(integration_file)
print(out_rds)

### Do direct matching from gtf
if(grepl('\\.gtf', gtf_file)){
  gtf <- read_gtf(gtf_file, feature_type = 'gene', atts_of_interest = c('gene_id', 'gene_name',  'gene_biotype'))
} else if(grepl('\\.RDS', gtf_file)){
  gtf <- readRDS(gtf_file)
}

### If the file exists and is not empty, read it
if(file.exists(exclude_file) & file.info(exclude_file)$size != 0){
  exclude_gene_names <- read.table(exclude_file)[,1]
} else{
  exclude_gene_names <- c()
}

batch_features <- ''
### Option 1. No integration was done.
if(integration_file == '' | is.na(integration_file) | !file.exists(integration_file)){
  print('Clustering on unintegrated data')
 
  sdat <- cluster_RNA(sdat, exclude_gene_names, RNA_resolution, min_clusters) ### Already renames clusters
  umap_name <- 'RNA_umap'
  plot_title <- 'Clustering (no integration)'
} else{
  qc_input <- read.table(integration_file, sep = ',', header = T)
  ### Use only the first line, even if other methods of integration were done
  method <- qc_input[1, 'method']
  batch_features <- unique(qc_input[, 'value'])
  if(toupper(method) == 'HARMONY'){
    print('Clustering on Harmony integrated data')
    ### For Harmony integration
    sdat <- FindNeighbors(sdat, assay = 'RNA', reduction =  'harmony', dims = 1:25)
    if('seurat_clusters' %in% colnames(sdat@meta.data)){
      sdat$seurat_clusters <- NULL
    }
    ### If n clusters is less than the minimum, adjust resolution and try again 
    nclust <- 0
    while(nclust < min_clusters){
      print(paste0(nclust, ' from resolution of ', RNA_resolution, ' is too few.'))
      sdat <- FindClusters(sdat, resolution = RNA_resolution, graph.name = 'RNA_snn')
      sdat$RNA_clusters <- sdat$seurat_clusters
      sdat$seurat_clusters <- NULL
      RNA_resolution <- RNA_resolution + 0.05
      nclust <- length(unique(sdat$RNA_clusters))
    }
    #sdat$clusters_harmony_integrated <- sdat$seurat_clusters
    sdat <- RunUMAP(sdat, dims = 1:26, reduction =  'harmony', assay = 'RNA', reduction.name = 'umap_harmony', reduction.key = 'umap_harmony_')
    umap_name <- 'umap_harmony'
    plot_title <- 'Clustering (harmony integration)'
    #DimPlot(sdat, group.by = 'seurat_clusters', label = T, reduction = 'umap_harmony') + ggtitle(paste0(feat, ', ', red))
  } else if(toupper(method) == 'SEURAT'){
    print('Clustering on Seurat integrated data')
    ### For Seurat integration
    sdat <- ScaleData(sdat, assay = 'integrated')
    sdat <- RunPCA(sdat, npcs = 30, assay = 'integrated', reduction.name = 'pca_seurat')
    sdat <- RunUMAP(sdat, reduction = "pca_seurat", dims = 1:30, assay = 'integrated', reduction.name = 'umap_seurat', reduction.key = 'umap_seurat')
    umap_name <- 'umap_seurat'
    sdat <- FindNeighbors(sdat, reduction = "pca_seurat", dims = 1:30, assay = 'integrated')
    nclust <- 0
    while(nclust < min_clusters){
      print(paste0(nclust, ' from resolution of ', RNA_resolution, ' is too few.'))
      sdat <- FindClusters(sdat, resolution = RNA_resolution, graph.name = 'integrated_snn')
      ### Rename clusters
      sdat$RNA_clusters <- sdat$seurat_clusters
      sdat$seurat_clusters <- NULL
      RNA_resolution <- RNA_resolution + 0.05
      nclust <- length(unique(sdat$RNA_clusters))
    }
    
    #sdat$clusters_seurat_integrated <- sdat$seurat_clusters
    plot_title <- 'Clustering (seurat integration)'
  }
}

ilisi_rna <- compute_lisi(sdat@reductions[[umap_name]]@cell.embeddings, sdat@meta.data[, 'RNA_clusters', drop = F], c('RNA_clusters'))

### Check if protein assay is included
if('prot' %in% names(res) ){
  DefaultAssay(sdat) <- 'prot'
  print('Doing protein clustering')
  ### UMAP, Find neighbors, cluster on sdat assay prot, slot 'data'
  sdat <- cluster_prot(sdat, res['prot'])
  ilisi_prot <- compute_lisi(sdat@reductions[['prot_umap']]@cell.embeddings, sdat@meta.data[, 'prot_clusters', drop = F], c('prot_clusters'))

  features_list <- row.names(sdat@assays$prot@data)
  p <- DotPlot(sdat, features = features_list, group.by = 'prot_clusters', scale.by = 'radius',
               scale = T, dot.scale = 3) +
    RotatedAxis() + 
    scale_x_discrete(labels = 'Protein markers') + labs(y = 'prot_clusters') + ggtitle('Protein markers/clusters')
  DefaultAssay(sdat) <- 'RNA'
  
}

print('save')
saveRDS(sdat, file = out_rds)
#sdat<- readRDS(out_rds)

### Get name of 'umap' reductions
reductions <- names(sdat@reductions)[grepl('umap', names(sdat@reductions))]
### Get names of column names that contain 'cluster'
clusters <- colnames(sdat@meta.data)[grepl('cluster', colnames(sdat@meta.data))]
### And don't begin with 'Subset'
clusters <- clusters[!grepl('^Subset', clusters)]

print(reductions)
print(clusters)
### Name them on the first '_' section, e.g. 'RNA_cluster' named 'RNA'
names(clusters) <- sapply(clusters, function(x) strsplit(x, '_')[[1]][1])
### Rename resolution array names
print('test1')
print(names(res))
names(res) <- clusters[names(res)]

print(clusters)

pdf(out_pdf)
print(p)
### Get quantification to be visualized on each reduction with FeaturePlot
quants <- c(paste0(c('nCount_', 'nFeature_'), rep(names(res), each = 2)), 'fraction_MT')
quants <- quants[quants %in% colnames(sdat@meta.data)]
features <- c('seurat_clusters', clusters, 'Sample_ID', 'Sample_Name', 'CR_ID', 'Time', 'Celltype', tests, batch_features)
features <- unique(features[features %in% colnames(sdat@meta.data)])
for(red in reductions){
  print(red)
  for(quant in quants){
    print(FeaturePlot(sdat, features = quant, reduction = red) + ggtitle(paste0(quant, ', ', red)))
  }
  # print(FeaturePlot(sdat, 'nCount_RNA', reduction = red))
  # print(FeaturePlot(sdat, 'nFeature_RNA', reduction = red))
  # print(FeaturePlot(sdat, 'fraction_MT', reduction = red))
  for(feat in features){
    print(feat)
    subtitle <- ''
    if(feat %in% names(res)){
      subtitle <- paste0('Resolution: ', res[feat])
    }
    p <- DimPlot(sdat, group.by = feat, label = T, reduction = red) + 
            ggtitle(paste0(feat, ', ', red), subtitle =  subtitle)
    #### Hide legend if it will make the plot illegibile
    if(length(unique(sdat@meta.data[, feat])) > 15){
      p <- p+ theme(legend.position = "none")      
    }
    print(p)
  }
}
#dev.off()
### Plot
# rna <- as.matrix(sdat@assays$RNA@data)
# calculate the average of each protein separately for each cluster 
# rnas = rownames(rna)
# adt_data = cbind(sdat@meta.data, as.data.frame(t(rna)))
# adt_plot = adt_data %>% 
#   group_by(seurat_clusters) %>% 
#   summarize_at(.vars = rnas, .funs = mean) %>% 
#   column_to_rownames("seurat_clusters")
# adt_plot <- t(adt_plot)
# 
# # plot a heatmap of the average rna values for each cluster
# p <- pheatmap::pheatmap(adt_plot, color = viridis::viridis(25, option = "B"), fontsize_row = 6, border_color = NA, 
#                         main = 'count averages within clusters')
# 
# adt_plot <- as.data.frame(adt_plot)
# adt_plot$colors <- ifelse(row.names(adt_plot) %in% negative_markers, 'red', 'black')
# 
# cols = adt_plot[order(match(row.names(adt_plot), p$gtable$grobs[[5]]$label)), ]$colors
# 
# p$gtable$grobs[[5]]$gp=gpar(col=cols, fontsize = 6)
# 
# ### Using ScaledData
# sdat <- ScaleData(sdat, assay = "RNA", features = row.names(sdat@assays$RNA@data))
# rna_relative <- as.matrix(sdat@assays$RNA@scale.data)
# rnas = rownames(rna_relative)
# adt_data = cbind(sdat@meta.data, as.data.frame(t(rna_relative)))
# adt_plot2 = adt_data %>% 
#   group_by(seurat_clusters) %>% 
#   summarize_at(.vars = rnas, .funs = mean) %>% 
#   column_to_rownames("seurat_clusters")
# adt_plot2 <- t(adt_plot2)
# # plot a heatmap of the average rna normalized values for each cluster
# p2 <- pheatmap::pheatmap(adt_plot2, color = viridis::viridis(25, option = "B"), fontsize_row = 6, border_color = NA, 
#                         main = 'Z-scaled count averages within clusters')
# adt_plot2 <- as.data.frame(adt_plot2)
# adt_plot2$colors <- ifelse(row.names(adt_plot2) %in% c('CD3', 'CD4', 'CD14'), 'red', 'black')
# cols = adt_plot2[order(match(row.names(adt_plot2), p2$gtable$grobs[[5]]$label)), ]$colors
# p2$gtable$grobs[[5]]$gp=gpar(col=cols, fontsize = 6)

### Specific gene expression dot plots
if(gene_file != '' & file.exists(gene_file) & file.info(gene_file)$size != 0){
  print(paste0('Violin plots for: ', gene_file))
  genes <- as.data.frame(read_excel(gene_file, sheet = 1))
  colnames(genes) <- c('gene_name', 'gene_id', 'category')
  ### Does not currently allow duplicate gene names
  genes <- complete_gene_table(genes, gtf)
  cats <- unique(genes$category)
  
  for(clust in clusters){
    p <- expression_plots2(sdat, genes, gtf, features, cats, group_by = clust)
    print(p)
  }

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

if(length(reductions) > 1){
  print(paste0('doing ilisi to compare ', toString(reductions[1:2])))
  plot_ilisi(sdat, reductions, 'Sample_ID')
}


dev.off()
