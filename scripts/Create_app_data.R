library('sys')
library('Seurat')
library('stringr')
library('pheatmap')
library('ggplot2')
library('umap')
#library('textshape')
library('dplyr')
library('biomaRt')
library('grid')
library('scales')
library('rsconnect')
library('readxl')
library('cowplot')

source('/data/vrc_his/douek_lab/snakemakes/sc_functions.R')
source('/data/vrc_his/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
source('/data/vrc_his/douek_lab/snakemakes/Utility_functions.R')

if(interactive()){
  # project <- '2021618_galt'
  # qc_name <- 'Go1'
  # additional <- ''
  # 
  project <- '2021600_kristin'
  qc_name <- 'Run2023-05-14'
  additional <- '/data/vrc_his/douek_lab/projects/RNASeq/2021600_kristin/SampleSheets/Additional.csv' # columns - cell id then more stuff
  gene_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/220921_transcriptome_analysis.xlsx')
  
  # project <- '2022620_857.1'
  # qc_name <- '2023-01-09'
  # additional <- '/data/vrc_his/douek_lab/wakecg/2020213_NHP857.1/preprocessing_2021-02-04_metadata.csv'
  
  sdat_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Mapped.RDS')
  #de_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Cluster_DE.RDS')
  additional <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/SampleSheets/Additional.csv')
  out_rdata <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/App/Data.RData')
  out_pdf <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Extra_plots.pdf')
  gene_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/genes.xlsx')
  gtf_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/data/gtf.RDS')
  exclude_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Excluded_genes.txt')
  
  cluster_compare <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Cluster_DE/Singles.txt')
  rds_files <- c(paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Cluster_DE.RDS'),
                 paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/SubClusters/RNA_clusters-0.RDS'),
                 paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/SubClusters/RNA_clusters-3.RDS'))

  username <- ''
  server <- 'rstudio-connect.niaid.nih.gov'
} else{
  args = commandArgs(trailingOnly=TRUE)
  out_rdata <- args[1] ### {output}
  out_pdf <- args[2]
  additional <- args[3]  ### {params}
  gene_file <- args[4]
  project <- args[5]
  username <- args[6]
  server <- args[7]
  sdat_file <- args[8] ### {input}
  de_file <- args[9]
  exclude_file <- args[10]
  gtf_file <- args[11]
  cluster_compare <- args[12]
  rds_files <- args[13:length(args)]
}

print(sdat_file)
#print(de_file)
print(gtf_file)

### Divide support RDS files into Cluster DE or sub clusters
cluster_de_files <- rds_files[grepl('/Cluster_DE', rds_files)]
sub_cluster_files <- rds_files[grepl('/SubClusters/', rds_files)]

### Read single cluster comparison file
comp_rds <- read.table(cluster_compare)[,1]
comp1 <- readRDS(comp_rds[1])

### Do direct matching from gtf
if(grepl('\\.gtf', gtf_file)){
  gtf <- read_gtf(gtf_file, feature_type = 'gene', atts_of_interest = c('gene_id', 'gene_name',  'gene_biotype'))
} else if(grepl('\\.RDS', gtf_file)){
  gtf <- readRDS(gtf_file)
} else{
  warning('No gtf file')
}
row.names(gtf) <- gtf$gene_id
### Determine whether 'gene_name' or 'gene_id' better matches those in the seurat object
gtf_cols <- c('gene_name', 'gene_id')
sdat <- readRDS(sdat_file)

if(project == '2020213_NHP857.1'){
  mis <- c('34941_wk2_3', '36186_wk2_2', '16C301_wk6_1')
  plates <- unique(sdat$plate_name)
  dat <- read.csv('/data/vrc_his/douek_lab/projects/RNASeq/2022620_857.1/Cell_sheet.csv')
  all(plates %in% dat$plate_name)
  mis %in% dat$plate_name
  table(dat$plate_name)
  table(dat[which(dat$plate_name %in% mis), 'plate_name'])
  for(plate in unique(sdat$plate_name)){
    cells <- row.names(sdat@meta.data[which(sdat$plate_name == plate),])
    dat <- sdat@assays$RNA@layers$counts[, cells]
    write.csv(dat, paste0('/data/vrc_his/douek_lab/projects/BCRSeq/2020213_NHP857.1/tempForDeposit/COUNTS-B-', plate, '.csv'), quote =F)
    gzip(paste0('/data/vrc_his/douek_lab/projects/BCRSeq/2020213_NHP857.1/tempForDeposit/COUNTS-B-', plate, '.csv'))
  } 
}

#id_type <- names(which.max(sapply(gtf_cols, function(col) sum(row.names(sdat) %in% gtf[, col]))))
id_type <- get_id_name(row.names(sdat), gtf)

### If the file exists and is not empty, read it
if(file.exists(exclude_file) & file.info(exclude_file)$size != 0){
  exclude_gene_names <- read.table(exclude_file)[,1]
} else{
  exclude_gene_names <- c()
}
features <- row.names(sdat)[which(!(row.names(sdat) %in% exclude_gene_names))]

### Save size by removing some unnecessary meta data
cols <- c('MT_sum', 'orig.ident', 'propmt', 'rna_size', 'prot_size', 'ngene', 'bc', 'Assignment_CR', 'Assignment_simple', 'cell_id', 'FCID', 'HTO_index_name', 'Cell_subset', 'Cell_Count',
          'Total_cells', 'cDNA_Concentration', 'Date_of_library_prep', 'Surface_protein_index_name', 'CSO_target_reads', 'GEx_index_name', 'VDJ_index_name', 'VDJ_target_reads', 
          'Date_Sequenced', 'well', 'plate', 'Number_of_Cells','Library_Preparation_Method', 'Index_Name', 'i5_index_ID', 'index2', 'Recipe', 'Notes', 'Sample_Project', 'Cell_Name', 'S_N', 
          'ssss_fq', 'b2f_fq', 'trim_log', 'trim_log_exists', 'trim_error', 'data_dir', 'fastq_dir', 'bam_dir', 'count_dir', 'count_files', 'bam_file', 'tab_files', 'star_files', 'trim_files')
cols <- cols[cols %in% colnames(sdat@meta.data)]
for(c in cols){
  sdat[[c]] <- NULL
}

if(additional != ''){
  adds <- read.csv(additional, sep = ',', header = T, row.names = 1)
  adds <- adds[colnames(sdat), ]
  ### 
  override <- colnames(adds)[colnames(adds) %in% colnames(sdat@meta.data)]
  if(length(override) > 0){
    print(paste0('Overriding columns ', toString(override)))
  }
  for(c in colnames(adds)){
    sdat[[c]] <- adds[,c]
  }
}

### Subcluster data
dat <- lapply(sub_cluster_files, function(file) readRDS(file))
names(dat) <- sapply(sub_cluster_files, function(file) gsub('.RDS', '', basename(file)))
### Add each dimensionality reduction and meta data column to the real cseq
for(sc in names(dat)){
  ### Add cluster info to the original object
  sdat@meta.data[, colnames(dat[[sc]]@meta.data)] <- dat[[sc]]@meta.data
  ### Add UMAP info to the original object
  ### Having these subcluster UMAPs breaks ability to subset sdat. Removing them fixes the problem.
  for(umap in names(dat[[sc]]@reductions)){
    sdat@reductions[paste0(umap, '_umap')] <- dat[[sc]]@reductions[umap]
    #### To access these slots, we provide the Embeddings, Loadings, and Stdev functions
    # ### Subcluster Embedding 
    # dat <- Embeddings(dat[[sc]], reduction = umap)
    # abs <- colnames(sdat)[which(!colnames(sdat) %in% row.names(dat))]
    # newdat <- as.data.frame(matrix(nrow = length(abs), ncol = length(colnames(dat))))
    # row.names(newdat) <- abs
    # colnames(newdat) <- colnames(dat)
    # dat <- rbind(dat, newdat)
    # ### Reorder to match sdat
    # dat <- dat[colnames(sdat), ]
    # ### Modified subcluser Emedding in the original seurat object
    # Embeddings(sdat, reduction = paste0(umap, '_umap')) <- dat
  }
  sdat_plot <- subset(x = sdat, subset = Celltype == "CD4")
  sdat@reductions[paste0(umap, '_umap')] <- NULL
}
rm(dat)

print(names(sdat@reductions))
### removing some assays and reductions
DefaultAssay(sdat) <- 'RNA'
umaps <- names(sdat@reductions)[grepl('umap', names(sdat@reductions))]
print(umaps)
### Rename 'umap' as 'RNA_umap'
# sdat@reductions$RNA_umap <- sdat@reductions$umap
# umaps[which(umaps == 'umap')] <- 'RNA_umap'

pdf(out_pdf)
##### Some basic plots
### Specific gene violin plots
if(gene_file != '' & file.exists(gene_file) & file.info(gene_file)$size != 0){
  genes <- as.data.frame(read_excel(gene_file, sheet = 1))
  colnames(genes) <- c('gene_name', 'category')
  ### In the case that a gene has been added to multiple categories
  dups <- names(which(table(genes$gene_name) > 1))
  print(paste0(toString(dups), ' appear in more than one category'))
  
  #id_type2 <- names(which.max(sapply(gtf_cols, function(col) sum(genes %in% gtf[, col]))))
  # id_type2 <- get_id_name(genes$gene_name, gtf)
  # if(id_type != id_type2){
  #   res <- gene_matching_v1(genes$gene_name, row.names(sdat), gtf, one2multi = 'exclude')
  #   ### Remove those without a match
  #   res <- res[which(!is.na(res[, id_type])),]
  #   genes <- genes[which(genes[, id_type2] %in% res[, id_type2]),]
  #   genes[, id_type] <- sapply(genes$gene_name, function(x) res[which(res$gene_name == x), id_type])
  # }

  ### For each category, plot genes together (up to 4 each)
  for(cat in unique(genes$category)){
    print(cat)
    genes_sub <- genes[which(genes$category == cat), ]
    #id_type2 <- get_id_name(genes_sub$gene_name, gtf)
    id_type2 <- 'gene_name'
    if(id_type != id_type2){
      res <- gene_matching_v1(genes_sub$gene_name, row.names(sdat), gtf, one2multi = 'exclude')
      ### In case of no matches at all
      if(!id_type2 %in% colnames(res)){
        res <- NULL
      }
      ### Remove those without a match
      res <- res[which(!is.na(res[, id_type])),]
      ### If no genes remaining, just skip to next category
      if(length(res) == 0){
        next
      }
      ### 
      genes_sub <- genes_sub[which(genes_sub[, id_type2] %in% res[, id_type2]),]
      genes_sub[, id_type] <- sapply(genes_sub$gene_name, function(x) res[which(res$gene_name == x), id_type])
    }
    
    ### Gene names of this category
    #gene_names <- genes_sub[which(genes_sub$category == cat), 'gene_name'][[1]]
    gene_group <- genes_sub[which(genes_sub$category == cat), ]
    row.names(gene_group) <- 1:dim(gene_group)[1]
    ### Split into groups if there are more than 4.
    if(length(row.names(gene_group)) > 4){
      #groups <- split(gene_names, cut(1:length(gene_names), length(gene_names) %/% 4 + 1))
      groups <- split(1:dim(gene_group)[1], cut(1:dim(gene_group)[1], dim(gene_group)[1] %/% 4 + 1))
      groups <- lapply(groups, function(x) gene_group[x,])
    } else{
      groups <- list(gene_group)
 
    }
    
    title <- ggdraw() + draw_label(cat, fontface = 'bold', x = 0, hjust = 0) + theme(plot.margin = margin(0, 0, 0, 250))
    ### For each group (may be only 1) of up to four genes, plot violin plots together
    for(gdat in groups){
      ### For each gene in the group
      vplots <- lapply(1:dim(gdat)[1], function(i) my_vln_plot(sdat, gdat[i, 'gene_id'], title = gdat[i, 'gene_name'], de = NA, groupby = 'seurat_clusters'))
      pmain <- plot_grid(plotlist = vplots)
      print(plot_grid(title, pmain, ncol = 1 , rel_heights = c(0.1, 1)))
    }
  }
}

# sdat@meta.data[which(is.na(sdat$TCR.status)), 'TCR.status'] <- 'Unlabeled'
# sdat@meta.data[which(is.na(sdat$TCR.status)), 'TCR.status_classification'] <- 'Unclassified'
Idents(sdat) <- 'seurat_clusters'

classification_options <- c('predicted_label_fine', 'seurat_clusters')
classification_options <- classification_options[which(classification_options %in% colnames(sdat@meta.data))]
umap_color <- classification_options[1]
#umap_color <- 'predicted_label_main'
pal <- hue_pal()(length(unique(sdat@meta.data[, umap_color])))
cells.highlight <- c()
cols.highlight <- pal[as.character(sdat@meta.data[, umap_color])]

#DimPlot(sdat, group.by = umap_color, label = T, reduction = umaps[1], cols = pal)
DimPlot(sdat, group.by = umap_color, label = F, reduction = umaps[1], cols = pal)

dev.off()


### App data
sdat <- DietSeurat(sdat, layers = c('counts'), dimreducs = umaps, assays = c('RNA'))
rna_by_rna <- readRDS(cluster_de_files[1])

save.image(out_rdata)

total_bytes <- as.numeric(object.size(sdat)[1]) + as.numeric(object.size(rna_by_rna)[1])
deploy_app <- T
if(total_bytes >= 3145728000){
  deploy_app <- F
  print(paste0('App will not be deployed because it is too large (', total_bytes, ' bytes)'))
}

if(project != '' & username != '' & server != '' & deploy_app == T){
  project <- gsub('\\.','_', project)
  app_dir <- paste0(dirname(sdat_file), '/App/')
  deployApp(appDir = app_dir, appName = project, account = username, server = server)
}

