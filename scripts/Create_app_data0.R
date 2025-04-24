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
  qc_name <- 'Run2022-11-14'
  additional <- '/data/vrc_his/douek_lab/projects/RNASeq/2021600_kristin/SampleSheets/Additional.csv' # columns - cell id then more stuff
  gene_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/220921_transcriptome_analysis.xlsx')
  
  # project <- '2022620_857.1'
  # qc_name <- '2022-11-01'
  #additional <- '/data/vrc_his/douek_lab/wakecg/2020213_NHP857.1/preprocessing_2021-02-04_metadata.csv' 
  
  sdat_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Mapped.RDS')
  de_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Cluster_DE.RDS')
  additional <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/SampleSheets/Additional.csv')
  out_rdata <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/App/Data.RData')
  out_pdf <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Extra_plots.pdf')
  gene_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/genes.xlsx')
  gtf_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/data/gtf.RDS')
  exclude_file <- paste0('/data/vrc_his/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Excluded_genes.txt')
  
  username <- 'wakecg'
  server <- 'rstudio-connect.niaid.nih.gov'
} else{
  args = commandArgs(trailingOnly=TRUE)
  sdat_file <- args[1]
  de_file <- args[2]
  exclude_file <- args[3]
  gtf_file <- args[4]
  additional <- args[5]
  gene_file <- args[6]
  out_rdata <- args[7]
  out_pdf <- args[8]
  project <- args[9]
  username <- args[10]
  server <- args[11]
}

print(sdat_file)
print(de_file)
print(gtf_file)

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

print(names(sdat@reductions))
### removing some assays and reductions
DefaultAssay(sdat) <- 'RNA'
umaps <- names(sdat@reductions)[grepl('umap', names(sdat@reductions))]
print(umaps)

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
    genes_sub <- genes[which(genes$category == cat), ]
    id_type2 <- get_id_name(genes_sub$gene_name, gtf)
    if(id_type != id_type2){
      res <- gene_matching_v1(genes_sub$gene_name, row.names(sdat), gtf, one2multi = 'exclude')
      ### Remove those without a match
      res <- res[which(!is.na(res[, id_type])),]
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

umap_color1 <- 'predicted_label_fine'
#umap_color1 <- 'predicted_label_main'
pal <- hue_pal()(length(unique(sdat@meta.data[, umap_color1])))
cells.highlight <- c()
cols.highlight <- pal[as.character(sdat@meta.data[, umap_color1])]

#DimPlot(sdat, group.by = umap_color1, label = T, reduction = umaps[1], cols = pal)
DimPlot(sdat, group.by = umap_color1, label = F, reduction = umaps[1], cols = pal)

dev.off()


### App data
sdat <- DietSeurat(sdat, counts = F, dimreducs = umaps, assays = c('RNA'))
rna_by_rna <- readRDS(de_file)
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

