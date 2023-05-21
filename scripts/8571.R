R.Version()
print(.libPaths())

library('sys')
library('Seurat')
library('stringr')
library('pheatmap')
library('ggplot2')
library('umap')
library('textshape')
library('plyr')
library('dplyr')
library('biomaRt')
library('grid')
library('scales')
library('rsconnect')
library('readxl')
library('cowplot')
library('RColorBrewer')
library('data.table')
library('viridis')
library('WriteXLS')
library('VennDiagram')
library('fgsea')
library('circlize')
library('ComplexHeatmap')
library('readxl')
library('svglite')

source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R')

#qc_name <- '2022-11-01'
#qc_name <- '2022-11-16'
#qc_name <- '2022-12-09'
qc_name <- '2023-01-09'
#qc_name <- '2023_res0.3'

font <- 'Arial'
fig_theme = theme_bw() + 
  theme(text = element_text(family = font, color = "black"),
        axis.title = element_text(family = font, size = 8, color = "black"),
        axis.text  = element_text(family = font, size = 8, color = "black"),
        legend.title = element_text(family = font, size = 7, color = "black"),
        legend.text  = element_text(family = font, size = 7, color = "black"),
        legend.spacing = unit(5,'points'), 
        legend.key.size = unit(8, 'points'),
        legend.box.spacing = unit(0, 'points'),
        axis.ticks.length = unit(1, 'points'))

load(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/App/Data.RData'))
additional1 <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/SampleSheets/merged_data_for_chris.tsv'
additional2 <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/SampleSheets/Additional.csv'
additional3 <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/SampleSheets/VRC857.1_B_cell_sort_MFIs.xlsx'

strats <- c('All', 'seurat_clusters-0', 'seurat_clusters-1', 'seurat_clusters-2', 'seurat_clusters-3')
strats <- c('All', 'seurat_clusters-0', 'seurat_clusters-1+3', 'seurat_clusters-2')
names(strats) <- c('All', 'LZ-like', 'AM', 'RM')

de_files <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/DE/timepoint/', strats, '/DE_results.tsv')
#gmt_file <- '/hpcdata/vrc/vrc1_data/douek_lab/wakecg/genesets/h.all.v2022.1.Hs.symbols.gmt'
gmt_file <- '/hpcdata/vrc/vrc1_data/douek_lab/wakecg/genesets/c2.cp.v2022.1.Hs.symbols.gmt'

geneset_set <- ''
#geneset_set <- '_C2CP'
fgsea_files <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/DE/timepoint/', strats, '/fgsea', geneset_set, '.tsv')
gene_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/genes.xlsx'
gtf_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/data/gtf.RDS'
exclude_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/Excluded_genes.txt')

b_genes <- c('BCL2A1', 'FABP5', 'CD40LG', 'IL21', 'TMSB4X', 'FTH1', 'PTPRCAP', 'DNAJB1', 'CCL4L1', 'CCL3', 'IL10', 'IFNG')

d_sets <- c('Cell signaling','Cell signaling','Cell signaling','Cell signaling','Cellular component', 'Cytokine signaling', 'Development','Development','Development',
            'Differentiation', 'Differentiation', 'DNA damage', 'Immune response', 'Immune response', 'Inflammatory response', 'Inflammatory response', 'Interferon response', 
            'Interferon response', 'Metabolism','Metabolism','Metabolism','Proliferation','Proliferation','Proliferation' )
names(d_sets) <- c('Kras Signaling Up', 'Kras Signaling Dn', 'Hedgehog Signaling', 'Estrogen Response Early', 'Apical Junction', 'IL6 Jak Stat3 Signaling', 'Pancreas Beta Cells',
            'Epithelial Mesenchymal Transition', 'Angiogenesis', 'Myogenesis', 'Myc Targets V1', 'UV Response Up', 'Coagulation', 'Allograft Rejection', 'TNFa Signaling via NFkB', 
            'Inflammatory Response', 'Interferon Alpha Response', 'Oxidative Phosphorylation', 'Hypoxia', 'Heme Metabolism', 'Mitotic Spindle', 'E2F Targets', 'Apoptosis')


adds <- read.csv(additional1, sep = '\t', header = T)
adds$cell_id <- sapply(adds$cell_id, function(x) gsub('-', '_', x))
adds$cell_name <- sapply(adds$cell_id, function(x) gsub('\\.', '_', x))
row.names(adds) <- adds$cell_name

cols <- c('clone_count', 'VH_fam', 'KL', 'isotype', 'VH_SHM_PCT', 'CDRH3_AA_LEN', 'expanded_clone')
adds <- adds[which(adds$cell_name %in% sdat$cell_name), ]
adds <- adds[, c('cell_name', cols)]
adds <- merge(sdat@meta.data[, c('Cell_ID', 'cell_name')], adds, by = 'cell_name', all = T)
row.names(adds) <- adds$Cell_ID
adds <- adds[row.names(sdat@meta.data), ]

for(c in cols){
  sdat[[c]] <- adds[,c]
}

adds <- read.csv(additional2, sep = ',', header = T, row.names = 1)
adds <- adds[colnames(sdat), ]
### 
override <- colnames(adds)[colnames(adds) %in% colnames(sdat@meta.data)]
if(length(override) > 0){
  print(paste0('Overriding columns ', toString(override)))
}
for(c in colnames(adds)){
  sdat[[c]] <- adds[,c]
}

sdat$epitope <- gate_collapse(sdat)
cols <- c('epitope', cols)

### MFI
adds <- read_excel(additional3)
adds <- as.data.frame(adds[, c('Sample ID', 'Plate', 'Well', 'Median (*G780 CD19 CY7PE-A)', 'Median (*U785 CD20 BUV 805-A)', 'Median (*V605 CD27 BV605-A)')])
colnames(adds) <- c('Subject_ID', 'Plate', 'Well', 'CD19_MFI', 'CD20_MFI', 'CD27_MFI')
adds$Plate <- sapply(adds$Plate, function(x) gsub('P0', '', x))
adds$Plate <- sapply(adds$Plate, function(x) gsub('a', '', x))
adds$Plate <- sapply(adds$Plate, function(x) gsub('b', '', x))
adds$Subject_ID <- sapply(adds$Subject_ID, function(x) gsub(' WK', '_wk', x))
adds$cell_id <- paste(adds$Subject_ID, adds$Plate, adds$Well, sep = '_')
row.names(adds) <- adds$cell_id

all(sdat$cell_name %in% row.names(adds))
adds <- adds[sdat$cell_name, ]
sdat[['CD19_MFI']] <- adds[,'CD19_MFI']
sdat[['CD20_MFI']] <- adds[,'CD20_MFI']
sdat[['CD27_MFI']] <- adds[,'CD27_MFI']
cols <- c('CD19_MFI', 'CD20_MFI','CD27_MFI', cols)

### If the file exists and is not empty, read it
if(file.exists(exclude_file) & file.info(exclude_file)$size != 0){
  exclude_gene_names <- read.table(exclude_file)[,1] 
} else{
  exclude_gene_names <- c()
}

### Do direct matching from gtf
if(grepl('\\.gtf', gtf_file)){
  gtf <- read_gtf(gtf_file, feature_type = 'gene', atts_of_interest = c('gene_id', 'gene_name',  'gene_biotype'))
} else if(grepl('\\.RDS', gtf_file)){
  gtf <- readRDS(gtf_file)
}
row.names(gtf) <- gtf$gene_id
### Determine whether 'gene_name' or 'gene_id' better matches those in the seurat object
gtf_cols <- c('gene_name', 'gene_id')
#id_type <- names(which.max(sapply(gtf_cols, function(col) sum(row.names(sdat) %in% gtf[, col]))))
id_type <- get_id_name(row.names(sdat), gtf)
all(b_genes %in% gtf$gene_name)
gtf$gene_plot_name <- sapply(row.names(gtf), function(x) ifelse(is.na(gtf[x, 'gene_name']), x, gtf[x, 'gene_name']))

### Read set of genes
DB <- gmtPathways(gmt_file)
names(DB) <- sapply(names(DB), function(x) gsub('HALLMARK_', '', x))


subject_key <- c('#1B9E77', '#D95F02', '#7570B3', '#E7298A', '#66A61E', '#E6AB02', '#A6761D', '#666666')
names(subject_key) <- c('16C222', '16C235', '16C237', '16C283', '16C301', '16C303', '34941', '36186')

time_key <- c('#A07178', '#AFC2D5', '#676F54')
names(time_key) <- c('baseline', 'wk2', 'wk6')
time_key <- time_key[c('wk2', 'wk6')]

epitope_key <- c('#1B9E77', '#D95F02', '#E7298A', '#66A61E', '#E6AB02')
epitope_key <- c('#D71D1D', '#FEAE61', '#ABDDA4', '#2A83BB', '#888888')
names(epitope_key) <- levels(sdat$epitope)

#cluster_key <- c('#A07178', '#AFC2D5', '#676F54')
#names(cluster_key) <- c('Unknown', 'Resting', 'Activated')
#cluster_key <- c('#AFC2D5', '#A07178')
#names(cluster_key) <- c('Resting', 'Activated')
#cluster_key <- c('#AFC2D5', '#1B9E77', '#A07178', '#D95F02')
#names(cluster_key) <- c('cluster 0', 'cluster 2', 'cluster 1', 'cluster 3')
#names(cluster_key) <- c('cluster 0', 'AM1', 'AM2 cluster1', 'AM2 cluster2')

cluster_key <- c('#88CCEE', '#CC6677', '#DDCC77', '#888888')
names(cluster_key) <- c('LZ-like', 'RM', 'AM cluster1', 'AM cluster2')

### Merged 1 and 3
cluster_key_merged <- c('#88CCEE', '#CC6677', '#DDCC77')
names(cluster_key_merged) <- c('LZ-like', 'RM', 'AM')

### Name clusters, AM separated
cluster_names <- names(cluster_key)
names(cluster_names) <- c('seurat_clusters-0', 'seurat_clusters-2', 'seurat_clusters-1', 'seurat_clusters-3')
sdat$cluster_names <- cluster_names[paste0('seurat_clusters-', sdat$seurat_clusters)]
sdat$cluster_names <- factor(sdat$cluster_names, levels = names(cluster_key))

### Name clusters, AM merged
cluster_names_merged <- c('LZ-like', 'RM', 'AM', 'AM')
names(cluster_names_merged) <- c('seurat_clusters-0', 'seurat_clusters-2', 'seurat_clusters-1', 'seurat_clusters-3')
sdat$cluster_names_merged <- cluster_names_merged[paste0('seurat_clusters-', sdat$seurat_clusters)]
sdat$cluster_names_merged <- factor(sdat$cluster_names_merged, levels = unique(cluster_names_merged))


### features from the analysis
features <- row.names(sdat)[which(!(row.names(sdat) %in% exclude_gene_names))]

Idents(sdat) <- sdat$cluster_names
# de <- FindMarkers(sdat, ident.1 = 'Resting', ident.2 = 'Activated', min.pct = 0, logfc.threshold = 0, min.cells.group = 0, test.use = 'LR', features = features)
# de[, 'gene_name'] <- gtf[row.names(de), 'gene_name']
# de$gene_id <- row.names(de)
de <- FindMarkers(sdat, ident.1 = 'AM2 cluster1', ident.2 = 'AM2 cluster2', min.pct = 0, logfc.threshold = 0, min.cells.group = 0, test.use = 'LR', features = features)
de[, 'gene_name'] <- gtf[row.names(de), 'gene_name']
de$gene_id <- row.names(de)

### Excel summary for cluster DE
rna_by_rna[, 'gene_name'] <- gtf[rna_by_rna$gene, 'gene_name']
#rna_by_rna <- rna_by_rna[, c('gene', 'gene_name', 'p_val', 'p_val_adj', 'avg_log2FC', 'diff', 'pct.1', 'pct.2', 'cluster')]
clusters <- unique(rna_by_rna$cluster)
des <- lapply(clusters, function(x) rna_by_rna[which(rna_by_rna$cluster == x),  c('gene', 'gene_name', 'p_val', 'p_val_adj', 'avg_log2FC', 'diff', 'pct.1', 'pct.2')])
names(des) <- paste0('Cluster ', clusters)

#des[['Resting-Activated']] <- de
des[['AM2_comparison']] <- de

WriteXLS(des, ExcelFileName = paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/Cluster_DE.xls'))

write.csv(sdat@meta.data[, c('cell_name', 'cluster_names')], file = paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/cluster_membership.csv'))
### Remove cluster 2
#sdat <- subset(sdat, cells = row.names(sdat@meta.data[which(sdat$seurat_clusters != 2),]))
#sdat$cluster_names <- factor(sdat$cluster_names, levels = c('Resting', 'Activated'))

### Blended FeaturePlot
genes <- c('ENSMMUG00000016579', 'ENSMMUG00000013289')
pdf_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/FeaturePlot.pdf')
pdf(pdf_file)
DimPlot(sdat)
ps <- FeaturePlot(sdat, features = genes, blend = T, combine = T)
p1 <- ps[[3]] + theme(legend.position = "none") + ggtitle(toString(gtf[genes, 'gene_plot_name']))
p2 <- ps[[4]] 
p1+p2
dev.off()

p3 <- p1+p2
png_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/FeaturePlot.png')
png(png_file, width = 700)
p3
dev.off()

##### Some basic plots
a <- rowSums(sdat@assays$RNA@data > 0) / length(colnames(sdat))
a <- a[features]
hist(a, xlab = 'Fraction of total cells', main = 'Features')
VlnPlot(sdat, features = 'nFeature_RNA', group.by = 'seurat_clusters')
sdat$N_reads <- as.numeric(sdat$N_reads)
VlnPlot(sdat, features = 'N_reads', group.by = 'seurat_clusters')
sdat[['CD19_MFI']] <- log(sdat[['CD19_MFI']] + 1)
sdat[['CD20_MFI']] <- log(sdat[['CD20_MFI']] + 1)
sdat[['CD27_MFI']] <- log(sdat[['CD27_MFI']] + 1)

pdf(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/CD27_MFI.pdf'))
VlnPlot(sdat, features = 'CD19_MFI', group.by = 'seurat_clusters')
FeaturePlot(sdat, 'CD19_MFI')
VlnPlot(sdat, features = 'CD20_MFI', group.by = 'seurat_clusters')
FeaturePlot(sdat, 'CD20_MFI')
VlnPlot(sdat, features = 'CD27_MFI', group.by = 'seurat_clusters')
FeaturePlot(sdat, 'CD27_MFI')
dev.off()

##### Expression plots, violin and bubble
### FindAllMarkers
pdf_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/FAM.pdf')
cats <- c('1')
pthresh <- 0.01
genes <- rna_by_rna[which(rna_by_rna$p_val_adj < pthresh & rna_by_rna$cluster == cats[1]), c('gene', 'cluster')]
topn <- 20
while(dim(genes)[1] > topn){
  pthresh <- pthresh/10
  genes <- rna_by_rna[which(rna_by_rna$p_val_adj < pthresh & rna_by_rna$cluster == cats[1]), c('gene', 'cluster')]
}
if(dim(genes)[1] == 0){
  pthresh <- pthresh * 10
  genes <- rna_by_rna[which(rna_by_rna$p_val_adj < pthresh & rna_by_rna$cluster == cats[1]), c('gene', 'cluster')]
}
colnames(genes) <- c('gene_id', 'category')
cats <- c('1')
expression_plots(pdf_file, sdat, genes, gtf, features, cats)

### Heatmap to evaluate clustering resolution by eye
sdat <- ScaleData(sdat, features = row.names(sdat))
# rna_by_rna %>%
#   filter(p_val_adj < 0.05) %>%
#   group_by(cluster) %>%
#   top_n(n = 5, wt = avg_log2FC) -> top10

topn <- 50
top10 <- unique(unlist(lapply(unique(rna_by_rna$cluster), function(c) {
  csub <- rna_by_rna[which(rna_by_rna$cluster == c & rna_by_rna$p_val_adj < 0.05),]
  csub <- csub[order(abs(csub$avg_log2FC), decreasing = T),]
  return(csub[1:topn, 'gene'])
})))

sum(top10 %in% row.names(rna_by_rna[which(rna_by_rna$cluster == 1),]))

pdf(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/FAM_heatmap.pdf'))
DoHeatmap(sdat, features = top10) #+ NoLegend()
dev.off()
### Resting-Acitvated DE
pdf_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/Resting-Activated_DE.pdf')
pthresh <- 0.01
genes <- de[which(de$p_val_adj < pthresh), c('gene_id'), drop = F]
topn <- 20
while(dim(genes)[1] > topn){
  pthresh <- pthresh/10
  genes <- de[which(de$p_val_adj < pthresh), c('gene_id'), drop = F]
}
if(dim(genes)[1] == 0){
  pthresh <- pthresh * 10
  genes <- de[which(de$p_val_adj < pthresh), c('gene_id'), drop = F]
}
genes$category <- 'Resting-Activated'
cats <- c('Resting-Activated')
expression_plots(pdf_file, sdat, genes, gtf, features, cats, paste0('Resting-Activated DE, adjp < ', pthresh))


### Manual additions
fig_genes <- c('CD19', 'CD79B', 'CCL3', 'BCL2A1', 'MYC', 'TOX2', 'FGR', 'ITGAX', 'CR2', 'LTB', 'TCF7')
### 
if(gene_file != '' & file.exists(gene_file) & file.info(gene_file)$size != 0){
  dotplot_genes <- as.data.frame(read_excel(gene_file, sheet = 1))
  colnames(dotplot_genes) <- c('gene_name', 'gene_id', 'category')
  ### Adds manual additions
  dotplot_genes <- rbind(dotplot_genes, data.frame('gene_name' = fig_genes[which(!fig_genes %in% dotplot_genes$gene_name)], 'gene_id' = NA, 'category' = 'fig3b'))
  dotplot_genes <- complete_gene_table(dotplot_genes, gtf)
  
  ### Tweak names
  #dotplot_genes[which(dotplot_genes$category == 'B'), 'category'] <- 'B cell'
  #id_type2 <- names(which.max(sapply(gtf_cols, function(col) sum(dotplot_genes %in% gtf[, col]))))
  pdf_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/B_cells.pdf')
  cats <- c('B cell', 'B activation', 'Activated B subtype')
  #dotplot_genes[which(dotplot_genes$category %in% cats), 'category'] <- paste0(dotplot_genes[which(dotplot_genes$category %in% cats), 'category'], '\nfeatures')
  expression_plots2(sdat, dotplot_genes, gtf, features, cats, 'cluster_names', pdf_file = pdf_file)
  
  pdf_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/T_cells.pdf')
  expression_plots2(sdat, dotplot_genes, gtf, features, c('T cell', 'T an B'), 'cluster_names', pdf_file = pdf_file)
  
  pdf_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/Chaim2.pdf')
  expression_plots2(sdat, dotplot_genes, gtf, features, c('XYZ'), 'cluster_names', pdf_file = pdf_file)
  expression_plots2(sdat, dotplot_genes, gtf, features, c('XYZ'), 'timepoint', pdf_file = pdf_file)
  
  pdf_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/Fig3B.pdf')
  #fig_genes <- c('CD19', 'CD79B', 'CCL3', 'FOSL1', 'FOSL2', 'MYC', 'TOX2', 'DTX1', 'ITGAX', 'CR2', 'VPREB3', 'TCF7')
  fig_genes <- c('CD19', 'CD79B', 'CCL3', 'BCL2A1', 'MYC', 'TOX2', 'FGR', 'ITGAX', 'CR2', 'LTB', 'TCF7')
  dotplot_genes <- dotplot_genes[which(dotplot_genes$gene_name %in% fig_genes), ] 
  dotplot_genes$category <- 'fig3b'
  row.names(dotplot_genes) <- dotplot_genes$gene_name
  dotplot_genes <- dotplot_genes[fig_genes, ]
  #dotplot_genes[which(dotplot_genes$gene_name == 'CR2'), 'gene_plot_name'] <- 'CD21'
  sdat$cluster_names <- factor(sdat$cluster_names, levels = names(cluster_key))
  
  expression_plots2(sdat, dotplot_genes, gtf, features, c('fig3b'), 'cluster_names', cat_label = F, scale = F, scale_by = 'size', pdf_file = pdf_file)
  fig3b <- expression_plots2(sdat, dotplot_genes, gtf, features, c('fig3b'), 'cluster_names', cat_label = F, scale = T, scale_by = 'size')
  fig3b <- fig3b + 
    fig_theme + 
    theme(axis.title = element_blank(),
           axis.ticks.length = unit(4, 'points'),
           axis.text.x = element_text(angle = 45, hjust = 0.5, margin = margin(t = 10, r = 0, b = 0, l = 0, unit = "point")),
           panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(),
           panel.background = element_blank(), 
           panel.border = element_blank(),
           axis.line = element_line(colour = "black"),
           legend.box.margin = margin(t = 3, r = 10, b = 1, l = 1, unit = "point"))
  fig3b
  ggsave(path = paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/svg/'), filename = 'Fig3B.svg', 
         plot = fig3b, device = 'svg', units = 'in', width = 4.0, height = 2.5)
  
}

### Frequency of clusters by timepoint 
pdf(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/Fig3C.pdf'))
plot_datC <- sdat@meta.data[, c('timepoint', 'cluster_names_merged')]
plot_datC <- as.data.frame(table(paste(sdat$timepoint, sdat$cluster_names_merged, sep = '_')))
plot_datC[, c('timepoint', 'cluster_names_merged')] <- as.data.frame(t(as.data.frame(lapply(plot_datC$Var1, function(x) strsplit(as.character(x), '_')[[1]]))))
plot_datC[which(plot_datC$timepoint == 'wk2'), 'Frac'] <- plot_datC[which(plot_datC$timepoint == 'wk2'), 'Freq']/sum(plot_datC[which(plot_datC$timepoint == 'wk2'), 'Freq'])
plot_datC[which(plot_datC$timepoint == 'wk6'), 'Frac'] <- plot_datC[which(plot_datC$timepoint == 'wk6'), 'Freq']/sum(plot_datC[which(plot_datC$timepoint == 'wk6'), 'Freq'])
plot_datC$Perc <- paste0(signif(plot_datC$Frac*100, digits = 3), ' %')
plot_datC$Frac <- signif(plot_datC$Frac, digits = 3)

plot_datC$cluster_names_merged <- factor(plot_datC$cluster_names_merged, levels = unique(cluster_names_merged))

clusts <- levels(plot_datC$cluster_names_merged)
clusts <- clusts[length(clusts):1]
for(i in 1:length(clusts)){
  clust <- clusts[i]
  ### wk2
  height0 <- plot_datC[which(plot_datC$timepoint == 'wk2' & plot_datC$cluster_names_merged == clust), 'Frac'] / 2
  if(i == 1){
    bottom_height <- 0
  } else{
    bottom_height <- sum(plot_datC[which(plot_datC$timepoint == 'wk2' & plot_datC$cluster_names_merged %in% clusts[1:(i-1)]), 'Frac'])
  }
  plot_datC[which(plot_datC$timepoint == 'wk2' & plot_datC$cluster_names_merged == clust), 'text_height'] <- height0 + bottom_height
  ### wk6
  height0 <- plot_datC[which(plot_datC$timepoint == 'wk6' & plot_datC$cluster_names_merged == clust), 'Frac'] / 2
  if(i == 1){
    bottom_height <- 0
  } else{
    bottom_height <- sum(plot_datC[which(plot_datC$timepoint == 'wk6' & plot_datC$cluster_names_merged %in% clusts[1:(i-1)]), 'Frac'])
  }
  plot_datC[which(plot_datC$timepoint == 'wk6' & plot_datC$cluster_names_merged == clust), 'text_height'] <- height0 + bottom_height
}
plot_datC$timepoint <- gsub('wk', 'week ', plot_datC$timepoint)
fig3c <- ggplot(data = plot_datC, aes(fill = cluster_names_merged, x = timepoint, y = Freq)) +
  fig_theme + 
  geom_bar(position = 'fill', stat="identity") +
  geom_text(aes(x = timepoint, y = text_height, label = Perc), size = 3, family = 'Arial') + 
  scale_fill_manual(values = cluster_key_merged, name = 'Cluster') +
  ylab('Frequency') +
  theme(panel.border = element_blank(),
        axis.title.x=element_blank()) 

fig3c 
ggsave(path = paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/svg/'), filename = 'Fig3C.svg', 
       plot = fig3c, device = 'svg', units = 'in', width = 3.1, height = 3.1)


dev.off()


### Frequency of clusters by monkey 
plot_dat <- sdat@meta.data[, c('subject', 'cluster_names')]
plot_dat <- as.data.frame(table(paste(sdat$subject, sdat$cluster_names, sep = '_')))
plot_dat[, c('subject', 'cluster_names')] <- as.data.frame(t(as.data.frame(lapply(plot_dat$Var1, function(x) strsplit(as.character(x), '_')[[1]]))))
p <- ggplot(data = plot_dat, aes(fill = cluster_names, x = subject, y = Freq)) +
  geom_bar(position = 'fill', stat="identity")
p <- p + scale_fill_manual(values = cluster_key)
p 
### Frequency of epitopes by timepoint
plot_dat <- sdat@meta.data[, c('timepoint', 'epitope')]
plot_dat <- as.data.frame(table(paste(sdat$timepoint, sdat$epitope, sep = '_')))
plot_dat[, c('timepoint', 'epitope')] <- as.data.frame(t(as.data.frame(lapply(plot_dat$Var1, function(x) strsplit(as.character(x), '_')[[1]]))))
### Remove NAs
plot_dat <- plot_dat[!(is.na(plot_dat[, 'epitope']) | plot_dat[, 'epitope'] == 'NA'), ]
p <- ggplot(data = plot_dat, aes(fill = epitope, x = timepoint, y = Freq)) +
  geom_bar(position = 'fill', stat="identity")
p <- p + scale_fill_manual(values = epitope_key)
p

cols_num <- cols[sapply(cols, function(c) {length(unique(sdat@meta.data[which(!is.na(sdat@meta.data[, c])), c])) >= 8})]
cols_factor <- cols[sapply(cols, function(c) {length(unique(sdat@meta.data[which(!is.na(sdat@meta.data[, c])), c])) <= 8})]
fig3_cols <- c('clone_count', 'VH_SHM_PCT')

pdf(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/Cluster_associations.pdf'))
sdat$isotype[which(sdat$isotype == 'No result')] <- NA
for(c in cols_factor){
  ### Frequency of epitopes by cluster_names
  plot_dat <- as.data.frame(table(paste(sdat$cluster_names, sdat@meta.data[, c], sep = '_')))
  plot_dat[, c('cluster_names', c)] <- as.data.frame(t(as.data.frame(lapply(plot_dat$Var1, function(x) strsplit(as.character(x), '_')[[1]]))))
  ### Remove NAs
  plot_dat <- plot_dat[!(is.na(plot_dat[, c]) | plot_dat[, c] == 'NA'), ]
  plot_dat$cluster_names <- factor(plot_dat$cluster_names, levels = cluster_names)
  
  n <- length(unique(plot_dat[, c]))
  
  color_key <- subject_key[1:n]
  names(color_key) <- unique(plot_dat[, c])
  
  colnames(plot_dat) <- c('Var1', 'Freq', 'cluster_names', 'Var2')
  p <- ggplot(data = plot_dat, aes(fill = Var2, x = cluster_names, y = Freq)) +
    geom_bar(position = 'fill', stat = "identity")
  p <- p + scale_fill_manual(values = color_key, name = c)
  print(p)
}

sdat$N_mapped <- as.numeric(sdat$N_mapped)
for(c in cols_num){
  plot_dat <- sdat@meta.data[, c('cluster_names', c)]
  plot_dat <- plot_dat[which(!is.na(sdat@meta.data[, c])), ]
  colnames(plot_dat) <- c('cluster', 'Var1')
  
  mu <- ddply(plot_dat, "cluster", summarise, grp.mean=mean(Var1))
  p <- ggplot(plot_dat, aes(x = Var1, color = cluster)) + 
    geom_density() + scale_color_manual(values = cluster_key) + labs(x = c) + 
    geom_vline(data = mu, aes(xintercept = grp.mean, color = cluster_key), linetype = "dashed")
  print(p)
}
dev.off()

### Epitope
plot_dat <- as.data.frame(table(paste(sdat$cluster_names_merged, sdat@meta.data[, 'epitope'], sep = '_')))
plot_dat[, c('cluster_names_merged', 'epitope')] <- as.data.frame(t(as.data.frame(lapply(plot_dat$Var1, function(x) strsplit(as.character(x), '_')[[1]]))))
### Remove NAs
plot_dat <- plot_dat[!(is.na(plot_dat[, 'epitope']) | plot_dat[, 'epitope'] == 'NA'), ]
plot_dat$cluster_names_merged <- factor(plot_dat$cluster_names_merged, levels = names(cluster_key_merged))

n <- length(unique(plot_dat[, 'epitope']))
color_key <- subject_key[1:n]
names(color_key) <- unique(plot_dat[, 'epitope'])

colnames(plot_dat) <- c('Var1', 'Freq', 'cluster_names_merged', 'Var2')
fig3d <- ggplot(data = plot_dat, aes(fill = Var2, x = cluster_names_merged, y = Freq)) + xlab('Cluster') + 
  fig_theme +
  geom_bar(position = 'fill', stat = "identity") + 
  scale_fill_manual(values = epitope_key, name = 'Epitope') +
  ylab('Frequency') +
  theme(panel.border = element_blank(),
        axis.title.x=element_blank())
  #theme(legend.justification = c(100,0.7))
  #theme(legend.justification = c('left', 'center'))
  

fig3d
ggsave(path = paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/svg/'), filename = 'Fig3D.svg', 
       plot = fig3d, device = 'svg', units = 'in', width = 3.1, height = 3.1)

pdf(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/Fig3D.pdf'))
print(fig3d)
dev.off()


### Clone count
lz <- as.numeric(sdat@meta.data[which(sdat@meta.data$cluster_names_merged == 'LZ-like' & !is.na(sdat@meta.data$clone_count)), 'clone_count'])
rm <- as.numeric(sdat@meta.data[which(sdat@meta.data$cluster_names_merged == 'RM' & !is.na(sdat@meta.data$clone_count)), 'clone_count'])
am <- as.numeric(sdat@meta.data[which(sdat@meta.data$cluster_names_merged == 'AM' & !is.na(sdat@meta.data$clone_count)), 'clone_count'])
ks.test(lz, rm, alternative = 'two.sided')
ks.test(lz, am, alternative = 'two.sided')


plot_datF <- as.data.frame(table(paste(sdat$cluster_names_merged, sdat@meta.data[, 'clone_count'], sep = '_')))
plot_datF[, c('cluster_names_merged', 'clone_count')] <- as.data.frame(t(as.data.frame(lapply(plot_datF$Var1, function(x) strsplit(as.character(x), '_')[[1]]))))
### Remove NAs
plot_datF <- plot_datF[!(is.na(plot_datF[, 'clone_count']) | plot_datF[, 'clone_count'] == 'NA'), ]
plot_datF$cluster_names_merged <- factor(plot_datF$cluster_names_merged, levels = names(cluster_key_merged))

n <- length(unique(plot_datF[, 'clone_count']))
color_key <- subject_key[1:n]
ns <- unique(plot_datF[, 'clone_count'])
ns <- ns[order(ns)]
names(color_key) <- ns

colnames(plot_datF) <- c('Var1', 'Freq', 'cluster_names_merged', 'Var2')

fig3f <- ggplot(data = plot_datF, aes(fill = Var2, x = cluster_names_merged, y = Freq)) + xlab('Cluster') + 
  fig_theme +
  geom_bar(position = 'fill', stat = "identity") + 
  scale_fill_manual(values = color_key, name = 'n cells per\nlineage') +
  ylab('Frequency') +
  theme(panel.border = element_blank()) +
  theme(legend.spacing = unit(0,'points'),
        legend.box.spacing = unit(0, 'points'),
        axis.ticks.length = unit(0.2, 'points'),
        axis.text.x = element_text(angle = 0),
        axis.title.x=element_blank())
  
fig3f
ggsave(path = paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/svg/'), filename = 'Fig3F.svg', 
       plot = fig3f, device = 'svg', units = 'in', width = 3.1, height = 3.1)


pdf(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/Fig3E.pdf'), width = 5)
print(fig3e)
dev.off()
### SHM
plot_datE <- sdat@meta.data[, c('cluster_names_merged', 'VH_SHM_PCT')]
plot_datE <- plot_datE[which(!is.na(sdat@meta.data[, 'VH_SHM_PCT'])), ]

fig3e <- ggplot(plot_datE, aes(x = cluster_names_merged, y = VH_SHM_PCT, fill = cluster_names_merged, linewidth = 0)) +
  fig_theme + 
  #stat_compare_means(aes(group = timepoint), label = "p.signif", hide.ns=TRUE, color='red') +
  geom_violin(draw_quantiles = 0.5, size = 0.3) +
  scale_fill_manual(values = cluster_key_merged, labels = names(cluster_key_merged), name = 'Cluster') +
  #scale_x_discrete(labels = names(cluster_key_merged)) +
  ylab('% nt SHM') +
  theme(axis.title.x = element_blank(), 
        legend.position = "none",
        #axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 0))

fig3e
ggsave(path = paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/svg/'), filename = 'Fig3E.svg', 
       plot = fig3e, device = 'svg', units = 'in', width = 3.1, height = 3.1)

pdf(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/Fig3F.pdf'))

fig3f

dev.off()


pdf(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/Timepoint_bars.pdf'))
for(c in cols_factor){
  ### Frequency of epitopes by cluster_names
  plot_dat <- as.data.frame(table(paste(sdat$timepoint, sdat@meta.data[, c], sep = '_')))
  plot_dat[, c('timepoint', c)] <- as.data.frame(t(as.data.frame(lapply(plot_dat$Var1, function(x) strsplit(as.character(x), '_')[[1]]))))
  ### Remove NAs
  plot_dat <- plot_dat[!(is.na(plot_dat[, c]) | plot_dat[, c] == 'NA'), ]
  plot_dat$timepoint <- factor(plot_dat$timepoint, levels = c('wk2', 'wk6'))
  n <- length(unique(plot_dat[, c]))
  
  color_key <- subject_key[1:n]
  names(color_key) <- unique(plot_dat[, c])
  
  colnames(plot_dat) <- c('Var1', 'Freq', 'timepoint', 'Var2')
  p <- ggplot(data = plot_dat, aes(fill = Var2, x = timepoint, y = Freq)) +
    geom_bar(position = 'fill', stat = "identity")
  p <- p + scale_fill_manual(values = color_key, name = c)
  print(p)
}

for(c in cols_num){
  plot_dat <- sdat@meta.data[, c('timepoint', c)]
  plot_dat <- plot_dat[which(!is.na(sdat@meta.data[, c])), ]
  colnames(plot_dat) <- c('timepoint', 'Var1')
  
  mu <- ddply(plot_dat, "timepoint", summarise, grp.mean=mean(Var1))
  p <- ggplot(plot_dat, aes(x = Var1, color = timepoint)) + 
    geom_density() + scale_color_manual(values = time_key) + labs(x = c) + 
    geom_vline(data = mu, aes(xintercept = grp.mean, color = timepoint), linetype = "dashed")
  print(p)
}
dev.off()

pdf(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/reads_by_timepoint.pdf'))
cls <- c('N_reads', 'N_mapped')
for(c in cls){
  plot_dat <- sdat@meta.data[, c('timepoint', c)]
  plot_dat <- plot_dat[which(!is.na(sdat@meta.data[, c])), ]
  colnames(plot_dat) <- c('timepoint', 'Var1')
  plot_dat$Var1 <- as.numeric(plot_dat$Var1)
  
  mu <- ddply(plot_dat, "timepoint", summarise, grp.mean = mean(Var1))
  p <- ggplot(plot_dat, aes(x = Var1, color = timepoint)) + 
    geom_density() + scale_color_manual(values = time_key) + labs(x = c) + 
    geom_vline(data = mu, aes(xintercept = grp.mean, color = timepoint), linetype = "dashed")
  print(p)
}
dev.off()

pdf(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/UMAP_SHM.pdf'))
FeaturePlot(sdat, features = 'VH_SHM_PCT', cols = c('pink', 'red'))
dev.off()

### Umap basics
DimPlot(sdat, group.by = 'cluster_names', label = F, reduction = umaps[1], cols = cluster_key, pt.size = 0.5) + labs(x= 'UMAP 1', y = 'UMAP 2', title = '')

pdf(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/Fig3A.pdf'))
fig3a <- DimPlot(sdat, group.by = 'cluster_names', label = F, reduction = umaps[1], cols = cluster_key, pt.size = 0.1, label.size = 1) + 
  fig_theme + 
  labs(x = 'UMAP 1', y = 'UMAP 2', title = '') +
  theme(plot.title=element_blank(),
        legend.title = element_blank(),
        legend.spacing = unit(5, 'points'), 
        legend.box.spacing = unit(0, 'points'),
        legend.key.size = unit(12, 'points'))
fig3a
ggsave(path = paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/svg/'), filename = 'Fig3A.svg', 
       plot = fig3a, device = 'svg', units = 'in', width = 3.5, height = 2.5)

DimPlot(sdat, group.by = 'cluster_names_merged', label = F, reduction = umaps[1], cols = cluster_key_merged, pt.size = 1.2) + labs(x= 'UMAP 1', y = 'UMAP 2', title = '')
dev.off()

umap_color1 <- 'timepoint'
DimPlot(sdat, group.by = umap_color1, label = F, reduction = umaps[1], cols = time_key)


umap_color1 <- 'predicted_label_main'
pal <- hue_pal()(length(unique(sdat@meta.data[, umap_color1])))
cells.highlight <- c()
cols.highlight <- pal[as.character(sdat@meta.data[, umap_color1])]
#DimPlot(sdat, group.by = umap_color1, label = T, reduction = umaps[1], cols = pal)
DimPlot(sdat, group.by = umap_color1, label = F, reduction = umaps[1], cols = pal)

table(sdat$predicted_label_fine, sdat$seurat_clusters)

to_display <- names(table(sdat$predicted_label_fine)[which(table(sdat$predicted_label_fine) > 10)])
sdat$predicted_label_fine<- ifelse(sdat$predicted_label_fine %in% to_display, sdat$predicted_label_fine, 'Other')
table(sdat$predicted_label_fine, sdat$seurat_clusters)

umap_color1 <- 'predicted_label_fine'
pal <- hue_pal()(length(unique(sdat@meta.data[, umap_color1])))
cells.highlight <- c()
cols.highlight <- pal[as.character(sdat@meta.data[, umap_color1])]
#DimPlot(sdat, group.by = umap_color1, label = T, reduction = umaps[1], cols = pal)
DimPlot(sdat, group.by = umap_color1, label = F, reduction = umaps[1], cols = pal)


### Prep fgsea data
### Read DE results
de_list <- lapply(de_files, function(de_file) read.table(de_file, header = T, check.names = F, stringsAsFactors = F, na.strings = c("", "NA"), sep = '\t'))
#names(de_list) <- ifelse(strats %in% names(cluster_names_merged), cluster_names_merged[strats], strats)
names(de_list) <- names(strats)
  
fgseas <- lapply(fgsea_files, function(fgsea_file) read.table(fgsea_file, header = T, check.names = F, stringsAsFactors = F, na.strings = c("", "NA"), sep = '\t'))
names(fgseas) <- names(strats)
### Count the x from this warnings:
### """There were x pathways for which P-values were not calculated properly due to unbalanced (positive and negative) gene-level statistic values. 
### For such pathways pval, padj, NES, log2err are set to NA. You can try to increase the value of the argument nPermSimple 
### (for example set it nPermSimple = 100000)"""
sapply(fgseas, function(fgsea) sum(is.na(fgsea$pval)))

### Remove 'All' (?)
de_list <- de_list[which(names(de_list) != 'All')]
fgseas <- fgseas[which(names(fgseas) != 'All')]
names(d_sets) <- sapply(names(d_sets), function(x) toupper(gsub(' ', '_', x)))
fgseas <- lapply(fgseas, function(fgsea) {
  fgsea$pathway <- sapply(fgsea$pathway, function(x) gsub('HALLMARK_', '', x))
  fgsea[, 'category'] <- d_sets[fgsea$pathway]
  fgsea$category <- ifelse(is.na(fgsea$category), 'Other', fgsea$category)
  return(fgsea)
})

pthresh <- 0.05
de_sig <- lapply(1:length(de_list), function(i) de_list[[i]][which(de_list[[i]]$p_val_adj < pthresh), 'gene_name'])
names(de_sig) <- names(de_list)
my_venn(names(de_list), de_sig, paste0('Timepoint DE N genes (adjp < ', pthresh, ')\n compared across cluster'))


fgsea_sig <- lapply(1:length(fgseas), function(i) fgseas[[i]][which(fgseas[[i]]$padj < pthresh), 'pathway'])
names(fgsea_sig) <- names(fgseas)
my_venn(names(fgseas), fgsea_sig, paste0('Timepoint GSEA N gene-sets (adjp < ', pthresh, ')\n compared across cluster'))

### GSEA dot plot
pthresh <- 0.05
pathways <- unique(unlist(lapply(fgseas, function(fgsea) fgsea[which(fgsea$padj < pthresh), 'pathway'])))
while(length(pathways) > 50){
  pthresh <- pthresh/10
  pathways <- unique(unlist(lapply(fgseas, function(fgsea) fgsea[which(fgsea$padj < pthresh), 'pathway'])))
}
### Correct if threshold reduced too much
if(length(pathways) == 0){
  pthresh <- pthresh * 10
  pathways <- unique(unlist(lapply(fgseas, function(fgsea) fgsea[which(fgsea$padj < pthresh), 'pathway'])))
}
plot_datG <- rbindlist(fgseas, idcol = 'cluster')
plot_datG <- plot_datG[which(plot_datG$pathway %in% pathways),]
plot_datG$cluster <- factor(plot_datG$cluster, levels = unique(cluster_names_merged))
### Remove insignificant
plot_datG <- plot_datG[which(plot_datG$padj < 0.05), ]

### Handle long pathway names

#plot_datG$pathway <- sapply(plot_datG$pathway, function(x) split_name(x, 80))
### manually eyeball gene names in the pathway names that aren't directly matched in the gtf file
gene_names <- c('TH17', 'RANKLRANK', 'IL12', 'TNFR2', 'NTHI', 'AP1', 'TAK1', 'IKK', 'SASP', 'TXA2',
                'STAT5', 'MTORC1', 'UV', 'TNFA', 'NFKB')
plot_datG$pathway <- sapply(plot_datG$pathway, function(x) abbrev_pathway_name(x, gtf, gene_names))
plot_datG <- as.data.frame(plot_datG)
### Order pathways
#a<-hclust(dist(plot_datG$NES))

pathways <- unique(plot_datG$pathway)
pathways <- pathways[order(sapply(pathways, function(pathway) mean(abs(plot_datG[which(plot_datG$pathway == pathway), 'NES']))), decreasing = T)]

plot_datG$pathway <- factor(plot_datG$pathway, levels = pathways)

#pathway<-plot_datG$pathway[which.max(nchar(plot_datG$pathway))]

legend_size <-7
x_size <- 10
y_size <- 9
#y_size <- rep(10, length(plot_datG$pathway))
#y_size[which(grepl('\n', plot_datG$pathway))] <- 6

fig3g <- ggplot(data = plot_datG, aes(y = pathway, x = cluster, size = -log10(padj), fill = NES)) +
  geom_point(shape = 21) + 
  scale_size_continuous(range = c(0,4)) +
  scale_fill_gradientn(colours = viridis(100)) +
  fig_theme +
  scale_x_discrete(position = "top") + 
  scale_y_discrete(position = "right") + 
  theme(axis.title.y = element_blank(),
        axis.title = element_blank()) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 0, size = x_size, color = "black"),
  #      axis.text.y = element_text(color = "black", size = y_size)) + 
  labs(size = '-log10 adj. p value', fill = 'Normalized\nEnrichment Score') +
  theme(legend.position = "bottom",
        legend.title = element_text(colour = "black"),#, size = legend_size),
        legend.spacing = unit(0, 'points'), legend.key.size = unit(8, 'points'), 
        legend.direction = 'horizontal', legend.box = 'horizontal',
        legend.box.spacing = unit(1, 'points'),
        legend.box.margin = margin(t = 0, r = -10, b = 0, l = 10, unit = "point")) +
  #theme(axis.text.y = element_text(size = 18, face = "bold", color = "black")) + 

  guides(shape = guide_legend(override.aes = list(size = 8))) + 
  guides(color = guide_legend(override.aes = list(size = 8))) + 
  theme(panel.spacing = unit(1, "lines")) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5), size = guide_legend(title.position = "top", title.hjust = 0.5))
  #  facet_grid(category, scales = "free", space = "free", switch = "both") +
  # theme(text = element_text(family = "Arial", color = "black"),
  #       strip.background = element_blank(),
  #       strip.text.y.left = element_text(angle = 0, size = 8, hjust=1, face = "bold"),
  #       strip.text.x = element_text(size = 8, face = "bold")
  #)

fig3g
ggsave(path = paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/svg/'), filename = 'Fig3G.svg', 
       plot = fig3g, device = 'svg', units = 'in', width = 3.3, height = 3.6)


pdf(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/Fig3G_legend.pdf'), width = 11)
fig3g
dev.off()

pdf(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/Fig3G.pdf'))
fig3g
dev.off()

png(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/Fig3G.png'))
fig3g
dev.off()

### Read counts
norm_counts <- as.data.frame(sdat@assays$RNA@data)
norm_counts <- mutate_all(norm_counts, function(x) as.numeric(x))

de_ids <- unique(unlist(lapply(names(de_list), function(x) row.names(de_list[[x]]))))
### Remove genes that weren't in all DE analyses
dat <-as.data.frame(lapply(names(de_list), function(x) de_ids %in% row.names(de_list[[x]])))
row.names(dat) <- de_ids          
de_ids <- de_ids[sapply(row.names(dat), function(rn) rowSums(dat[rn, ]) == length(de_list))]

pdf(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/Pathway_LogFC.pdf'))
### Timepoint GSEA results
for(geneset_name in pathways){
  print(paste0(geneset_name, ' Heatmap'))
  genes <- DB[[geneset_name]]
  genes <- genes[which(genes %in% gtf$gene_name)]
  column_title <- 'wk2 v wk6'
  row_title <- paste0(geneset_name, '\n(', length(genes), ' of ', length(DB[[geneset_name]]), ' genes)')
  ### Set names as the ID (Remove the "[1]" if better handled multiple genes per ID)
  names(genes) <- sapply(genes, function(x) row.names(gtf[which(gtf$gene_name == x)[1],]))
  ### Only keep those in the DE results (all)
  genes <- genes[which(names(genes) %in% de_ids)]
  print(length(genes) / length(DB[[geneset_name]]))
  ### fgsea
  if(geneset_name %in% fgseas[[1]]$pathway){
    ps <- sapply(names(de_list), function(i) signif(fgseas[[i]][which(fgseas[[i]]$pathway == geneset_name), 'padj'], digits = 2))
    column_labels <- paste0(names(de_list),'\n(adj p = ', ps, ')')
    ### Abbreivate if label is too long
    #column_labels <- ifelse(nchar(column_labels) > 15, gsub('cluster', '', column_labels), column_labels)
  } else{
    column_labels <- names(de_list)
  }
  fontsize <- heatmap_N_to_size(length(genes))
  ### genes is an array of gene_names (for plotting) with names as IDs as appear in de_list
  heatp <- my_heatmap(de_list, genes, categories = NA, fontsize = fontsize, column_title, row_title, pheatmap = F, col_clust = F, column_labels = column_labels)
  draw(heatp, show_heatmap_legend = T) 
  
  ### Dynamic pdf height
  pdf_height = heatmap_N_to_pdfsize(length(genes))
  pdf(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/pathways/', geneset_name,'.pdf'), height = pdf_height)
  p <- DotPlot(sdat, features = names(genes), group.by = 'timepoint') +
    RotatedAxis() + scale_x_discrete(labels = genes) + labs(y = 'timepoint') + coord_flip()#+ xlab(paste(cat, ' features')) + ggtitle(plot_title)
  print(p)
  dev.off()
}
dev.off()

### Individual pathways, expression split by cluster and timepoint. Would be great to display timepoint DE logFC somehow
for(geneset_name in pathways){
  genes <- DB[[geneset_name]]
  genes <- genes[which(genes %in% gtf$gene_name)]
  ### Set names as the ID (Remove the "[1]" if better handled multiple genes per ID)
  names(genes) <- sapply(genes, function(x) row.names(gtf[which(gtf$gene_name == x)[1],]))
  ### Only keep those in the DE results (all)
  genes <- genes[which(names(genes) %in% de_ids)]
  ### Dynamic pdf height
  pdf_height = heatmap_N_to_pdfsize(length(genes))
  pdf(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/pathways/', geneset_name,'.pdf'), height = pdf_height)
  p <- DotPlot(sdat, features = names(genes), group.by = 'cluster_names', split.by = 'timepoint') +
    RotatedAxis() + scale_x_discrete(labels = genes) + labs(y = 'timepoint') + coord_flip()#+ xlab(paste(cat, ' features')) + ggtitle(plot_title)
  print(p)
  dev.off()
}

### To look at individual genes epxression as violin plots per cluster, with overlayed DE p and logFC
gene_names <- c('CNP', 'ENO2', 'EPHX1', 'ATF3', 'TNF')
ids <- row.names(gtf[which(gtf$gene_name %in% gene_names), ])
#ids <- row.names(de_list[[1]])[c(1,3)]

pdf(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/Violins.pdf'))
for(id in ids){
  ptext <- paste0('adjp: ', sapply(names(de_list), function(x) signif(de_list[[x]][id, 'p_val_adj'], digits = 3)))
  names(ptext) <- names(de_list)
  fctext <- paste0('FC: ',sapply(names(de_list), function(x) signif(de_list[[x]][id, 'avg_log2FC'], digits = 3)))
  names(fctext) <- names(de_list)
  ymax <- max(sdat@assays$RNA@data[id, ]) * 1.08
  ytext1 <- max(sdat@assays$RNA@data[id, ]) * 1.01
  ytext2 <- max(sdat@assays$RNA@data[id, ]) * 1.04
  
  p <- VlnPlot(sdat, features = id, group.by = 'cluster_names', split.by = 'timepoint', y.max = ymax) + ggtitle(gtf[id, 'gene_plot_name']) + ylab('TPM') +
    annotate('text', x = 1:length(ptext), y = ytext1, label = ptext) +
    annotate('text', x = 1:length(fctext), y = ytext2, label = fctext)
  print(p)
}
dev.off()

id <- row.names(gtf[which(gtf$gene_name %in% c('TNF')), ])
p <- VlnPlot(sdat, features = id, group.by = 'cluster_names_merged', split.by = 'timepoint', y.max = ymax) + ggtitle(gtf[id, 'gene_plot_name']) + ylab('TPM')
p
expression_plots2(sdat, dotplot_genes, gtf, features, c('fig3b'), 'cluster_names', cat_label = F, scale = T, scale_by = 'size')

p <- DotPlot(sdat, features = id, group.by = 'cluster_names_merged', scale = T, scale.by = 'radius', split.by = 'timepoint') +
  RotatedAxis() + scale_x_discrete(labels = 'TNF') + coord_flip() + 
  theme(axis.title = element_blank())
p

sdat$Cluster_TP <- paste0(sdat$cluster_names_merged, '_', sdat$timepoint)
sdat$Cluster_TP <- gsub('_wk', ' week ', sdat$Cluster_TP)
sdat$Cluster_TP <- factor(sdat$Cluster_TP, levels = c('AM week 6', 'AM week 2', 'RM week 6', 'RM week 2', 'LZ-like week 6', 'LZ-like week 2'))

p <- DotPlot(sdat, features = id, group.by = 'Cluster_TP', scale = T, scale.by = 'radius') +
  RotatedAxis() + 
  scale_x_discrete(labels = 'TNF') + 
  #coord_flip() + 
  theme(text = element_text(family = font, color = "black"),
        axis.title = element_blank(),
        #axis.title = element_text(family = font, size = 8, color = "black"),
        axis.text  = element_text(family = font, size = 12, color = "black"),
        legend.title = element_text(family = font, size = 12, color = "black"),
        legend.text  = element_text(family = font, size = 10, color = "black"),
        legend.spacing = unit(5,'points'),
        legend.key.size = unit(12, 'points'))
        #legend.box.spacing = unit(0, 'points'),
        #axis.ticks.length = unit(1, 'points')) + 
p
ggsave(path = paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/svg/'), filename = 'Supplemental.svg', 
       plot = p, device = 'svg', units = 'in', width = 5, height = 3)

ggplot(data = input_all, aes(y = hallmark_pathway, x = cell_type, size = -log10(padj), fill = NES)) +
  geom_point(shape = 21) + 
  scale_size_continuous(range = c(0,4)) +
  scale_fill_gradientn(colours=viridis(100)) +
  theme_bw() +
  scale_x_discrete(position = "top") + 
  scale_y_discrete(position = "right") + 
  theme(axis.text.x=element_text(angle = 45, hjust = 0)) + 
  theme(axis.title.y = element_blank()) +
  theme(legend.position = "bottom") + 
  labs(size = '-log10(adjusted P value)', fill = 'Normalized Enrichment Score') +
  theme(legend.title = element_text(face = "bold",colour = "black", size = 10)) +
  #theme(axis.text.y = element_text(size = 18, face = "bold", color = "black")) + 
  theme(axis.text.y = element_text(color = "black")) + 
  theme(axis.text.x = element_text(size = 6, color = "black")) + 
  guides(shape = guide_legend(override.aes = list(size = 5))) + 
  guides(color = guide_legend(override.aes = list(size = 5))) + 
  theme(legend.title = element_text(size = 6), legend.text = element_text(size = 6),
        panel.spacing = unit(0, "lines"))+
  facet_grid(category~timepoints, scales = "free", space = "free", switch = "both") +
  theme(text = element_text(family = "Arial", color="black"),
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 10, hjust=1, face = "bold"),
        strip.text.x = element_text(size = 8, face = "bold")
  )
# #pathway <- 'UV_RESPONSE_UP'
# pathway <- 'GSE14769_UNSTIM_VS_120MIN_LPS_BMDM_DN'
# #cell_type <- 'Resting'
# cell_type <- 'Activated'
# 
# fgsea <- fgseas[[cell_type]]
# gene_names <- strsplit(fgsea[which(fgsea$pathway == pathway), 'leadingEdgeStr'], ',')[[1]]
# gene_ids <- gtf[which(gtf$gene_name %in% gene_names), 'gene_id']
# ### Plot logFC of these genes
# 
# ### Plot norm_counts of these genes divided by test
# sdat_sub <- subset(sdat, cells = row.names(sdat@meta.data[which(sdat$cluster_names == cell_type),]))
# sdat_sub[[pathway]] <- sapply(colnames(sdat_sub), function(c) mean(sdat_sub@assays$RNA@data[gene_ids, c]))
# VlnPlot(sdat_sub, features = )

DimPlot(sdat, group.by = 'cluster_names', label = F, reduction = umaps[1], cols = cluster_key, pt.size = 1.2) + labs(x= 'UMAP 1', y = 'UMAP 2', title = '')
expression_plots2(sdat, dotplot_genes, gtf, features, c('fig3b'), 'cluster_names', cat_label = F, scale = T, scale_by = 'size')
fig3c + theme(legend.position = "bottom")
fig3d + theme(legend.position = "bottom")
fig3e
fig3f
fig3g

#dotplot theme
# F -  check y label isn't cut off
# GGh4x has "force panel size"
# Supplemental - Vln plot or Dot plot of TNF expression by cluster and timepoint
# Reread method section



# ggsave(path = paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/svg/'), filename = 'Fig3A.svg', 
#        plot = fig3a, device = 'svg', units = 'in', width = 3.1, height = 2.5)
# ggsave(path = paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/svg/'), filename = 'Fig3B.svg', 
#        plot = fig3b, device = 'svg', units = 'in', width = 4.0, height = 2.5)
# ggsave(path = paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/svg/'), filename = 'Fig3C.svg', 
#        plot = fig3c, device = 'svg', units = 'in', width = 3.0, height = 3.0)
# ggsave(path = paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/svg/'), filename = 'Fig3D.svg', 
#        plot = fig3d, device = 'svg', units = 'in', width = 2.7, height = 3.0)
# ggsave(path = paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/svg/'), filename = 'Fig3E.svg', 
#        plot = fig3e, device = 'svg', units = 'in', width = 2.22, height = 2.8)
# ggsave(path = paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/svg/'), filename = 'Fig3F.svg', 
#        plot = fig3f, device = 'svg', units = 'in', width = 2.45, height = 2.65)
# ggsave(path = paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/svg/'), filename = 'Fig3G.svg', 
#        plot = fig3g, device = 'svg', units = 'in', width = 2.77, height = 2.9)
# 

# svg(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/svg/Fig3A.svg'), width = 3.1, height = 2.5)
# fig3a
# dev.off()
# 
# svg(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/svg/Fig3B.svg'), width = 4, height = 2.5)
# fig3b
# dev.off()
# 
# svg(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/svg/Fig3C.svg'), width = 3, height = 3)
# fig3c
# dev.off()
# 
# svg(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/svg/Fig3D.svg'), width = 2.7, height = 3)
# fig3d
# dev.off()
# 
# svg(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/svg/Fig3E.svg'), width = 2.22, height = 2.8)
# fig3e
# dev.off()
# 
# svg(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/svg/Fig3F.svg'), width = 2.45, height = 2.65)
# fig3f
# dev.off()
# 
# svg(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/svg/Fig3G.svg'), width = 2.77, height = 2.9)
# fig3g
# dev.off()
# 
# svg(paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2022620_857.1/results/', qc_name, '/svg/Fig3G_legend.svg'), width = 7, height = 4)
# fig3g
# dev.off()
