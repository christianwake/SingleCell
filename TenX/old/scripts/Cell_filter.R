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
library('vegan')
library('data.table')

source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/Utility_functions.R')

args = commandArgs(trailingOnly=TRUE)

sdat_file <- args[1]
gtf_file <- args[2]
key_words <- args[3]
fraction_cells <- args[4]
nFeature_RNA <- args[5]
MT_prop <- args[6]
nCount_RNA <- args[7]
Shannon_diversity <- args[8]
out_rds <- args[9]
out_pdf <- args[10]
out_txt1 <- args[11]
out_txt2 <- args[12]

# sdat_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021618_galt/data/Seurat_raw.RDS'
# gtf_file <- '/hpcdata/vrc/vrc1_data/douek_lab/reference_sets/tenX/Homo_sapiens.GRCh38.93/GRCh38_protein_coding_only/genes/genes.gtf'
# key_words <- '^IG,^RPS'
# fraction_cells <- NA
# nFeature_RNA <- NA
# MT_prop <- NA
# nCount_RNA <- NA
# Shannon_diversity <- NA
# out_rds <- "/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021618_galt/results/RNA_cell_filtered.RDS"
# out_pdf <- "/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021618_galt/results/MT_nFeatures.pdf"
# out_txt1 <- "/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021618_galt/results/Excluded_genes.txt"
# out_txt2 <-  "/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021618_galt/results/Filters_Ncells.txt"

### Set defaults
if(is.na(fraction_cells) | fraction_cells == ''){
  fraction_cells <- '<0.03'
}
if(is.na(nFeature_RNA) | nFeature_RNA == ''){
  nFeature_RNA <- '1.5(SD),3(SD)'
}
if(is.na(nCount_RNA) | nCount_RNA == ''){
  nCount_RNA <- '1.5(SD),3(SD)'
}
if(is.na(MT_prop) | MT_prop == ''){
  MT_prop <- '0,0.4'
}
if(is.na(Shannon_diversity) | Shannon_diversity == ''){
  Shannon_diversity <- '3(SD),3(SD)'
}
print(key_words)
print(fraction_cells)
print(nFeature_RNA)
print(MT_prop)

sdat <- readRDS(sdat_file)
DefaultAssay(sdat) <- 'RNA'

all_names <- row.names(sdat@assays$RNA@counts)
key_words <- strsplit(key_words, ',')[[1]]
exclude_gene_names <- all_names[sapply(all_names, function(x) any(sapply(key_words, function(y) grepl(pattern = y, x = x, perl = T))))]
exclude_gene_names <- sort(exclude_gene_names)
write.table(exclude_gene_names, out_txt1, quote = F, sep = ',', row.names = F, col.names = F)

### Determine which genes are mitochondrial by gtf if it is input. Otherwise, by gene name pattern '^MT'.
if(!is.na(gtf_file) & gtf_file != ''){
  gtf <- read_gtf(gtf_file, 'gene', 'gene_name')
  sdat[['MT_sum']] <- sum_chromosome(sdat@assays$RNA@counts, chr ='MT', gtf, 'gene_name')
} else {
  mts <- row.names(sdat@assays$RNA)[grepl('^MT-', row.names(sdat@assays$RNA))]
  sdat[['MT_sum']] <- colSums(sdat@assays$RNA@counts[mts, ])
}
sdat[['MT_prop']] <- sdat[['MT_sum']]/sdat[['nCount_RNA']]
### Shannon diversity
sdat[['Shannon_diversity']] <- diversity(x = sdat@assays$RNA@counts, index = 'shannon', MARGIN = 2)

pdf(out_pdf)
VlnPlot(sdat, features = "nFeature_RNA")
VlnPlot(sdat, features = 'MT_prop')
VlnPlot(sdat, features = 'nCount_RNA')
VlnPlot(sdat, features = 'Shannon_diversity')

cor(sdat$nFeature_RNA, sdat$MT_prop)
plot(sdat$nFeature_RNA, sdat$MT_prop, xlab = 'N genes', ylab = 'MT sum / nCount_RNA')
#FeatureScatter(sdat, feature1 = "nFeature_RNA", feature2 = "MT_sum", group.by = "sample")

### Get array of lower,upper bound thresholds from the comma-delimited string. Options are flat values, or if in values of standard deviation, with a trailing (SD)
MT_thresh <- threshold_string_seurat(sdat, MT_prop, 'MT_prop')
NFeat_thresh <- threshold_string_seurat(sdat, nFeature_RNA, 'nFeature_RNA')
NCount_thresh <- threshold_string_seurat(sdat, nCount_RNA, 'nCount_RNA')
div_thresh<- threshold_string_seurat(sdat, Shannon_diversity, 'Shannon_diversity')

plot(sdat$nFeature_RNA, sdat$MT_prop, xlab = 'N genes', ylab = 'Mitochondrial proportion') +
  abline(v = NFeat_thresh[1], col = 'red') +
  abline(v = NFeat_thresh[2], col = 'red') +
  abline(h = MT_thresh[2], col = 'red')
plot(sdat$nCount_RNA, sdat$Shannon_diversity, xlab = 'N reads', ylab = 'Shannon Diversity') +
  abline(v = NCount_thresh[1], col = 'red') +
  abline(v = NCount_thresh[2], col = 'red') +
  abline(h = div_thresh[1], col = 'red') +
  abline(h = div_thresh[2], col = 'red')

txt <- as.data.frame(matrix(ncol = 2, nrow = 9))
colnames(txt) <- c('filter', 'n')
txt[, 'filter'] <- c('Original', 'nFeature_low', 'nFeature_high', 'MT_low', 'MT_high', 'nCount_low', 'nCount_high', 'Div_low', 'Div_high')
row.names(txt) <- txt[, 'filter']

sdat_full <- sdat
txt['Original', 'n'] <- dim(sdat_full)[2]
sdat <- subset(sdat, subset = nFeature_RNA >  NFeat_thresh[1])
txt['nFeature_low', 'n'] <- dim(sdat)[2]
sdat <- subset(sdat, subset = nFeature_RNA <  NFeat_thresh[2])
txt['nFeature_high', 'n'] <- dim(sdat)[2]

sdat <- subset(sdat, subset = MT_prop >  MT_thresh[1])
txt['MT_low', 'n'] <- dim(sdat)[2]
sdat <- subset(sdat, subset = MT_prop < MT_thresh[2])
txt['MT_high', 'n'] <- dim(sdat)[2]

sdat <- subset(sdat, subset = nCount_RNA >  NCount_thresh[1])
txt['nCount_low', 'n'] <- dim(sdat)[2]
sdat <- subset(sdat, subset = nCount_RNA < NCount_thresh[2])
txt['nCount_high', 'n'] <- dim(sdat)[2]

sdat <- subset(sdat, subset = Shannon_diversity >  div_thresh[1])
txt['Div_low', 'n'] <- dim(sdat)[2]
sdat <- subset(sdat, subset = Shannon_diversity < div_thresh[2])
txt['Div_high', 'n'] <- dim(sdat)[2]


plot(sdat$nFeature_RNA, sdat$MT_prop, xlab = 'N genes', ylab = 'MT proportion', main = 'post-filtering')
plot(sdat$nCount_RNA, sdat$Shannon_diversity, xlab = 'N reads', ylab = 'Shannon Diversity', main = 'post-filtering')

VlnPlot(sdat, features = 'nCount_RNA', group.by = 'tissue')
VlnPlot(sdat, features = 'nCount_RNA', group.by = 'celltype')
#VlnPlot(sdat, features = 'nCount_RNA', group.by = 'Lane')
VlnPlot(sdat, features = 'nFeature_RNA', group.by = 'tissue')
VlnPlot(sdat, features = 'nFeature_RNA', group.by = 'celltype')
#VlnPlot(sdat, features = 'nFeature_RNA', group.by = 'Lane')
VlnPlot(sdat, features = 'MT_prop', group.by = 'tissue')
VlnPlot(sdat, features = 'MT_prop', group.by = 'celltype')
#VlnPlot(sdat, features = 'MT_prop', group.by = 'Lane')
dev.off()

write.table(txt, out_txt2, sep = ',', quote = F, row.names = F, col.names = T)

#sdat_removed <- subset(sdat_full, cells = colnames(sdat_full)[!colnames(sdat_full) %in% colnames(sdat)])

### Fraction kept by each sample. The first two lines are necessary in the case that an entire sample is filtered
# temp <- table(sdat$sample)
# temp[unique(sdat_full$sample)[which(!unique(sdat_full$sample) %in% unique(sdat$sample))]] <- 0
# temp/table(sdat_full$sample)[names(temp)]

saveRDS(sdat, file = out_rds)

