library('sys')
library('ggplot2')
library('stringr')
library('GenomicRanges')
library('dplyr')
library('readr')
library('effsize')
library('data.table')
library('biomaRt')
library('fgsea')
library('GSEABase')
library('readxl')

source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/DE_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')

if(interactive()){
  project <- '2021600_kristin'
  qc_name <- 'Run2023-05-14'
  test <- 'RNA_clusters'
  value1 <- '2'
  value2 <- '3'
  strat <- 'All'
  strat_values <- 'All'
  species <- 'hsapiens'
  sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Filtered_clustered.RDS')
  
  # project <- '2022620_857.1'
  # qc_name <- '2022-12-09'
  # test <- 'timepoint'
  # strat <- 'seurat_clusters-1'
  # species <- 'mmulatta'
  # sdat_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Filtered.RDS')
  
  # res_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', 
  #                    qc_name, '/DE/', test, '/', strat, '/DE_results.tsv')
  res_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', 
                     qc_name, '/DE/', test, '/', strat, '/', test, '-', value1, '-', value2, 
                     '_DE_within_', strat, '-', strat_values, '.tsv')
  
  #res_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/DE/', test, '/', strat, '/DE_results.tsv')
  exclude_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/Excluded_genes.txt')
  gtf_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/data/gtf.RDS')
  #gmt_file <- '/hpcdata/vrc/vrc1_data/douek_lab/wakecg/genesets/c2.cp.v7.2.symbols.gmt'
  #gmt_file <- '/hpcdata/vrc/vrc1_data/douek_lab/wakecg/genesets/h.all.v2022.1.Hs.symbols.gmt'
  gmt_file <- '/hpcdata/vrc/vrc1_data/douek_lab/wakecg/genesets/c7.all.v2022.1.Hs.symbols.gmt'
  out_tsv <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/DE/', test, '/', strat, '/fgsea.tsv')
  out_pdf <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, '/results/', qc_name, '/DE/', test, '/', strat, '/fgsea.pdf')
  custom_sets <- NA
} else{
  args = commandArgs(trailingOnly=TRUE)
  sdat_file <- args[1]
  res_file <- args[2]
  exclude_file <- args[3]
  gtf_file <- args[4]
  gmt_file <- args[5]
  species <- args[6]
  out_tsv <- args[7]
  #out_pdf <- args[8]
  custom_sets <- args[8]
}

minSize <- 5
### If the file exists and is not empty, read it
if(file.exists(exclude_file) & file.info(exclude_file)$size != 0){
  exclude_gene_names <- read.table(exclude_file)[,1]
} else{
  exclude_gene_names <- c()
}

### Read counts
sdat <- readRDS(sdat_file)
norm_counts <- as.data.frame(sdat@assays$RNA@data)
norm_counts <- norm_counts[which(!(row.names(norm_counts) %in% exclude_gene_names)),]
rm(sdat)

### Convert to numeric
#norm_counts <- mutate_all(norm_counts, function(x) as.numeric(x))

### Do direct matching from gtf
if(grepl('\\.gtf', gtf_file)){
  gtf <- read_gtf(gtf_file, feature_type = 'gene', atts_of_interest = c('gene_id', 'gene_name',  'gene_biotype'))
} else if(grepl('\\.RDS', gtf_file)){
  gtf <- readRDS(gtf_file)
} else{
  warning('No gtf file')
}
row.names(gtf) <- gtf$gene_id

### Read DE results
res <- read.table(res_file, header = T, check.names = F, stringsAsFactors = F, na.strings = c("", "NA"), sep = '\t')

### Run biomaRt to get gene name information from ensembl IDs.
res$ensembl <- row.names(res)
#colnames(res) <- gsub('gene_name', 'gene_name_gtf', colnames(res))
if(species != 'hsapiens'){
  ### BiomaRt, for gene names
  id_key <- run_biomaRt(res, to_merge = F, species = species)
  ### (For now) Only keep one-to-one orthologs.
  id_key <- id_key[which(id_key$ortholog_type == 'ortholog_one2one'),]
  id_key[id_key == ''] <- NA
  row.names(id_key) <- id_key$ens_short
  colSums(is.na(id_key))
  ### Prioritize Gene names from biomart human orthologs
  res[row.names(id_key), 'gene_name'] <- id_key$gene_name
  ### But when that is NA, use gene name from the other species' gtf
  res[which(is.na(res$gene_name) & !is.na(res$gene_name_gtf)), 'gene_name'] <- res[which(is.na(res$gene_name) & !is.na(res$gene_name_gtf)), 'gene_name_gtf']
  res[, 'gene_name_plot'] <- res$gene_name
  res[which(is.na(res$gene_name)), 'gene_name_plot'] <- row.names(res[which(is.na(res$gene_name)), ])
  ### Adjust it if there are duplicates in 'gene_name_plot'
  dups <- names(table(res$gene_name_plot))[which(table(res$gene_name_plot) > 1)]
  res[which(res$gene_name_plot %in% dups), 'gene_name_plot'] <- 
    paste0(res[which(res$gene_name_plot %in% dups), 'gene_name_plot'], '(', row.names(res[which(res$gene_name_plot %in% dups), ]), ')')
} else{
  id_key <- gtf[, c('gene_id', 'gene_name')]
  colnames(id_key) <- c('gene_id', 'hsapiens_name')
  if('gene_name_gtf' %in% colnames(res)){
    res$gene_name <- res$gene_name_gtf
  } else{
    res$gene_name <- id_key[row.names(res), 'hsapiens_name']
    res$gene_name_plot <- res$gene_name
  }
}

# biomart_results_file <-  '/hpcdata/vrc/vrc1_data/douek_lab/wakecg/2020213_NHP857.1/Biomart_output.RDS'
# ids <- row.names(sdat)[which(!(row.names(sdat) %in% exclude_gene_names))]
# id_key <- readRDS(biomart_results_file)
# ### (For now) Only keep one-to-one orthologs.
# id_key <- id_key[which(id_key$ortholog_type == 'ortholog_one2one'),]
# id_key[id_key == ''] <- NA
# row.names(id_key) <- id_key$ens_short
# colSums(is.na(id_key))
# 
# ### Add gene_name's from id_key (biomart)
# res[ , 'gene_name'] <- id_key[row.names(res), 'hsapiens_name']
# res <- res[which(!is.na(res$gene_name)),]
# 
# rnk <- seurat_wilcox_Res2Rank(res)
# ranks <- setNames(rnk[, 'logp'], rnk$ID)

### Remove those without gene name annotation
res <- res[which(!(is.na(res$gene_name))),] 
dup_names <- names(table(res$gene_name))[which(table(res$gene_name) != 1)]
### Remove duplicate gene names
res <- res[which(!res$gene_name %in% dup_names),]
#rnk <- deseqRes2Rnk(res, by = '-log(p)', id = 'gene_name')
#rnk <- rnk[order(rnk$NegLogP, decreasing = T),]
#ranks <- setNames(rnk[, 'NegLogP'], rnk$ID)
rnk <- seurat_wilcox_Res2Rank(res)
rnk <- rnk[order(rnk$logp, decreasing = T),]
ranks <- setNames(rnk[, 'logp'], rnk$ID)
### Remove and send warning if any are NA, INF or -INF
if(any(is.na(ranks))){
  warning(paste0(sum(is.na(ranks)), ' features are NA logp and must be removed'))
  ranks <- ranks[!is.na(ranks)]
}
if(any(ranks == '-Inf')){
  warning(paste0(sum(ranks == '-Inf'), ' features are NA logp and must be removed'))
  ranks <- ranks[ranks != '-Inf']
}
if(any(ranks == 'Inf')){
  warning(paste0(sum(ranks == 'Inf'), ' features are NA logp and must be removed'))
  ranks <- ranks[ranks != 'Inf']
}
### MSigDB
genesets <- getGmt(gmt_file)
DB <- gmtPathways(gmt_file)
### Custom sets
goal_set <- res$gene_name_gtf[!is.na(res$gene_name_gtf)]
greeks <- c('Α', 'Β', 'Ε', 'Γ', 'Κ')

if(custom_sets != '' & !is.na(custom_sets)){
  sheets <- excel_sheets(path = custom_sets)
  ### When biomart is NA, it runs biomaRt_aliases(cust) to get it
  cus <- lapply(sheets, function(set_name) get_custom_geneset(custom_sets, set_name, goal_set, res, biomart = NA)[['genes']])
  names(cus) <- sheets
  ### Add custom sets to the DB object
  DB <- c(DB, cus)
} 

### Genes from all  gene sets
#genes <- unique(unlist(lapply(genesets, function(x) geneIds(x))))
genes <- unique(unlist(DB))

### Mitigate sample source bias by making sure that we only have genes with at least one count
a_ids <- row.names(norm_counts)[which(rowSums(norm_counts) > 0)]
a_names <- id_key[a_ids, 'hsapiens_name']
a_names <- a_names[which(!is.na(a_names))]
rm(norm_counts)
paste0('Of ', length(a_ids), ' ids, ', length(a_names), ' are in ', species, '-human biomaRt, and ', sum(a_names %in% genes), ' of those are in human ', gmt_file)

#DB <- DB[grepl('BIOCARTA', names(DB))]
# redundancy <- unlist(fgseaRes[which(fgseaRes$pathway %in% c('BIOCARTA_P38MAPK_PATHWAY', 'BIOCARTA_TCR_PATHWAY')), 'leadingEdge'])
# le <- unlist(fgseaRes[, 'leadingEdge'])
# a <- sapply(redundancy, function(x) sum(le == x))
# b <- sapply(unique(le), function(x) sum(le == x))
# mean(a)
# mean(b[!names(b) %in% names(a)])

print('Beginning fgsea')
fgseaRes <- fgsea(pathways = DB, stats = ranks, minSize = minSize, maxSize = 500, nPermSimple = 10000)
print('test1')
fgseaRes <- as.data.frame(fgseaRes)
print('test2')
fgseaRes <- fgseaRes[order(fgseaRes$pval), ]
print('test3')
### Convert the leadingEdge list of gene names into a comma delimited string, in order to print. 
fgseaRes[, 'leadingEdgeStr'] <- apply(fgseaRes, 1, function(x) gsub(', ', ',', toString(x['leadingEdge'][[1]])))
### Subset significant results
print('test4')
sig <- fgseaRes[which(fgseaRes$padj < 0.05), ]
print('test5')
min(fgseaRes$pval)
num_sig <- length(sig$padj)
#print(paste0(num_sig, ' padj 0.05 significant GSEA results for -log(p) sorted ', names(res_file)))
print('test6')

### Print significant GSEA results to a file
write.table(fgseaRes[,c("pathway","pval","padj","ES","NES","size","leadingEdgeStr")], out_tsv, sep = '\t', row.names = F, quote = F)

#fgseaRes <- fgseaRes[, which(colnames(fgseaRes) != 'leadingEdge')]
#WriteXLS(ExcelFileName = out_xls, x = 'fgseaRes', row.names = F, col.names = T)

### Plot of the most significant result
# top <- fgseaRes[order(fgseaRes$pval, decreasing = F),][1, 'pathway']
# ptext <- as.character(signif(fgseaRes[order(fgseaRes$pval, decreasing = F),][2, 'padj'][[1]], digits = 3 ))
# pdf(out_pdf)
# ### Top 10 most significant results
# for(i in 1:10){
#   ### Plot of the result
#   top <- fgseaRes[order(fgseaRes$pval, decreasing = F),][i, 'pathway']
#   ptext <- as.character(signif(fgseaRes[order(fgseaRes$pval, decreasing = F),][i, 'padj'][[1]], digits = 3 ))
#   print(plotEnrichment(DB[[top]], ranks) + labs(title=top, subtitle = paste0(names(res_file),', padj: ',ptext)))
# }
# dev.off() 
