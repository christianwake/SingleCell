### This script assumes you've already completed bcl2fastq, Trimmomatic, STAR and maybe BALDR
library('sys')
library('viridis')
library('data.table')
library('PKI')
library('stringr')
library('stringi')
library('biomaRt')
library('fgsea')
library('GSEABase')

source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/Utility_functions.R')

if(interactive()){
  species <- 'hsapiens'
  sdat_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021600_kristin/results/Run2022/PostQC3.RDS'
  gtf_file <- '/hpcdata/vrc/vrc1_data/douek_lab/reference_sets/tenX/Homo_sapiens.GRCh38.93/GRCh38_protein_coding_only/genes/genes.gtf'
  gmt_file <- '/hpcdata/vrc/vrc1_data/douek_lab/wakecg/c2.cp.v7.2.symbols.gmt'
  gene_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021600_kristin/220921_transcriptome_analysis.xlsx'
  out_file <- '/hpcdata/vrc/vrc1_data/douek_lab/wakecg/data/Gene_names.RDS'
} else{
  args = commandArgs(trailingOnly=TRUE)
  ### snakemake input
  species <- args[1]
  sdat_file <- args[2]
  gtf_file <- args[3]
  gmt_file <- args[4]
  gene_file <- args[5]
  out_file <- args[6]
}
greeks <- c('Α', 'Β', 'Ε', 'Γ', 'Κ')

### Read Seurat
sdat <- readRDS(sdat_file)
DefaultAssay(sdat) <- 'RNA'

### Read gtf info
if(grepl('.gtf', gtf_file)){
  gtf <- read_gtf(gtf_file, feature_type = 'gene', atts_of_interest = c('gene_id', 'gene_name',  'gene_biotype'))
} else if(grepl('.RDS', gtf_file)){
  gtf <- readRDS(gtf_file)
}

### Read a gene set file in gmt format
genesets <- getGmt(gmt_file)
DB <- gmtPathways(gmt_file)
### Read cutom gene set info
genes <- read_excel(gene_file, sheet = 1)
colnames(genes) <- c('gene_name', 'category')

### Instantiate results by set2, a dataframe of the gene names that need a match to those in sdat
set2 <- unique(c(genes$gene_name, unique(unlist(DB))))
res <- as.data.frame(set2)
colnames(res) <- c('gene_name')
row.names(res) <- res[,1]

### What they need to match to
set1 <- row.names(sdat)
rm(sdat)

### First, does set2 match to set1 by the gtf file
key <- gtf$gene_id
names(key) <- gtf$gene_name
res[,'match'] <- res[]
row.names(res) %in% gtf$gene_name
### Get gene_names from gtf
row.names(gtf) <- gtf$gene_id
res$gtf_name <- gtf[res$gtf_id, 'gene_name']



### Alternatively, defining the results by set1
### Do directy matching from gtf
res <- gene_matching_v1(row.names(sdat), set2, gtf)

### Define a name to be used for plotting
res$plot_name <-NA
res[which(res$gtf_name %in% set2), 'plot_name'] <- res[which(res$gtf_name %in% set2), 'gtf_name']

### Step2, look for aliases from the gtf_name
custom <- res[which(is.na(res$plot_name)), 'gtf_name']
custom <- custom[which(!is.na(custom))]
goal_set <- set2
### Fraction of custom set in the goal_set BEFORE matching
sum(custom %in% goal_set)/length(custom)
### Make dataframe 
mat <- match_sets(custom, goal_set)
### Fraction of custom set in the goal_set AFTER matching
sum(!is.na(mat$match))/length(custom)
mat <- mat[which(!is.na(mat$match)),]
#mat <- mat[which(row.names(mat) == mat$original),]

mat[res[which(is.na(res$plot_name) & (res$gtf_name %in% mat$original)), 'gtf_name'], 'plot_name'] <- 
  mat[res[which(is.na(res$plot_name) & (res$gtf_name %in% mat$original)), 'gtf_name'], 'match']

saveRDS(res, out_file)
