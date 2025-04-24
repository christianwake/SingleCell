

if(interactive()){
  in_file <- '/data/vrc_his/douek_lab/projects/RNASeq/2021600_kristin/data/bam/Nreads.csv'
  out_file <- '/data/vrc_his/douek_lab/projects/RNASeq/2021600_kristin/data/Post_filter.csv'
  
  in_file <- '/data/vrc_his/douek_lab/projects/RNASeq/2022620_857.1/data/bam/Nreads.csv'
  out_file <- '/data/vrc_his/douek_lab/projects/RNASeq/2022620_857.1/data/Post_filter.csv'
  
} else{
  args = commandArgs(trailingOnly=TRUE)
  in_file <- args[1]
  #filt <- args[2]
  out_file <- args[2] 
}

### Minimum number of reads (column N)
filt <- 0
covs <- read.csv(in_file, stringsAsFactors = F, header = F, check.names = F)
colnames(covs) <- c('ID', 'N')
row.names(covs) <- covs$ID

covs <- covs[which(covs$N > filt),]
write.table(covs, out_file, sep = ',', quote = F, row.names = F, col.names = T)



