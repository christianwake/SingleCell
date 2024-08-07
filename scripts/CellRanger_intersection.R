library('sys')

if(interactive()){
  project <- '2021614_21-002'
  qc_name <- '2024-01-20'
  method <- 'Custom'
  calls_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, 
                       '/data/Dehash_calls_', method, '.tsv')
  out_file <- paste0('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/', project, 
                     '/results/', qc_name, '/Dehash_calls_', method, '_CRint.tsv')
}else{
  args <- commandArgs(trailingOnly=TRUE)
  calls_file <- args[1]
  method <- args[2]
  out_file <- args[3]
}

if(!file.exists(calls_file)){
  calls_file <- paste0('data/Dehash_calls_', method, 'tsv')
}

calls <- read.table(calls_file, sep = '\t', header = T, quote = '')
calls[, 'Assignment'] <- ifelse(calls$Assignment_CR == 'Background', 'Negative', calls$Assignment)
calls[, 'Assignment_simple'] <- ifelse(calls$Assignment_CR == 'Background', 'Negative', calls$Assignment_simple)

write.table(calls, out_file, sep = '\t', quote = F, row.names = F, col.names = T)
