library('sys')
library('Seurat')
library('SeuratDisk')
library('stringr')
library('pheatmap')
library('ggplot2')
library('umap')
library('textshape')
library('dplyr')
library('biomaRt')
library('grid')
library('scales')
library('celldex')
library('scRNAseq')
library('SingleR')
library('Azimuth')

source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/sc_functions.R')
source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/CITESeq_functions.R')

if(interactive()){
  cseq_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/Go2/Processed_data.RDS'
  exclude_file <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/Go2/Excluded_genes.txt'
  out_rds <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/Go2/Mapped.RDS'
  out_pdf <- '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/Go2/Mapping.pdf'
} else{
  args = commandArgs(trailingOnly=TRUE)
  cseq_file <- args[1]
  out_rds <- args[2]
  out_pdf <- args[3]
}

print('Reading data')
cseq <- readRDS(cseq_file)
DefaultAssay(cseq) <- 'RNA'

#make_test_rds(cseq_file)
#cseq <- read_test_rds(cseq_file)
### Normalize with same method used on the query data
print('SC transformation')
cseq <- SCTransform(cseq, assay = 'RNA', verbose = FALSE)
print(names(cseq@assays))
print('reading reference data')
reference <- LoadH5Seurat("/hpcdata/vrc/vrc1_data/douek_lab/reference_sets/pbmc/pbmc_multimodal.h5seurat")
print(reference)
print('Finding anchors')
anchors <- FindTransferAnchors(
  reference = reference,
  query = cseq,
  normalization.method = "SCT",
  reference.reduction = 'spca',
  dims = 1:50
)

### prediction.score.celltype.l1,prediction.score.celltype.l2 and predicted_ADT become assays, containing a feature per reference cell type.
### The values per cell are in the data slot, and are prediction scores ranging from 0 to 1.
print('mapping query to reference')
cseq <- MapQuery(
  anchorset = anchors,
  query = cseq,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

### predicted.celltype.l1 and 2 contain the predicted cell type from the reference
### predicted.celltype.l1.score and 2.score contain the ... maximum score???
### 
print(table(cseq$predicted.celltype.l2))
Idents(cseq) <- 'predicted.celltype.l2'
#treg_markers <- FindMarkers(cseq, ident.1 = "Treg", only.pos = TRUE, logfc.threshold = 0.1)
#print(head(treg_markers))
DefaultAssay(cseq) <- 'RNA'
saveRDS(cseq, out_rds)

pdf(out_pdf)
DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
### ref.umap is defined by the reference
DimPlot(cseq, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
DimPlot(cseq, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()

### Cell type scores, to compare to the cell type declarations
FeaturePlot(cseq, features = c("pDC", "CD16 Mono", "Treg"),  reduction = "ref.umap", cols = c("lightgrey", "darkred"), ncol = 3) & theme(plot.title = element_text(size = 10))

VlnPlot(cseq, features = c("CLEC4C", "LILRA4"), sort = TRUE) + NoLegend()

DefaultAssay(cseq) <- 'predicted_ADT'
# see a list of proteins: rownames(cseq)
FeaturePlot(cseq, features = c("CD3-1", "CD45RA", "IgD"), reduction = "ref.umap", cols = c("lightgrey", "darkgreen"), ncol = 3)


### "De novo" visualization, to catch if there are cell types in query that are not within reference
#merge reference and query
reference$id <- 'reference'
cseq$id <- 'query'
refquery <- merge(reference, cseq)
refquery[["spca"]] <- merge(reference[["spca"]], cseq[["ref.spca"]])
refquery <- RunUMAP(refquery, reduction = 'spca', dims = 1:50)
DimPlot(refquery, group.by = 'id', shuffle = TRUE)

dev.off()

### To create a test set for interactive use
# csub <- subset(x = cseq, downsample = 10000)
# saveRDS(csub, file = '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/Go2/Test_subset.RDS')

### To create a test set for interactive use
# csub <- subset(x = cseq, downsample = 1000)
# print(dim(csub))
# saveRDS(csub, file = '/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/Go2/Test_suberset.RDS')

# cols <- c('MT_sum', 'orig.ident', 'propmt', 'rna_size', 'ngene', 'bc')
# for(c in cols){
#   cseq[[c]] <- NULL
# }
