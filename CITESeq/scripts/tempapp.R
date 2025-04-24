### Note: As of 9/2020, the RStudio Connect server requires R version 3.5.0, 3.6.2 or 4.0.2

### Note: this line needs to be run before compiling, so that the compiler has access to the Bioconductor repository for BiocGenerics installation which is not already available on the Connect server.
options(repos = BiocManager::repositories())

library('shiny')
library('DT')
library('ggplot2')
library('dplyr')
library('Seurat')
library('scales')
#library('textshape')
library('grid')
options(bitmapType='cairo')
### Note: This line is necessary because X11 forwarding has some unknown issue causing a plot error

### These are for compiling a shiny
library('RJSONIO')
library('PKI')
library('packrat')
library('rsconnect')
library('XML')
library('MASS')

# source('/data/vrc_his/douek_lab/snakemakes/sc_functions.R')
# load('/data/vrc_his/douek_lab/wakecg/CITESeq/snakemake/results/50Ab/App/Data.RData')
# cseq_ab <- cseq
# load('/data/vrc_his/douek_lab/wakecg/CITESeq/snakemake/results/Biolegend/App/Data.RData')
# cseq_b <- cseq
#load('/data/vrc_his/douek_lab/wakecg/CITESeq/Flu/App/Data.RData')
#cseq0 <- cseq
#load('/data/vrc_his/douek_lab/projects/RNASeq/2021614_21-002/results/Go2/App/Data.RData')
# dim(cseq_ab@assays$prot)
# dim(cseq_b@assays$RNA)
# dim(cseq0@assays$RNA)


load('Data.RData')
source('sc_functions.R')

#source('/data/vrc_his/douek_lab/wakecg/sc_functions.R')
cluster_names = c('RNA_clusters', 'prot_clusters', 'wnn_clusters', 'predicted.celltype.l1', 'predicted.celltype.l2')
names(cluster_names) <- c('RNA', 'Protein', 'WNN', 'Celltype prediction', 'Celltype prediction refined')
cluster_names <- cluster_names[cluster_names %in% colnames(cseq@meta.data)]

wnn_umap <- names(cseq@reductions)[grepl('wnn_umap', names(cseq@reductions))]

my_vln_plot <- function(cseq, de, gene_names, assay = 'prot', groupby = 'prot_clusters', lone_cluster = NA, log = F){
  subtitle <- ''
  ### If DE was done on a lone cluster use that info in subtitles
  if(!is.na(lone_cluster)){
    if(gene_names %in% row.names(de)){
      subtitle <- paste0('DE ', lone_cluster, ' v. others adjp: ', signif(de[gene_names, 'p_val_adj']), ', logFC: ', signif(de[ids, 'avg_log2FC']))
    } 
  }
  if(log){
    vplot <- VlnPlot(cseq, group.by = groupby, features = gene_names, assay = assay) + 
      ggtitle(gene_names, subtitle = subtitle) +  theme(legend.position = "none") + scale_y_log10()
  } else{
    vplot <- VlnPlot(cseq, group.by = groupby, features = gene_names, assay = assay) + 
      ggtitle(gene_names, subtitle = subtitle) +  theme(legend.position = "none")
  }
  return(vplot)
}

ui <- fluidPage(
  titlePanel("Genes by cluster"),
  mainPanel(
    headerPanel('RNA UMAP'),
    selectInput('umap_color1', 'UMAP color', cluster_names, selected = 'RNA_clusters'),
    uiOutput("dynamic_n1"),
    plotOutput('umap_plot1', height = 700)
  ),
  mainPanel(
    headerPanel('Protein UMAP'),
    selectInput('umap_color2', 'UMAP color', cluster_names, selected = 'prot_clusters'),
    uiOutput("dynamic_n2"),
    plotOutput('umap_plot2', height = 700)
  ),
  mainPanel(
    headerPanel('Weighted Nearest Neighbor UMAP'),
    selectInput('umap_color3', 'UMAP color', cluster_names, selected = 'wnn_clusters'),
    uiOutput("dynamic_n3"),
    plotOutput('umap_plot3', height = 700)
  ),
  ### Scatter plots
  mainPanel(
    headerPanel('Scatter Plots'),
    selectInput('data_type1', 'Data type', c('RNA', 'Protein'), selected = 'Protein'),
    uiOutput("searchbar1"),
    uiOutput("searchbar2"),
    selectInput('scatter_group', 'Color by', c(cluster_names, 'Assignment'), selected = 'prot_clusters'),
    plotOutput('scatter_plot', height = 700)
  ),
  
  mainPanel(
    headerPanel('Protein heatmap'),
    selectInput('cluster_type1', 'Cluster type', cluster_names, selected = 'prot_clusters'),
    selectInput('scaled', 'Scaled?', c(F, T), selected = F),
    selectizeInput('prot_search', 'Highlight Protein', choices = row.names(cseq@assays$prot@data), multiple = T),
    plotOutput('heatmap1', height = 700)
  ),
  mainPanel(
    headerPanel('Differential expression of clusters (Seurats FindAllMarkers function)'),
    selectInput('data_type2', 'Data type', c('RNA', 'Protein'), selected = 'RNA'),
    selectInput('cluster_type', 'Cluster type', cluster_names[1:3], selected = 'RNA_clusters'),
    uiOutput("dynamic_n4")
  ),
  mainPanel(
    DT::dataTableOutput("table")
  ),
  mainPanel(
    #selectInput('gene_name', 'Gene Name', fam$gene_name)
    #selectizeInput('gene_name', 'Gene Name', choices = fam$gene_name, multiple = T)
    uiOutput("searchbar")
  ),
  mainPanel(
    plotOutput('violin_plot', height = 700)
  ),
  mainPanel(
    headerPanel('Weighted Nearest Neighbor RNA weights'),
    selectInput('cluster_type2', 'Cluster type', cluster_names, selected = 'RNA_clusters'),
    plotOutput('RNA_weight_plot', height = 700)
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  ### Options to highlight specific clusters on the UMAPs
  output$dynamic_n1 <- renderUI({
    selected_clusters <- unique(cseq@meta.data[, input$umap_color1])
    checkboxGroupInput("highlight1", h3("Highlight cluster?"), choices = selected_clusters[order(selected_clusters)], inline = T) 
  })
  output$dynamic_n2 <- renderUI({
    selected_clusters <- unique(cseq@meta.data[, input$umap_color2])
    checkboxGroupInput("highlight2", h3("Highlight cluster?"), choices = selected_clusters[order(selected_clusters)], inline = T) 
  })
  output$dynamic_n3 <- renderUI({
    selected_clusters <- unique(cseq@meta.data[, input$umap_color3])
    checkboxGroupInput("highlight3", h3("Highlight cluster?"), choices = selected_clusters[order(selected_clusters)], inline = T) 
  })
  output$dynamic_n4 <- renderUI({
    selected_clusters <- unique(cseq@meta.data[, input$cluster_type])
    selected_clusters <- as.list(selected_clusters[order(selected_clusters)])
    checkboxGroupInput("table_cluster", h3("Cluster"), choices = selected_clusters, selected = selected_clusters, inline = T) 
  })
  
  ### UMAP plots
  output$umap_plot1 <- renderPlot({
    pal <- hue_pal()(length(unique(cseq@meta.data[, input$umap_color1])))
    names(pal) <- unique(cseq@meta.data[, input$umap_color1])
    par(mar = c(5.1, 4.1, 0, 0.1))
    if(length(input$highlight1) == 0){
      cells.highlight <- c()
      cols.highlight <- pal[as.character(cseq@meta.data[, input$umap_color1])]
      
      DimPlot(cseq, group.by = input$umap_color1, label = T, reduction = 'RNA_umap', cols = pal)
    } else if(length(input$highlight1) > 0 ){
      cells.highlight <- lapply(as.character(input$highlight1), function(x) colnames(cseq)[which(cseq@meta.data[, input$umap_color1] == x)])
      names(cells.highlight) <- as.character(input$highlight1)
      cols.highlight <- sapply(as.character(input$highlight1), function(x) pal[x])
      names(cols.highlight) <- as.character(input$highlight1)
      
      DimPlot(cseq, group.by = input$umap_color1, label = T, reduction = 'RNA_umap', 
              cells.highlight = cells.highlight,
              cols.highlight = cols.highlight) + 
        theme(legend.position = "none")
    }
  })
  output$umap_plot2 <- renderPlot({
    pal <- hue_pal()(length(unique(cseq@meta.data[, input$umap_color2])))
    names(pal) <- unique(cseq@meta.data[, input$umap_color2])
    par(mar = c(5.1, 4.1, 0, 0.1))
    if(length(input$highlight2) == 0){
      cells.highlight <- c()
      cols.highlight <- pal[as.character(cseq@meta.data[, input$umap_color2])]
      
      DimPlot(cseq, group.by = input$umap_color2, label = T, reduction = 'prot_umap', cols = pal)
    } else if(length(input$highlight2) > 0 ){
      cells.highlight <- lapply(as.character(input$highlight2), function(x) colnames(cseq)[which(cseq@meta.data[, input$umap_color2] == x)])
      names(cells.highlight) <- as.character(input$highlight2)
      cols.highlight <- sapply(as.character(input$highlight2), function(x) pal[x])
      names(cols.highlight) <- as.character(input$highlight2)
      
      DimPlot(cseq, group.by = input$umap_color2, label = T, reduction = 'prot_umap', 
              cells.highlight = cells.highlight,
              cols.highlight = cols.highlight) + 
        theme(legend.position = "none")
    }
  })
  output$umap_plot3 <- renderPlot({
    pal <- hue_pal()(length(unique(cseq@meta.data[, input$umap_color3])))
    names(pal) <- unique(cseq@meta.data[, input$umap_color3])
    par(mar = c(5.1, 4.1, 0, 0.1))
    if(length(input$highlight3) == 0){
      cells.highlight <- c()
      cols.highlight <- pal[as.character(cseq@meta.data[, input$umap_color3])]
      
      DimPlot(cseq, group.by = input$umap_color3, label = T, reduction = wnn_umap, cols = pal)
    } else if(length(input$highlight3) > 0 ){
      cells.highlight <- lapply(as.character(input$highlight3), function(x) colnames(cseq)[which(cseq@meta.data[, input$umap_color3] == x)])
      names(cells.highlight) <- as.character(input$highlight3)
      cols.highlight <- sapply(as.character(input$highlight3), function(x) pal[x])
      names(cols.highlight) <- as.character(input$highlight3)
      
      DimPlot(cseq, group.by = input$umap_color3, label = T, reduction = wnn_umap, 
              cells.highlight = cells.highlight,
              cols.highlight = cols.highlight) + 
        theme(legend.position = "none")
    }
  })
  
  output$searchbar1 <- renderUI({
    assay_names <- c('RNA', 'prot')
    names(assay_names) <- c('RNA', 'Protein')
    assay <- assay_names[input$data_type1]
    features <- row.names(cseq@assays[[assay]]@data)
    selectizeInput('scatter_feature1', 'Gene Name', choices = features, selected = features[1], multiple = F)
  })
  output$searchbar2 <- renderUI({
    assay_names <- c('RNA', 'prot')
    names(assay_names) <- c('RNA', 'Protein')
    assay <- assay_names[input$data_type1]
    features <- row.names(cseq@assays[[assay]]@data)
    selectizeInput('scatter_feature2', 'Gene Name', choices = features, selected = features[2], multiple = F)
  })
  output$scatter_plot <- renderPlot({
    assay_names <- c('RNA', 'prot')
    names(assay_names) <- c('RNA', 'Protein')
    assay <- assay_names[input$data_type1]
    DefaultAssay(cseq) <- assay
    par(mar = c(5.1, 4.1, 0, 0.1))
    FeatureScatter(cseq, feature1 = input$scatter_feature1, feature2 = input$scatter_feature2, 
                   group.by = input$scatter_group)
  })
  
  output$heatmap1 <- renderPlot({
    dsb <- as.matrix(cseq@assays$prot@data)
    if(input$scaled == T){
      dsb <- dsb/rowSums(dsb)
    }
    # calculate the average of each protein separately for each cluster 
    prots = rownames(dsb)
    adt_data <- cbind(cseq@meta.data, as.data.frame(t(dsb)))
    adt_plot <- adt_data %>% group_by(eval(as.name(input$cluster_type1))) 
    adt_plot <- adt_plot[, which(colnames(adt_plot) != input$cluster_type1)]
    colnames(adt_plot)[length(colnames(adt_plot))] <- 'group_var'
    adt_plot <- adt_plot %>% 
      summarize_at(.vars = prots, .funs = mean) %>% 
      column_to_rownames('group_var')
    adt_plot <- t(adt_plot)
    # plot a heatmap of the average dsb normalized values for each cluster
    p <- pheatmap::pheatmap(adt_plot, color = viridis::viridis(25, option = "B"), fontsize_row = 6, border_color = NA)
    adt_plot <- as.data.frame(adt_plot)
    adt_plot$colors <- ifelse(row.names(adt_plot) %in% input$prot_search, 'red', 'black')
    cols = adt_plot[order(match(row.names(adt_plot), p$gtable$grobs[[5]]$label)), ]$colors
    p$gtable$grobs[[5]]$gp=gpar(col=cols, fontsize = 6)
    p
  })
  
  fam <- reactive({
    if(input$data_type2 == 'RNA' & input$cluster_type == cluster_names['RNA']){
      res <- rna_by_rna
    } else if(input$data_type2 == 'RNA' & input$cluster_type == cluster_names['Protein']){
      res <- rna_by_prot
    } else if(input$data_type2 == 'RNA' & input$cluster_type == cluster_names['WNN']){
      res <- rna_by_wnn
    } else if(input$data_type2 == 'Protein' & input$cluster_type == cluster_names['RNA']){
      res <- prot_by_rna
    } else if(input$data_type2 == 'Protein' & input$cluster_type == cluster_names['Protein']){
      res <- prot_by_prot
    } else if(input$data_type2 == 'Protein' & input$cluster_type == cluster_names['WNN']){
      res <- prot_by_wnn
    }
    res
  })
  
  output$table <- DT::renderDataTable(DT::datatable({
    data <- fam()
    data <- data[which(data$cluster %in% input$table_cluster), c('gene', 'cluster', 'p_val', 'p_val_adj', 'avg_log2FC', 'pct.1', 'pct.2', 'diff')]
    #data <- data[, c('gene', 'cluster', 'p_val', 'p_val_adj', 'avg_log2FC', 'pct.1', 'pct.2', 'diff')]
    
    data
  }))
  
  output$searchbar <- renderUI({
    data <- fam()
    selectizeInput('gene_name', 'Gene Name', choices = data$gene, selected = data$gene[1], multiple = T)
  })
  
  output$violin_plot <- renderPlot({
    assay_names <- c('RNA', 'prot')
    names(assay_names) <- c('RNA', 'Protein')
    assay <- assay_names[input$data_type2]
    par(mar = c(5.1, 4.1, 0, 0.1))
    my_vln_plot(cseq, fam(), input$gene_name, assay = assay, groupby = input$cluster_type)
  })
  
  output$RNA_weight_plot <- renderPlot({
    assay_names <- c('RNA', 'prot')
    names(assay_names) <- c('RNA', 'Protein')
    assay <- assay_names[input$data_type2]
    par(mar = c(5.1, 4.1, 0, 0.1))
    my_vln_plot(cseq, fam(), input$gene_name, assay = assay, groupby = input$cluster_type)
    VlnPlot(cseq, features = "RNA.weight", group.by = input$cluster_type2, sort = TRUE, pt.size = 0.1) + NoLegend()
    
  })
}

shinyApp(ui = ui, server = server)
