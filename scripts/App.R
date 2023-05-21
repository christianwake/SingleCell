### Note: As of 9/2020, the RStudio Connect server requires R version 3.5.0, 3.6.2 or 4.0.2

### Note: this line needs to be run before compiling, so that the compiler has access to the Bioconductor repository for BiocGenerics installation which is not already available on the Connect server.
options(repos = BiocManager::repositories())

library('shiny')
library('DT')
library('ggplot2')
library('dplyr')
library('Seurat')
library('scales')
library('textshape')
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

# load('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021617_mis-c/results/Dropout_mitigated/App/Data.RData')
# load('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021618_galt/results/Go1/App/Data.RData')

load('Data.RData')
#source('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/sc_functions.R')
source('sc_functions.R')

### Find all cluster options from the serurat object and choose from our list, the first in our order of preference
cluster_names <- c('RNA_clusters', 'clusters_seurat_integrated', 'clusters_harmony_integrated', 'cluster_unintegrated', 'seurat_clusters', 'predicted.celltype.l1', 'predicted.celltype.l2', 'predicted_label_main', 'predicted_label_fine')
#cluster_names <- cluster_names[cluster_names %in% colnames(sdat@meta.data)[grepl('clusters', colnames(sdat@meta.data))]]
cluster_names <- cluster_names[cluster_names %in% colnames(sdat@meta.data)]
extra_annots <- c('Assignment', 'Condition', 'CR_ID')
extra_annots <- extra_annots[extra_annots %in% colnames(sdat@meta.data)]

umap_name <- c('umap_seurat', 'umap_harmony', 'umap')
umap_name <- umap_name[umap_name %in% names(sdat@reductions)[grepl('umap', names(sdat@reductions))]]
umap_name <- umap_name[1]

ui <- fluidPage(
  titlePanel("Genes by cluster"),
  mainPanel(
    headerPanel('RNA UMAP'),
    selectInput('umap_color1', 'UMAP color', cluster_names, selected = cluster_names[1]),
    uiOutput("dynamic_n1"),
    plotOutput('umap_plot1', height = 700)
  ),
  ### Scatter plots
  mainPanel(
    headerPanel('Scatter Plots'),
    uiOutput("searchbar1"),
    uiOutput("searchbar2"),
    selectInput('scatter_group', 'Color by', c(cluster_names, extra_annots), selected = cluster_names[1]),
    plotOutput('scatter_plot', height = 700)
  ),
  mainPanel(
    headerPanel('Differential expression of clusters (Seurats FindAllMarkers function)'),
    selectInput('cluster_type', 'Cluster type', cluster_names, selected = cluster_names[1]),
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
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  ### Options to highlight specific clusters on the UMAPs
  output$dynamic_n1 <- renderUI({
    selected_clusters <- unique(sdat@meta.data[, input$umap_color1])
    checkboxGroupInput("highlight1", h3("Highlight cluster?"), choices = selected_clusters[order(selected_clusters)], inline = T) 
  })
  output$dynamic_n4 <- renderUI({
    selected_clusters <- unique(sdat@meta.data[, input$cluster_type])
    selected_clusters <- as.list(selected_clusters[order(selected_clusters)])
    checkboxGroupInput("table_cluster", h3("Cluster"), choices = selected_clusters, selected = selected_clusters, inline = T) 
  })
  ### UMAP plots
  output$umap_plot1 <- renderPlot({
    pal <- hue_pal()(length(unique(sdat@meta.data[, input$umap_color1])))
    names(pal) <- unique(sdat@meta.data[, input$umap_color1])
    par(mar = c(5.1, 4.1, 0, 0.1))
    if(length(input$highlight1) == 0){
      cells.highlight <- c()
      cols.highlight <- pal[as.character(sdat@meta.data[, input$umap_color1])]
      
      DimPlot(sdat, group.by = input$umap_color1, label = T, reduction = umap_name, cols = pal)
    } else if(length(input$highlight1) > 0 ){
      cells.highlight <- lapply(as.character(input$highlight1), function(x) colnames(sdat)[which(sdat@meta.data[, input$umap_color1] == x)])
      names(cells.highlight) <- as.character(input$highlight1)
      cols.highlight <- sapply(as.character(input$highlight1), function(x) pal[x])
      names(cols.highlight) <- as.character(input$highlight1)
      
      DimPlot(sdat, group.by = input$umap_color1, label = T, reduction = umap_name, 
              cells.highlight = cells.highlight,
              cols.highlight = cols.highlight) + 
        theme(legend.position = "none")
    }
  })
  
  output$searchbar1 <- renderUI({
    features <- row.names(sdat@assays[['RNA']]@data)
    selectizeInput('scatter_feature1', 'Gene Name', choices = features, selected = features[1], multiple = F)
  })
  output$searchbar2 <- renderUI({

    features <- row.names(sdat@assays[['RNA']]@data)
    selectizeInput('scatter_feature2', 'Gene Name', choices = features, selected = features[2], multiple = F)
  })
  output$scatter_plot <- renderPlot({
    DefaultAssay(sdat) <- 'RNA'
    par(mar = c(5.1, 4.1, 0, 0.1))
    FeatureScatter(sdat, feature1 = input$scatter_feature1, feature2 = input$scatter_feature2, 
                   group.by = input$scatter_group)
  })

  output$table <- DT::renderDataTable(DT::datatable({
    data <- rna_by_rna
    data <- data[which(data$cluster %in% input$table_cluster), c('gene', 'cluster', 'p_val', 'p_val_adj', 'avg_log2FC', 'pct.1', 'pct.2', 'diff')]
    #data <- data[, c('gene', 'cluster', 'p_val', 'p_val_adj', 'avg_log2FC', 'pct.1', 'pct.2', 'diff')]
    
    data
  }))
  
  output$searchbar <- renderUI({
    data <- rna_by_rna
    selectizeInput('gene_name', 'Gene Name', choices = data$gene, selected = data$gene[1], multiple = T)
  })
  
  output$violin_plot <- renderPlot({
    par(mar = c(5.1, 4.1, 0, 0.1))
    my_vln_plot(sdat, fam(), input$gene_name, assay = 'RNA', groupby = input$umap_color1)
  })
  

}

shinyApp(ui = ui, server = server)
