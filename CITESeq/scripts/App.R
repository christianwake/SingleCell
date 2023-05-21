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

#source('/hpcdata/vrc/vrc1_data/douek_lab/snakemakes/sc_functions.R')
#load('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/snakemake/results/50Ab/App/Data.RData')
#load('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/snakemake/results/Biolegend/App/Data.RData')
#load('/hpcdata/vrc/vrc1_data/douek_lab/wakecg/CITESeq/Flu/App/Data.RData')
#load('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/2021614_21-002/results/Go2/App/Data.RData')

load('Data.RData')
source('sc_functions.R')

### Kristin-specific
data_types <- c('RNA_umap')
names(data_types) <- c('RNA')
cluster_names <- c('seurat_clusters')
names(cluster_names) <- c('RNA')

# data_type_patterns <- list(c('RNA_umap'), c('prot_umap'), c('wnn_umap'))
# names(data_type_patterns) <- c('RNA', 'Protein' , 'Weighted Nearest Neighbor')
# ### Grepl possible patterns for each possible data type
# data_types <- lapply(data_type_patterns, function(dt) ### For each possible data type
#   sapply(dt, function(dt_name) ### Loop over its associated grepl patterns
#     names(sdat@reductions)[sapply(names(sdat@reductions), function(x) grepl(dt_name, x))]
# data_types <- data_types[which(data_types != 'character(0)')]
# ### Old, specific to citseq
# #names(data_types) <- c('RNA', 'Protein')
# #data_types <- c('RNA_umap', 'prot_umap', names(sdat@reductions)[grepl('wnn_umap', names(sdat@reductions))])
# #names(data_types) <- c('RNA', 'Protein', 'Weighted Nearest Neighbor')
# 
# cluster_name_patterns <- list(c('seurat_clusters', 'RNA_clusters'), c('prot_clusters'), c('wnn_clusters'), c('predicted.celltype.l1'), c('predicted.celltype.l2'))
# names(cluster_name_patterns) <- c('RNA', 'Protein', 'WNN', 'Celltype prediction', 'Celltype prediction refined')
# 
# cluster_names <- sapply(names(cluster_name_patterns), function(cn) cluster_name_patterns[[cn]][sapply(cluster_name_patterns[[cn]], function(cn_name) any(grepl(cn_name, colnames(sdat@meta.data))))])
# names(cluster_names) <- names(cluster_name_patterns)
# cluster_names <- cluster_names[which(cluster_names != 'character(0)')]
# 
# ### Old, specific to citseq
# # cluster_names = c('RNA_clusters', 'prot_clusters', 'wnn_clusters', 'predicted.celltype.l1', 'predicted.celltype.l2')
# # names(cluster_names) <- c('RNA', 'Protein', 'WNN', 'Celltype prediction', 'Celltype prediction refined')
# cluster_names <- cluster_names[cluster_names %in% colnames(sdat@meta.data)]

ui <- fluidPage(
  titlePanel("Genes by cluster"),
  mainPanel(
    headerPanel('UMAP'),
    selectInput('umap_data', 'UMAP data', names(data_types), selected = 'RNA'), # CS-specific
    selectInput('umap_color', 'UMAP color', cluster_names, selected = 'RNA_clusters'),
    uiOutput("dynamic_n1"),
    plotOutput('umap_plot', height = 700),
    DT::dataTableOutput("freq_table"), 
    conditionalPanel( # Sub-clustering specific
      #mainPanel(
      #condition = "input.highlight == 0",
      condition = "Object.keys(input.highlight).length",
      headerPanel('sub UMAP based on above highlighted clusters'),
      uiOutput("dynamic_umap_subdata"),
      uiOutput("dynamic_umap_subcolor"),
      #selectInput('umap_subdata', 'sub-UMAP data', names(data_types), selected = 'RNA),
      #selectInput('umap_subcolor', 'sub-UMAP color', cluster_names, selected = 'RNA_clusters'),
      plotOutput('umap_subplot', height = 700),
      ### Frequency table based on selected clusters (from above)
      DT::dataTableOutput("sub_freq_table")
    )
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
    selectizeInput('prot_search', 'Highlight Protein', choices = row.names(sdat@assays$prot@data), multiple = T),
    plotOutput('heatmap1', height = 700)
  ),
  mainPanel(
    headerPanel('Differential expression of clusters (Seurats FindAllMarkers function)'),
    selectInput('data_type2', 'Data type', c('RNA', 'Protein'), selected = 'RNA'),
    selectInput('cluster_type', 'Cluster type', cluster_names[1:3], selected = 'RNA_clusters'),
    uiOutput("dynamic_n2")
  ),
  mainPanel(
    DT::dataTableOutput("DE_table")
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
    selected_clusters <- unique(sdat@meta.data[, input$umap_color])
    checkboxGroupInput("highlight", h3("Highlight cluster?"), choices = selected_clusters[order(selected_clusters)], inline = T) 
  })
  output$dynamic_n2 <- renderUI({
    selected_clusters <- unique(sdat@meta.data[, input$cluster_type])
    selected_clusters <- as.list(selected_clusters[order(selected_clusters)])
    checkboxGroupInput("table_cluster", h3("Cluster"), choices = selected_clusters, selected = selected_clusters, inline = T) 
  })
  ### dynamic_umap_subdata is the name refered to by uiOutput (above)
  output$dynamic_umap_subdata <- renderUI({
    ### Add dimensionality reduction options based on the highlighted subset of cells
    data_options <- c(data_types, paste0('Subset_', input$umap_color, '_', input$highlight, '_', c('RNA', 'prot')))
    names(data_options) <- c(names(data_types), paste0('Subset_', input$umap_color, '_', input$highlight, '_', c('RNA', 'prot')))
    ### But remove them if they aren't in the Seurat object
    data_options <- data_options[data_options %in% names(sdat@reductions)]
    ###umap_subdata is the name refered to in input$
    selectInput('umap_subdata', 'sub-UMAP data', names(data_options), selected = input$umap_data)
  })
  output$dynamic_umap_subcolor <- renderUI({
    ### Add clusters options based on the highlighted subset of cells
    cluster_options <- c(cluster_names, paste0('Subset_', input$umap_color, '_', input$highlight, '_', c('RNA', 'prot'), '_clusters'))
    ### But remove them if they aren't in the Seurat object
    cluster_options <- cluster_options[cluster_options %in% colnames(sdat@meta.data)]
    selectInput('umap_subcolor', 'sub-UMAP color', cluster_options, selected = input$umap_color)
  })
  ### UMAP plots
  output$umap_plot <- renderPlot({
    pal <- hue_pal()(length(unique(sdat@meta.data[, input$umap_color])))
    names(pal) <- unique(sdat@meta.data[, input$umap_color])
    par(mar = c(5.1, 4.1, 0, 0.1))
    if(length(input$highlight) == 0){
      cells.highlight <- c()
      cols.highlight <- pal[as.character(sdat@meta.data[, input$umap_color])]
      
      DimPlot(sdat, group.by = input$umap_color, label = T, reduction = data_types[input$umap_data], cols = pal)
    } else if(length(input$highlight) > 0 ){
      cells.highlight <- lapply(as.character(input$highlight), function(x) colnames(sdat)[which(sdat@meta.data[, input$umap_color] == x)])
      names(cells.highlight) <- as.character(input$highlight)
      cols.highlight <- sapply(as.character(input$highlight), function(x) pal[x])
      names(cols.highlight) <- as.character(input$highlight)
      
      DimPlot(sdat, group.by = input$umap_color, label = T, reduction = data_types[input$umap_data], 
              cells.highlight = cells.highlight,
              cols.highlight = cols.highlight) + 
        theme(legend.position = "none")
    }
  })
  
  output$freq_table <- DT::renderDataTable(DT::datatable({
    data <- as.data.frame(table(sdat@meta.data[, input$umap_color]))
    data[,3] <- data$Freq/sum(data$Freq)
    colnames(data) <- c(input$umap_color, 'N', 'Fraction')
    data
  }))
  
  ### Optional UMAP of the highlighted subset of cells
  output$umap_subplot <- renderPlot({
    ### Use the value associated with the input name, or if there is no associated value, use the name itself as the value
    reduc <- data_types[input$umap_subdata]
    if(is.na(reduc)){
      reduc <- input$umap_subdata
    }
    cells <- colnames(sdat)[sdat@meta.data[, input$umap_color] %in% c(input$highlight)]
    sdatsub <- subset(sdat, cells = cells)
    pal <- hue_pal()(length(unique(sdatsub@meta.data[, input$umap_subcolor])))
    names(pal) <- unique(sdatsub@meta.data[, input$umap_subcolor])
    par(mar = c(5.1, 4.1, 0, 0.1))

    cells.highlight <- c()
    cols.highlight <- pal[as.character(sdatsub@meta.data[, input$umap_subcolor])]
      
    DimPlot(sdatsub, group.by = input$umap_subcolor, label = T, reduction = reduc, cols = pal)
  })
  output$sub_freq_table <- DT::renderDataTable(DT::datatable({
    data <- as.data.frame(table(sdat@meta.data[, input$umap_subcolor]))
    data[,3] <- data$Freq/sum(data$Freq)
    colnames(data) <- c(input$umap_subcolor, 'N', 'Fraction')
    data
  }))
  ### Scatterplot inputs
  output$searchbar1 <- renderUI({
    assay_names <- c('RNA', 'prot')
    names(assay_names) <- c('RNA', 'Protein')
    assay <- assay_names[input$data_type1]
    features <- row.names(sdat@assays[[assay]]@data)
    selectizeInput('scatter_feature1', 'Gene Name', choices = features, selected = features[1], multiple = F)
  })
  output$searchbar2 <- renderUI({
    assay_names <- c('RNA', 'prot')
    names(assay_names) <- c('RNA', 'Protein')
    assay <- assay_names[input$data_type1]
    features <- row.names(sdat@assays[[assay]]@data)
    selectizeInput('scatter_feature2', 'Gene Name', choices = features, selected = features[2], multiple = F)
  })
  output$scatter_plot <- renderPlot({
    assay_names <- c('RNA', 'prot')
    names(assay_names) <- c('RNA', 'Protein')
    assay <- assay_names[input$data_type1]
    DefaultAssay(sdat) <- assay
    par(mar = c(5.1, 4.1, 0, 0.1))
    FeatureScatter(sdat, feature1 = input$scatter_feature1, feature2 = input$scatter_feature2, 
                   group.by = input$scatter_group)
  })
  
  output$heatmap1 <- renderPlot({
    dsb <- as.matrix(sdat@assays$prot@data)
    if(input$scaled == T){
      dsb <- dsb/rowSums(dsb)
    }
    # calculate the average of each protein separately for each cluster 
    prots = rownames(dsb)
    adt_data <- cbind(sdat@meta.data, as.data.frame(t(dsb)))
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
  ### Find all Markers DE data
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
  
  output$DE_table <- DT::renderDataTable(DT::datatable({
    data <- fam()
    data <- data[which(data$cluster %in% input$table_cluster), c('gene', 'cluster', 'p_val', 'p_val_adj', 'avg_log2FC', 'pct.1', 'pct.2', 'diff')]
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
    my_vln_plot(sdat, fam(), input$gene_name, assay = assay, groupby = input$cluster_type)
  })
  
  output$RNA_weight_plot <- renderPlot({
    assay_names <- c('RNA', 'prot')
    names(assay_names) <- c('RNA', 'Protein')
    assay <- assay_names[input$data_type2]
    par(mar = c(5.1, 4.1, 0, 0.1))
    my_vln_plot(sdat, fam(), input$gene_name, assay = assay, groupby = input$cluster_type)
    VlnPlot(sdat, features = "RNA.weight", group.by = input$cluster_type2, sort = TRUE, pt.size = 0.1) + NoLegend()
  })
}

shinyApp(ui = ui, server = server)
