require(tidyverse)
library(shiny)
require(Seurat)
require(patchwork)
library(BiocManager)
library(AnnotationDbi)
library(DESeq2)
library(org.Hs.eg.db)
library(pathview)
library(gage)
library(gageData)
require(msigdbr)
require(fgsea)
library(ComplexHeatmap)
options(repos = BiocManager::repositories())
require(DT)
require(shinydashboard)
require(shinyFiles)
require(reshape2)
require(Matrix)
require(shinyjs)
require(RColorBrewer)
require(igraph)
require(shinyDarkmode)
require(shinyWidgets)
require(shinyBS)
options(shiny.maxRequestSize=500*1024^2) 
options(warn=-1)
graphics.off()

shinyServer(function(input, output, session){
  
  #################### Hide toggle from sidebar ############
  session$sendCustomMessage(type = "manipulateMenuItem1", message = list(action = "hide",toggle = "offcanvas", role = "button"))
  ######################################################
  
  
  ###################### Dark mode ####################
  darkmode_toggle(inputid = 'togglemode')
  #####################################################
  
  ####################################### single cell ######################################
  
  #################################### error checking ######################################
  output$files_sanity_checking <- reactive({
    error_catching <- function(){
      mat_n <- grep('^matrix.mtx', input$mydata[3:6])
      bar_n <- grep('^barcodes.tsv', input$mydata[3:6])
      gene_n <- grep('^genes.tsv', input$mydata[3:6])
      num_files <- length(mat_n) + length(bar_n) + length(gene_n)
      if( length(input$mydata) != 6 | num_files != 3){
        FALSE
      } else {
        TRUE
      }
    }
    error_catching()
  })
  outputOptions(output, 'files_sanity_checking', suspendWhenHidden=FALSE)
  
  reading_10x <- reactive({
      files <- input$fileIn
      files <- files[order(files[,1]), ]
      # path <- files$datapath[1]
      # path <- str_extract(path, "/tmp.*(?=[0-2]{1}.)")
      # setwd(path)
      # file.rename(list.files(pattern="0.gz|0.tsv"), "barcodes.tsv")
      # file.rename(list.files(pattern="1.gz|1.tsv"), "genes.tsv")
      # file.rename(list.files(pattern="2.gz|2.mtx"), "matrix.tsv")
      # files$datapath[1]
      ReadMtx_enhanced <- function(
        mtx,
        cells,
        features,
        cell.column = 1,
        feature.column = 2,
        cell.sep = "\t",
        feature.sep = "\t",
        skip.cell = 0,
        skip.feature = 0,
        mtx.transpose = FALSE,
        unique.features = TRUE,
        strip.suffix = FALSE
      ) {
        all.files <- list(
          "expression matrix" = mtx,
          "barcode list" = cells,
          "feature list" = features
        )
        cell.barcodes <- read.table(
          file = all.files[['barcode list']],
          header = FALSE,
          sep = cell.sep,
          row.names = NULL,
          skip = skip.cell
        )
        feature.names <- read.table(
          file = all.files[['feature list']],
          header = FALSE,
          sep = feature.sep,
          row.names = NULL,
          skip = skip.feature
        )
        # read barcodes
        bcols <- ncol(x = cell.barcodes)
        if (bcols < cell.column) {
          stop(
            "cell.column was set to ",
            cell.column,
            " but ",
            cells,
            " only has ",
            bcols,
            " columns.",
            " Try setting the cell.column argument to a value <= to ",
            bcols,
            "."
          )
        }
        cell.names <- cell.barcodes[, cell.column]
        if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
          cell.names <- as.vector(x = as.character(x = sapply(
            X = cell.names,
            FUN = ExtractField,
            field = 1,
            delim = "-"
          )))
        }
        # read features
        fcols <- ncol(x = feature.names)
        if (fcols < feature.column) {
          stop(
            "feature.column was set to ",
            feature.column,
            " but ",
            features,
            " only has ",
            fcols, " column(s).",
            " Try setting the feature.column argument to a value <= to ",
            fcols,
            "."
          )
        }
        if (any(is.na(x = feature.names[, feature.column]))) {
          na.features <- which(x = is.na(x = feature.names[, feature.column]))
          replacement.column <- ifelse(test = feature.column == 2, yes = 1, no = 2)
          if (replacement.column > fcols) {
            stop(
              "Some features names are NA in column ",
              feature.column,
              ". Try specifiying a different column.",
              call. = FALSE
            )
          } else {
            warning(
              "Some features names are NA in column ",
              feature.column,
              ". Replacing NA names with ID from column ",
              replacement.column,
              ".",
              call. = FALSE
            )
          }
          feature.names[na.features, feature.column] <- feature.names[na.features, replacement.column]
        }
        feature.names <- feature.names[, feature.column]
        if (unique.features) {
          feature.names <- make.unique(names = feature.names)
        }
        data <- readMM(file = all.files[['expression matrix']])
        if (mtx.transpose) {
          data <- t(x = data)
        }
        if (length(x = cell.names) != ncol(x = data)) {
          stop(
            "Matrix has ",
            ncol(data),
            " columns but found ", length(cell.names),
            " barcodes. ",
            ifelse(
              test = length(x = cell.names) > ncol(x = data),
              yes = "Try increasing `skip.cell`. ",
              no = ""
            ),
            call. = FALSE
          )
        }
        if (length(x = feature.names) != nrow(x = data)) {
          stop(
            "Matrix has ",
            nrow(data),
            " rows but found ", length(feature.names),
            " features. ",
            ifelse(
              test = length(x = feature.names) > nrow(x = data),
              yes = "Try increasing `skip.feature`. ",
              no = ""
            ),
            call. = FALSE
          )
        }

        colnames(x = data) <- cell.names
        rownames(x = data) <- feature.names
        data <- as(data, Class = "CsparseMatrix")
        return(data)
      }

      myturn <- ReadMtx_enhanced(mtx = files$datapath[3],
                                 cells = files$datapath[1],
                                 features = files$datapath[2],
                                 cell.column = input$barcodescol,
                                 feature.column =input$featurescol)

      proj_data <- CreateSeuratObject(counts = myturn,
                                      project = input$projectname,
                                      min.cells = 3,
                                      min.features = 200)

  })

  #Table of meta data
  output$metadata <- renderDT({
    if(input$start_seurat){
    req(reading_10x())
      data <- isolate(reading_10x())
      data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
      meta <- data@meta.data
      meta$percent.mt <- round(meta$percent.mt, 2)
      meta %>%
        datatable(
          fillContainer = F,
          options = list(scrollX = TRUE, pageLength = 5),
          extensions = "AutoFill",
          style = "auto",
          rownames = TRUE
        )
    }
  })

#   # QC plots
  output$qcplot0 <- renderPlot({
    if(input$qcplotaction){
      req(reading_10x())
      data <- isolate(reading_10x())
      data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
      VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, cols = "#873600", pt.size = 1.1)
    }
  })
#   ###############################################
  output$qcplot1 <- renderPlot({
    if(input$qcplotaction){
      req(reading_10x())
      meta <- isolate(reading_10x()@meta.data)
      ggplot(meta, aes(x= nCount_RNA, color = orig.ident, fill = orig.ident)) +
        geom_density(alpha= 0.7, color='#FAE5D3',fill='#7E5109') +
        scale_x_log10() +
        ylab("Cell density") +
        geom_vline(xintercept = 500, linetype='dotted', col = 'red') +
        theme_classic() +
        theme(axis.text = element_text(family = "Times",size = 13 , colour = "black"),
              axis.text.x = element_text(family = "Times",colour = "black", size = 13),
              axis.text.y = element_text(family = "Times",colour = "black"),
              plot.subtitle = element_text(family = "Times",size = 20, colour = "black", hjust = 0.5),
              axis.title.y = element_text(family = "Times", size = rel(1.4), angle = 90),
              axis.title.x = element_text(family = "Times", size = rel(1.4), angle = 00))
    }
  })
#   #########################################
  output$qcplot2 <- renderPlot({
    if(input$qcplotaction){
      req(reading_10x())
      meta <- isolate(reading_10x()@meta.data)
      ggplot(meta, aes(x= nFeature_RNA, color = orig.ident, fill = orig.ident)) +
        geom_density(alpha= 0.7, color='#FAE5D3',fill='#7E5109') +
        scale_x_log10() +
        ylab("Cell density") +
        geom_vline(xintercept = 500, linetype='dotted', col = 'red') +
        theme_classic() +
        theme(axis.text = element_text(family = "Times",size = 13 , colour = "black"),
              axis.text.x = element_text(family = "Times",colour = "black", size = 13),
              axis.text.y = element_text(family = "Times",colour = "black"),
              plot.subtitle = element_text(family = "Times",size = 20, colour = "black", hjust = 0.5),
              axis.title.y = element_text(family = "Times", size = rel(1.4), angle = 90),
              axis.title.x = element_text(family = "Times", size = rel(1.4), angle = 00))
    }
  })
#   ##########################################
  output$qcplot3 <- renderPlot({
    if(input$qcplotaction){
      req(reading_10x())
      data <- isolate(reading_10x())
      data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
      meta <- data@meta.data
      ggplot(meta, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
        geom_point() +
        scale_color_gradient(low = "#D7DBDD", high = "#186A3B") +
        geom_smooth(method = lm, color = "#F8C471") +
        scale_x_log10() +
        scale_y_log10() +
        ylab("Cell density") +
        theme_classic() +
        theme(axis.text = element_text(family = "Times",size = 13 , colour = "black"),
              axis.text.x = element_text(family = "Times",colour = "black", size = 13),
              axis.text.y = element_text(family = "Times",colour = "black"),
              plot.subtitle = element_text(family = "Times",size = 20, colour = "black", hjust = 0.5),
              axis.title.y = element_text(family = "Times", size = rel(1.4), angle = 90),
              axis.title.x = element_text(family = "Times", size = rel(1.4), angle = 00))

    }
  })
#   #########################################
  output$qcplot4 <- renderPlot({
    if(input$qcplotaction){
      req(reading_10x())
      data <- isolate(reading_10x())
      data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
      FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = "#D68910", pt.size = 2)
    }
  })
#   ##########################################
  output$qcplot5 <- renderPlot({
    if(input$qcplotaction){
      req(reading_10x())
      data <- isolate(reading_10x())
      FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = "#D68910", pt.size = 2)
    }
  })
#   ###########################################
  output$qcplot6 <- renderPlot({
    if(input$qcplotaction){
      req(reading_10x())
      data <- isolate(reading_10x())
      data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
      mto <- data.frame(mt = sort(data@meta.data$percent.mt))
      ggplot(mto, aes(x = 1:nrow(mto), y=mt)) +
        geom_point() +
        xlab("Cells sorted by percentage mitochondrial counts") +
        ylab("Percentage of mitochondrial counts") +
        geom_hline(yintercept = 5, linetype='dotted', col = 'red') +
        geom_hline(yintercept = 10, linetype='dotted', col = 'red') +
        geom_hline(yintercept = 15, linetype='dotted', col = 'red') +
        theme_classic() +
        theme(axis.text = element_text(family = "Times",size = 13 , colour = "black"),
              axis.text.x = element_text(family = "Times",colour = "black", size = 13),
              axis.text.y = element_text(family = "Times",colour = "black"),
              plot.subtitle = element_text(family = "Times",size = 20, colour = "black", hjust = 0.5),
              axis.title.y = element_text(family = "Times", size = rel(1.4), angle = 90),
              axis.title.x = element_text(family = "Times", size = rel(1.4), angle = 00))
    }
  })
#   ################################# Filtering ####################################
  filtering <- reactive({
    if(input$filteringaction){
      req(reading_10x())
      data <- isolate(reading_10x())
      data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
      data <- subset(data,
                     subset = nFeature_RNA > input$filteringlessfeature &
                     nFeature_RNA < input$filteringmorefeature &
                     percent.mt < input$filteringmtgenes)
      counts <- GetAssayData(object = data, slot = "counts")
      nonzero <- counts > 0
      keep_genes <- Matrix::rowSums(nonzero) >= input$filteringgeneswithzero
      filtered_counts <- counts[keep_genes, ]
      data <- CreateSeuratObject(filtered_counts, meta.data = data@meta.data)
    }
  })
  
  
  filteringresult_reactive <- reactive({
    if(input$filteringaction){
      req(reading_10x())
      req(filtering())
      data_b <- isolate(reading_10x())
      data_a <- isolate(filtering())
      # before
      Before <- data.frame(
        Total_cells = nrow(data_b@meta.data),
        Mean_n_genes = mean(data_b@meta.data$nFeature_RNA),
        SD_n_genes = sd(data_b@meta.data$nFeature_RNA),
        Max_n_genes = max(data_b@meta.data$nFeature_RNA),
        Min_n_genes = min(data_b@meta.data$nFeature_RNA),
        Mean_UMI = mean(data_b@meta.data$nCount_RNA),
        SD_UMI = sd(data_b@meta.data$nCount_RNA),
        Max_UMI = max(data_b@meta.data$nCount_RNA),
        Min_UMI = min(data_b@meta.data$nCount_RNA),
        Total_genes = nrow(data_b)) %>%
        mutate_all(function(x) round(x, 2))
      ##### after
      After <- data.frame(
        Total_cells = nrow(data_a@meta.data),
        Mean_n_genes = mean(data_a@meta.data$nFeature_RNA),
        SD_n_genes = sd(data_a@meta.data$nFeature_RNA),
        Max_n_genes = max(data_a@meta.data$nFeature_RNA),
        Min_n_genes = min(data_a@meta.data$nFeature_RNA),
        Mean_UMI = mean(data_a@meta.data$nCount_RNA),
        SD_UMI = sd(data_a@meta.data$nCount_RNA),
        Max_UMI = max(data_a@meta.data$nCount_RNA),
        Min_UMI = min(data_a@meta.data$nCount_RNA),
        Total_genes = nrow(data_a)) %>%
        mutate_all(function(x) round(x, 2))
      
      stat_ab <- cbind(t(Before), t(After))
      colnames(stat_ab) <- c("Before", "After")
      stat_ab %>%
        datatable(
          fillContainer = F,
          options = list(scrollX = F, pageLength = 10, dom = "t"),
          extensions = "AutoFill",
          style = "auto",
          rownames = TRUE
        )
    }
  })
#   # Table of stats before and after filtering
  output$stats_before_after <- renderDT({
    if(input$filteringaction){
      filtering <- isolate(filteringresult_reactive())
      filtering
    }
  })
#   ########################### Normalization ###############################
  normalized_data <- reactive({
    if(input$normalizationaction){
      req(filtering())
      filtered_data <- isolate(filtering())
      if(input$normalizationmargin == "Features"){
        margin = 1
      } else{
        margin = 2
      }
      normalized_data <- NormalizeData(filtered_data,
                                       normalization.method = input$normalizationmethod,
                                       scale.factor = input$normalization_scalefactor,
                                       margin = margin,
                                       verbose = TRUE)

    }
  })
  # print the head of normalized data
  output$head_normalized_data <- renderDT({
    if(input$normalizationaction){
      req(normalized_data())
      normalized_data <- isolate(normalized_data())
      df <- data.frame(normalized_data[["RNA"]]@data)
      df <- df[order(df[,1], decreasing = T), ]
      df <- df[1:7, 1:10] %>%
        mutate_all(function(x) round(x, 2))
      df %>%
      datatable(
        fillContainer = F,
        options = list(scrollX = T, pageLength = 10, dom = "t"),
        extensions = "AutoFill",
        style = "auto",
        rownames = TRUE
      )
    }
  })
   
   # output$head_normalized_data <- renderPrint({
   #   if(input$normalizationaction){
   #     print('Normalization has successfully completed.')
   #   }
   # })
   
#   ############################## Feature selection ###############################
  feature_selection <- reactive({
    if(input$featureselectionaction){
      req(normalized_data())
      normalized_data <- isolate(normalized_data())
      FindVariableFeatures(normalized_data,
                           selection.method = input$featureselectionmethod,
                           nfeatures = input$featureselectionnfeatures)
    }
  })
  ## plotting features
  output$feature_selection_plots <- renderPlot({
    if(input$featureselectionaction){
      feature_selection()
      feature_selection <- isolate(feature_selection())
      top10 <- head(VariableFeatures(feature_selection), 10)
      plot1 <- VariableFeaturePlot(feature_selection,cols = c("#95A5A6","#28B463"))
      plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
      plot1 + plot2
    }
  })
#   ############################## Dimensional reduction ##############################
  dim_reduction <- reactive({
    if(input$applydimreduction){
      req(feature_selection())
      dim_reduction <- isolate(feature_selection())
      # scaling
      all.genes <- rownames(dim_reduction)
      dim_reduction <- ScaleData(dim_reduction, features = all.genes)
      # PCA
      dim_reduction <- RunPCA(dim_reduction, features = VariableFeatures(object = dim_reduction))
    }
  })
  # decide about clusters
  output$elbowplot <- renderPlot({
    if(input$applydimreduction){
      req(dim_reduction())
      dim_reduction <- isolate(dim_reduction())
      ElbowPlot(dim_reduction)
      }
  })
  # Clusters reactive
  clustering_reactive <- reactive({
    if(input$clusteringaction){
      req(dim_reduction())
      clusters <- isolate(dim_reduction())
      clusters <- FindNeighbors(clusters, dims = 1:input$number_components)
      if(input$clustringalgorithm == "Louvain"){
        algorithm <- 1
      }
      else if(input$clustringalgorithm == "LouvainMR"){
        algorithm <- 2
      }
      else if(input$clustringalgorithm == "SLM"){
        algorithm <- 3
      }
      else if(input$clustringalgorithm == "Leiden"){
        algorithm <- 4
      }
      clusters <- FindClusters(clusters,
                                    resolution = input$clusteringresolution,
                                    algorithm = algorithm,
                                    method = input$leidenmethod)
      clusters <- RunUMAP(clusters, dims = 1:input$number_components)
    }
  })
  # Cluster plot
  
  
  clusterplot_output <- reactive({
    
    if(!is.null(clustering_reactive())){
      cluster_plot <- isolate(clustering_reactive())
      
      
      qual_colals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
      col_vector <- unlist(mapply(brewer.pal, qual_colals$maxcolors, rownames(qual_colals)))
      col_vector <- col_vector[c(-3, -7, -8, -12, -14,-17,-24,-40,-41,-44,-45)]
      p <- DimPlot(cluster_plot, 
                   reduction = input$reductionmethod,
                   cols = col_vector,
                   label.size = 8,
                   pt.size = 1.5)
      ggsave("plots/Clusters.jpeg", width = 13, height = 9)
      p
      
    }
  })
  
  output$clusterplot <- renderPlot({
    if(input$clusteringaction){
      cplot <- isolate(clusterplot_output())
      cplot
    }
  })
  
  # Download clusters plot
  output$download_Cluster_Plot <- renderUI({
    if(input$clusteringaction) {
      downloadButton('clusters_plot', 'Save plot')
    }
  })
  output$clusters_plot <- downloadHandler(
    
    filename = function(){
      paste0(Sys.Date(),"_Cluster_Plot.jpeg")
    },
    content = function(file){
      file.copy("plots/Clusters.jpeg",
              file)
    }
  )
  
  
#   ########################### Find Markers ########################################
  findingmarkers_all <- reactive({
    if(!is.null(clustering_reactive()) && input$findmarkersaction){
      markers <-isolate(clustering_reactive())
      DEGs <- FindAllMarkers(markers,
                             min.pct = input$findmarkersminpct,
                             logfc.threshold = input$findmarkerslogfcthreshold,
                             test.use = input$findmarkerstestuse)
    }
  })
  
  ### Table of up regulated genes
  DEGs_up_reactive <- reactive({
    if(input$findmarkersaction) {
      DEGs <- isolate(findingmarkers_all())
      DEGs <- DEGs[which(DEGs$p_val_adj < 0.05 & abs(DEGs$avg_log2FC) > 1),]
      up <- DEGs[which(DEGs$avg_log2FC > 0), ]
    }
  })
  
  DEGs_down_reactive <- reactive({
    if(input$findmarkersaction) {
      DEGs <- isolate(findingmarkers_all())
      DEGs <- DEGs[which(DEGs$p_val_adj < 0.05 & abs(DEGs$avg_log2FC) > 1),]
      down <- DEGs[which(DEGs$avg_log2FC < 0), ]
    }
  })
  
  
  output$download_DEGs_up <- renderUI({
    if(input$findmarkersaction) {
      downloadButton('DEGs_up_table', 'Save table')
    }
  })
  
  output$DEGs_up_table <- downloadHandler(
    filename = function(){
      paste0(Sys.Date(),"_DEGs_up.csv")
    },
    content = function(file){
      DEGs_up_reactive() %>%
        write_csv(file)
    }
  ) 
  output$download_DEGs_down <- renderUI({
    if(input$findmarkersaction) {
      downloadButton('DEGs_down_table', 'Save table')
    }
  })
  
  output$DEGs_down_table <- downloadHandler(
    filename = function(){
      paste0(Sys.Date(),"_DEGs_down.csv")
    },
    content = function(file){
      DEGs_down_reactive() %>%
        write_csv(file)
    }
  ) 
  
  
  
  
  
  
  output$up_regulated <- renderDT({
    if(input$findmarkersaction){
    up <- isolate(DEGs_up_reactive())
    lastcol <- NCOL(up)
    up[, -lastcol] %>%
      datatable(extensions = "Buttons",
                options = list(autoWidth = FALSE,scrollX = TRUE,
                               buttons =c("excel", "csv"),
                               dom = "lBfrtip"),
                filter = list(position = 'top', clear = FALSE),
                rownames = TRUE)
    }
  })
  # Dot plot for up regulated
  output$dotplot_Up <- renderPlot({
    if(input$dotplot_up_action){
      findingmarkers_UP <- isolate(DEGs_up_reactive())
      clustering_reactive <- isolate(clustering_reactive())
      top10 <- findingmarkers_UP %>%
        group_by(cluster) %>%
        slice_max(avg_log2FC, n= 10) %>%
        ungroup()
      d <- DotPlot(
        clustering_reactive,
        assay = NULL,
        unique(top10$gene),
        cols = c("blue", "red"),
        col.min = -2.5,
        col.max = 2.5,
        group.by = NULL,
        split.by = NULL,
        scale.by = "radius",
        scale.max = NA,
        scale.min = NA
      ) + theme_classic() +
        theme(axis.text = element_text(family = "Times",size = 13 , colour = "black"),
              axis.text.x = element_text(family = "Times",colour = "black", size = 13, angle = 90, hjust = 1, vjust = .5),
              axis.text.y = element_text(family = "Times",colour = "black"),
              plot.subtitle = element_text(family = "Times",size = 20, colour = "black", hjust = 0.5),
              axis.title.y = element_text(family = "Times", size = rel(1.4), angle = 90),
              axis.title.x = element_text(family = "Times", size = rel(1.4), angle = 00))
      ggsave("plots/dotplot_upregulated.jpeg", width = 13, height = 9)
      d

    }
  })
  comparemarker_reactive <- reactive({
    if(input$markersvis_action){
      req(clustering_reactive())
      clustering_reactive <- isolate(clustering_reactive())
      qual_colals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
      col_vector <- unlist(mapply(brewer.pal, qual_colals$maxcolors, rownames(qual_colals)))
      col_vector <- col_vector[c(-3, -7, -8, -12, -14,-17,-24,-40,-41,-44,-45)]
      V <- VlnPlot(clustering_reactive, features = c(input$marker1, input$marker2), cols = col_vector)
      ggsave("plots/Compare_Markers.jpeg", width = 13, height = 9)
      V
    }
  })
  output$markervis_up_down <- renderPlot({
    if(input$markersvis_action){
   markers <- isolate(comparemarker_reactive())
   markers
    }
  })
  # Download compare markers
  output$download_Compare_Markers <- renderUI({
    if(input$markersvis_action) {
      downloadButton('Compare_Markers', 'Save plot')
    }
  })
  output$Compare_Markers <- downloadHandler(
    filename = function(){
      paste0(Sys.Date(),"_Compare_Markers.jpeg")
    },
    content = function(file){
      file.copy("plots/Compare_Markers.jpeg",
                file)
    }
  )
  
  
  ## Down regulated
  output$down_regulated <- renderDT({
    if(input$findmarkersaction){
      down <- isolate(DEGs_down_reactive())
      lastcol <- NCOL(down)
      down[, -lastcol] %>%
        datatable(extensions = "Buttons",
                  options = list(autoWidth = FALSE,scrollX = TRUE,
                                 buttons =c("excel", "csv"),
                                 dom = "lBfrtip"),
                  filter = list(position = 'top', clear = FALSE),
                  rownames = TRUE)
    }
  })
  # Dot plot for down regulated
  output$dotplot_Down <- renderPlot({
    if(input$dotplot_down_action){
      findingmarkers_DOWN <- isolate(DEGs_down_reactive())
      clustering_reactive <- isolate(clustering_reactive())
      top10 <- findingmarkers_DOWN %>%
        group_by(cluster) %>%
        slice_max(avg_log2FC, n= 10) %>%
        ungroup()
      d <- DotPlot(
        clustering_reactive,
        assay = NULL,
        unique(top10$gene),
        cols = c("blue", "red"),
        col.min = -2.5,
        col.max = 2.5,
        group.by = NULL,
        split.by = NULL,
        scale.by = "radius",
        scale.max = NA,
        scale.min = NA
      ) + theme_classic() +
        theme(axis.text = element_text(family = "Times",size = 13 , colour = "black"),
              axis.text.x = element_text(family = "Times",colour = "black", size = 13, angle = 90, hjust = 1, vjust = .5),
              axis.text.y = element_text(family = "Times",colour = "black"),
              plot.subtitle = element_text(family = "Times",size = 20, colour = "black", hjust = 0.5),
              axis.title.y = element_text(family = "Times", size = rel(1.4), angle = 90),
              axis.title.x = element_text(family = "Times", size = rel(1.4), angle = 00))
      ggsave("plots/dotplot_downregulated.jpeg", width = 13, height = 9)
      d

    }
  })
  
  ### Download dotplots (Up & down)
  output$download_dotplot_Up <- renderUI({
    if(input$dotplot_up_action) {
      downloadButton('Dotplot_Up', 'Save plot')
    }
  })
  output$Dotplot_Up <- downloadHandler(
    
    filename = function(){
      paste0(Sys.Date(),"_dotplot_upregulated.jpeg")
    },
    content = function(file){
      file.copy("plots/dotplot_upregulated.jpeg",
                file)
    }
  )
  
  output$download_dotplot_down <- renderUI({
    if(input$dotplot_down_action) {
      downloadButton('dotplot_down', 'Save plot')
    }
  })
  output$dotplot_down <- downloadHandler(
    
    filename = function(){
      paste0(Sys.Date(),"_dotplot_downregulated.jpeg")
    },
    content = function(file){
      file.copy("plots/dotplot_downregulated.jpeg",
                file)
    }
  )
  
  
  
  ### Volcano plot and down stream reactive.
  up_down_genes <- reactive({
    if(!is.null(findingmarkers_all())){
      deg_table <- isolate(findingmarkers_all())
      deg_table <- deg_table[, c(7, 2, 5, 1, 6)]
      colnames(deg_table) <- c("ID", "logFC", "adj.P.Val", "P.Value", "cluster")
      deg_table
    }
  })

  volcanoplot_reactive <- reactive({
    if(input$volcanoaction) {
      deg_table <- isolate(up_down_genes())
      significance <- input$significance
      zero_padj <- which(deg_table$adj.P.Val == 0)
      if(length(zero_padj) > 0) {
        minv <- min(deg_table[-zero_padj, "adj.P.Val"])
        deg_table$adj.P.Val <- deg_table$adj.P.Val + minv
        significance <- significance + minv
      }
      zero_pval <- which(deg_table$P.Value  == 0)
      if(length(zero_pval) > 0) {
        minv <- min(deg_table[-zero_pval, "P.Value"])
        deg_table$P.Value <- deg_table$P.Value + minv
        significance <- significance + minv
      }
      
      volcano <- deg_table %>%
        mutate(Pvalue = -log10(P.Value)) %>%
        mutate(AdjustedPvalue = -log10(adj.P.Val)) %>%
        mutate(DEGs = "NotSgnificance") %>%
        mutate(DEGs = ifelse(AdjustedPvalue > -log10(significance) & logFC > input$logfcup, "UpRegulated", DEGs)) %>%
        mutate(DEGs = ifelse(AdjustedPvalue > -log10(significance) & logFC < input$logfcdown, "DownRegulated", DEGs)) %>%
        mutate(labels = ifelse(DEGs != "NotSgnificance", ID, "")) %>%
        mutate(DEG = "NotSgnificance") %>%
        mutate(DEG = ifelse(Pvalue > -log10(significance) & logFC > input$logfcup , "UpRegulated", DEG)) %>%
        mutate(DEG = ifelse(Pvalue > -log10(significance) & logFC < input$logfcdown , "DownRegulated", DEG)) %>%
        mutate(label = ifelse(DEG != "NotSgnificance", ID, "")) %>%
        ggplot(aes_string(x = "logFC", y= input$pval_adjpval, label=ifelse(input$pval_adjpval == "Pvalue", "label", "labels"))) +
        labs(x= 'log2FC', y= ifelse(input$pval_adjpval == "Pvalue", "-log10(P-value)","-log10(Adjusted P-value)")) +
        geom_point(aes_string(color=ifelse(input$pval_adjpval == "Pvalue", "DEG", "DEGs"),
                              fill=ifelse(input$pval_adjpval == "Pvalue", "DEG", "DEGs"))) + scale_color_manual(values = c("#B03A2E", "#B2BABB", "#27AE60")) +
        theme_classic() +
        geom_text(check_overlap = TRUE, vjust = 0, nudge_y = 0.1) +
        theme(axis.text = element_text(family = "Times",size = 13 , colour = "black"),
              axis.text.x = element_text(family = "Times",colour = "black", size = 13),
              axis.text.y = element_text(family = "Times",colour = "black"),
              plot.subtitle = element_text(family = "Times",size = 20, colour = "black", hjust = 0.5),
              axis.title.y = element_text(family = "Times", size = rel(1.4), angle = 90),
              axis.title.x = element_text(family = "Times", size = rel(1.4), angle = 00)) +
        labs(subtitle = "Volcano plot   ")
      ggsave("plots/VolcanoPlots.jpeg", width = 13, height = 9)
      volcano
    }
  })
  
  output$volcano_plot <- renderPlot({
    if(input$volcanoaction) {
    volcano <- isolate(volcanoplot_reactive())
    volcano
    }
  })
  
  # Download compare markers
  output$download_VolcanoPlots <- renderUI({
    if(input$volcanoaction) {
      downloadButton('VolcanoPlots', 'Save plot')
    }
  })
  output$VolcanoPlots <- downloadHandler(
    filename = function(){
      paste0(Sys.Date(),"_VolcanoPlots.jpeg")
    },
    content = function(file){
      file.copy("plots/VolcanoPlots.jpeg",
                file)
    }
  )
  
  
#   ## Count matrix for heatmap 
  count_matrix <- reactive({
    if(!is.null(clustering_reactive()) & input$findmarkersaction){
      count_matrix <- isolate(clustering_reactive())
      count_matrix <- data.frame(count_matrix[["RNA"]]@scale.data)
    }
  })

#   # cluster choices for heatmap
  observeEvent(input$findmarkersaction, {
      req(clustering_reactive())
    clustering_reactive <- isolate(clustering_reactive())
      maxclust <- length(levels(Idents(clustering_reactive)))
      updateSliderInput(
        session,
        "clusterchiecheatmap", value = 2, min = 0, max = maxclust , step = 1

      )
  })
  DEGs_total <- reactive({
    if(input$findmarkersaction){
      req(up_down_genes())
      diftable <- isolate(up_down_genes())
      DEGs <- diftable[which(diftable$adj.P.Val < 0.05 & abs(diftable$logFC) > 1), ]
      up <- DEGs[which(DEGs$logFC > 0), ]
      down <- DEGs[which(DEGs$logFC < 0), ]
      up_down <- rbind(up, down)
    }
    })

#   # Choice of number of DEGs to include
  observeEvent(input$findmarkersaction, {
    req(DEGs_total())
    DEGs_total <- isolate(DEGs_total())
    maxdegs <- nrow(DEGs_total)
    updateSliderInput(
      session,
      "heatmap_degsnumber", value = 100, min = 10, max = maxdegs , step = 1

    )
  })
  heatmap_reactive <- reactive({
    if(input$heatmapaction){
      req(up_down_genes())
      req(clustering_reactive())
      req(count_matrix())
      clustering_reactive <- isolate(clustering_reactive())
      nclust <- which(Idents(clustering_reactive) == input$clusterchiecheatmap)
      cellsname <- names(Idents(clustering_reactive)[nclust])
      num_cellsname <- which(colnames(clustering_reactive) %in% cellsname)
      ex_set <- isolate(count_matrix())
      ex_set <- ex_set[, num_cellsname]
      
      deg_table <- isolate(up_down_genes())
      significance <- input$heatmapsignificance
      
      zero_padj <- which(deg_table$adj.P.Val == 0)
      if(length(zero_padj) > 0) {
        minv <- min(deg_table[-zero_padj, "adj.P.Val"])
        deg_table$adj.P.Val <- deg_table$adj.P.Val + minv
        significance <- significance + minv
      }
      zero_pval <- which(deg_table$P.Value  == 0)
      if(length(zero_pval) > 0) {
        minv <- min(deg_table[-zero_pval, "P.Value"])
        deg_table$P.Value <- deg_table$P.Value + minv
        significance <- significance + minv
      }
      
      deg_table <- deg_table %>%
        mutate(AdjustedPvalue = -log10(adj.P.Val)) %>%
        mutate(DEGs = "NotSgnificance") %>%
        mutate(DEGs = ifelse(AdjustedPvalue > -log10(significance) & logFC > input$heatmaplogfcup, "UpRegulated", DEGs)) %>%
        mutate(DEGs = ifelse(AdjustedPvalue > -log10(significance) & logFC < input$heatmaplogfcdown, "DownRegulated", DEGs)) %>%
        mutate(labels = ifelse(DEGs != "NotSgnificance", ID, ""))
      
      up <- deg_table[which(deg_table$DEG == "UpRegulated"),]
      up <-  up[order(up$logFC, decreasing = TRUE),1][1:input$heatmap_degsnumber]
      down <- deg_table[which(deg_table$DEG == "DownRegulated"),]
      down <-  down[order(down$logFC, decreasing = FALSE),1][1:input$heatmap_degsnumber]
      genelist <- melt(cbind(up, down))[,3]
      num <- which(rownames(ex_set) %in% genelist)
      network_df <- data.frame(ex_set[num,])
      rows <- rownames(network_df)
      network_df <- as.data.frame(sapply(network_df, as.numeric))
      network_df <- as.matrix(network_df)
      rownames(network_df) <- rows
      type <- colnames(network_df)
      Heatmap(network_df, name = "Expression", km = 1,
              bottom_annotation= NULL,
              use_raster = TRUE,
              show_heatmap_legend = TRUE,
              row_title = "",
              show_row_names = TRUE,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE,
              heatmap_legend_param = list(legend_direction = "vertical",
                                          legend_height = unit(50, "mm"),
                                          grid_width = unit(6, "mm"),
                                          grid_height = unit(50, "cm")
              ),
              border_gp = gpar(col = "black"),
              row_names_gp = gpar(fontsize = input$heatmapgenelabesize))
    }
  })

   # Heatmap
  output$heatmap_analysis <- renderPlot({
    if(input$heatmapaction){
      heatmap <- isolate(heatmap_reactive())
      heatmap
    }
  })
  
  output$downloadHeatmap <- downloadHandler(
    filename = function() { "Heatmap.jpeg" },
    content = function(file) {
      jpeg(file, width = 700, height = 800, units = "px", pointsize = 300, quality = 100)
      print(heatmap_reactive())
      dev.off()
    })
  
  
#   # Network
  
  genenetwork_reactive <- reactive({
    if(input$genenetworkeaction){
      clustering_reactive <- isolate(clustering_reactive())
      nclust <- which(Idents(clustering_reactive) == input$clustergenenetwork)
      cellsname <- names(Idents(clustering_reactive)[nclust])
      num_cellsname <- which(colnames(clustering_reactive) %in% cellsname)
      ex_set <- isolate(count_matrix())
      ex_set <- ex_set[, num_cellsname]

      deg_table <- isolate(up_down_genes())
      significance <- input$genenetworksignificance

      zero_padj <- which(deg_table$adj.P.Val == 0)
      if(length(zero_padj) > 0) {
        minv <- min(deg_table[-zero_padj, "adj.P.Val"])
        deg_table$adj.P.Val <- deg_table$adj.P.Val + minv
        significance <- significance + minv
      }
      zero_pval <- which(deg_table$P.Value  == 0)
      if(length(zero_pval) > 0) {
        minv <- min(deg_table[-zero_pval, "P.Value"])
        deg_table$P.Value <- deg_table$P.Value + minv
        significance <- significance + minv
      }

      deg_table <- deg_table %>%
        mutate(AdjustedPvalue = -log10(adj.P.Val)) %>%
        mutate(DEGs = "NotSgnificance") %>%
        mutate(DEGs = ifelse(AdjustedPvalue > -log10(significance) & logFC > input$genenetworklogfcup, "UpRegulated", DEGs)) %>%
        mutate(DEGs = ifelse(AdjustedPvalue > -log10(significance) & logFC < input$genenetworklogfcdown, "DownRegulated", DEGs)) %>%
        mutate(labels = ifelse(DEGs != "NotSgnificance", ID, ""))

      up <- deg_table[which(deg_table$DEG == "UpRegulated"),]
      up <-  up[order(up$logFC, decreasing = TRUE),1][1:input$genenetwork_degsnumber]
      down <- deg_table[which(deg_table$DEG == "DownRegulated"),]
      down <-  down[order(down$logFC, decreasing = FALSE),1][1:input$genenetwork_degsnumber]
      genelist <- melt(cbind(up, down))[,3]
      num <- which(rownames(ex_set) %in% genelist)
      network_df <- data.frame(ex_set[num,])
      rows <- rownames(network_df)
      network_df <- as.data.frame(sapply(network_df, as.numeric))
      network_df <- as.matrix(network_df)
      rownames(network_df) <- rows
      if(NROW(network_df) > 0){
      set.seed(123)
      g <- graph.adjacency(
        as.matrix(as.dist(cor(t(network_df), method=input$genenetworkcormethod))),
        mode=input$genenetworkmode,
        weighted = TRUE,
        diag= FALSE
      )

      g <- simplify(g,
                    remove.multiple = TRUE,
                    remove.loops = TRUE)
      # Colour negative correlation edges as blue
      E(g)[which(E(g)$weight<0)]$color <- "darkblue"
      # Colour positive correlation edges as red
      E(g)[which(E(g)$weight>0)]$color <- "darkred"
      E(g)$weight <- abs(E(g)$weight)
      # Remove edges below absolute Pearson correlation 0.8
      g <- delete_edges(g, E(g)[which(E(g)$weight< input$genenetworkcordegree)])
      g <- delete.vertices(g, degree(g) < 3)
      V(g)$name <- V(g)$name
      V(g)$shape <- "sphere"
      V(g)$color <- "skyblue"
      V(g)$vertex.frame.color <- "white"
  

        #Scale the size of the vertices to be proportional to the level of expression of each gene represented by each vertex
        #Multiply scaled vales by a factor of 10
        scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
        vSizes <- (scale01(apply(network_df, 1, mean)) + 1.0) * 10
        # Amplify or decrease the width of the edges
        edgeweights <- E(g)$weight * 2
        # Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
        mst <- mst(g, algorithm="prim")
        mst.communities <- cluster_edge_betweenness(mst, weights=NA, directed=FALSE)
        mst.clustering <- make_clusters(mst, membership=mst.communities$membership)
        V(mst)$color <- mst.communities$membership + 1
        if(input$networklayout == "layout.auto"){
          layout <- layout.auto
        }
        else if(input$networklayout == "layout.random"){
          layout <- layout.random
        }
        else if(input$networklayout == "layout.circle"){
          layout <- layout.circle
        }
        else if(input$networklayout == "layout.fruchterman.reingold"){
          layout <- layout.fruchterman.reingold
        }
        else if(input$networklayout == "layout.kamada.kawai"){
          layout <- layout.kamada.kawai
        }
        else if(input$networklayout == "layout.spring"){
          layout <- layout.spring
        }
        else if(input$networklayout == "layout.reingold.tilford"){
          layout <- layout.reingold.tilford
        }
        else if(input$networklayout == "layout.fruchterman.reingold.grid"){
          layout <- layout.fruchterman.reingold.grid
        }
        else if(input$networklayout == "layout.lgl"){
          layout <- layout.lgl
        }
        else if(input$networklayout == "layout.graphopt"){
          layout <- layout.graphopt
        }
        else if(input$networklayout == "layout.mds"){
          layout <- layout.mds
        }
        else if(input$networklayout == "layout.svd"){
          layout <- layout.svd
        }
        #mst.clustering,
      p <- plot(
        mst,
        layout=layout,
        edge.curved= TRUE,
        vertex.size=vSizes,
        vertex.label.dist= - input$genenetworklabledist,
        vertex.label.color="black",
        asp=FALSE,
        vertex.label.cex= input$genenetworklablesize,
        edge.width=edgeweights,
        edge.arrow.mode=0,
        main="")
      p
      }
    }
  })
  
  output$genenetwork_analysis <- renderPlot({
    if(input$genenetworkeaction){
      network <- isolate(genenetwork_reactive())
      network

    }
  })
  
  output$downloadNetwork <- downloadHandler(
    filename = function() { "GeneNetwork.jpeg" },
    content = function(file) {
      jpeg(file, width = 700, height = 800, units = "px", pointsize = 300, quality = 100)
      print(genenetwork_reactive())
      dev.off()
    })
#   ################################# Pathway analysis ####################################
  pathway <- reactive({
    if(input$pathwayaction){
      req(up_down_genes())
      diftable <- isolate(up_down_genes())
      DEGs <- diftable[which(diftable$adj.P.Val < input$pathwaysignificance & abs(diftable$logFC) > 1), ]
      up <- DEGs[which(DEGs$logFC > 0), ]
      down <- DEGs[which(DEGs$logFC < 0), ]
      up_down <- rbind(up, down)

      up_down$entrez <- mapIds(org.Hs.eg.db,
                               keys=up_down$ID,
                               column="ENTREZID",
                               keytype="SYMBOL",
                               multiVals="first")

      up_down$name <- mapIds(org.Hs.eg.db,
                             keys=up_down$ID,
                             column="GENENAME",
                             keytype="SYMBOL",
                             multiVals="first")
      # Load KEGG data
      data(kegg.sets.hs)
      data(sigmet.idx.hs)
      kegg.sets.hs <- kegg.sets.hs[sigmet.idx.hs]

      foldchanges <- up_down$logFC
      names(foldchanges) <- up_down$entrez
      kegg_up_down <- gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
    }
  })
  # KEGG for up
  KEGG_Pathaway_reactive <- reactive({
    if(input$pathwayaction){
    req(pathway())
    kegg_up_down <- isolate(pathway())
    pathway_up <- data.frame(kegg_up_down$greater)
    pathway_up$ID <- str_extract(rownames(pathway_up), "hsa[0-9]+(?= )")
    pathway_up$Pathway <- gsub("hsa[0-9]+ ", "", rownames(pathway_up))
    pathway_up <- pathway_up[order(pathway_up$q.val, decreasing = FALSE), c(7, 8, 5, 4)]
    colnames(pathway_up) <- c("ID", "Pathway", "n_Gene", "FDR")
    rownames(pathway_up) <- NULL
    
    pathway_up$mode <- "Upregulated"
    pathway_up <- pathway_up[which(pathway_up$FDR < input$pathwaysignificance), ]
    
    pathway_down <- data.frame(kegg_up_down$less)
    pathway_down$ID <- str_extract(rownames(pathway_down), "hsa[0-9]+(?= )")
    pathway_down$Pathway <- gsub("hsa[0-9]+ ", "", rownames(pathway_down))
    pathway_down <- pathway_down[order(pathway_down$q.val, decreasing = FALSE), c(7, 8, 5, 4)]
    colnames(pathway_down) <- c("ID", "Pathway", "n_Gene", "FDR")
    rownames(pathway_down) <- NULL
    pathway_down$mode <- "Downregulated"
    pathway_down <- pathway_down[which(pathway_down$FDR < input$pathwaysignificance), ]
    pathways <- rbind(pathway_up[1:input$numberoftoppathways, ], pathway_down[1:input$numberoftoppathways, ])
    pathways <- na.omit(pathways)
    g <- pathways %>%
      ggplot(aes(reorder(Pathway, n_Gene, sum), n_Gene)) +
      geom_col(aes(color = FDR, fill = FDR), width = 0.4) +
      coord_flip() + xlab("Pathway") + ylab("Enriched genes") +
      theme_classic() +
      scale_color_gradient(low = "#D7DBDD", high = "#186A3B") +
      scale_fill_gradient(low = "#D7DBDD", high = "#186A3B") +
      theme(axis.text = element_text(family = "Times",size = 10 , colour = "black", angle = 00),
            axis.text.x = element_text(family = "Times",colour = "black", size = 14),
            axis.text.y = element_text(family = "Times",colour = "black", size = 14),
            plot.subtitle = element_text(family = "Times",size = 24, colour = "black", hjust = 0.5),
            axis.title.y = element_text(family = "Times", size = 18, angle = 90, hjust = 0.5),
            axis.title.x = element_text(family = "Times", size = 18, angle = 00)) +
      labs(subtitle = "")
    ggsave(filename = paste0("plots/KEGG_Pathways.jpeg"),
           width = 10,
           height = 7)
    g
  }
  })
  
  
  output$kegg_up_down_pathways <- renderPlot({
    if(input$pathwayaction){
      pathways <- isolate(KEGG_Pathaway_reactive())
      pathways
    }
  })
  
  # Download Pathway plot
  output$download_KEGGPathway <- renderUI({
    if(input$pathwayaction) {
      downloadButton('KEGGPathway', 'Save plot')
    }
  })
  output$KEGGPathway <- downloadHandler(
    filename = function(){
      paste0(Sys.Date(),"_KEGG_Pathways.jpeg")
    },
    content = function(file){
      file.copy("plots/KEGG_Pathways.jpeg",
                file)
    }
  )
  
  ### Reactive choices for up ides
  choices_up_kegg_reactive <- reactive({
      req(pathway())
      kegg_up_down <- isolate(pathway())
      pathway_up <- data.frame(kegg_up_down$greater)
      pathway_up$ID <- str_extract(rownames(pathway_up), "hsa[0-9]+(?= )")
      pathway_up$Pathway <- gsub("hsa[0-9]+ ", "", rownames(pathway_up))
      pathway_up <- pathway_up[order(pathway_up$q.val, decreasing = FALSE), c(7, 8, 5, 4)]
      colnames(pathway_up) <- c("ID", "Pathway", "n_Gene", "FDR")
      rownames(pathway_up) <- NULL
      kegg_up_ids <- pathway_up[1:5,1]
  })

#   # Choices of top 5 for up kegg
  observeEvent(input$pathwayaction,{
    kegg_up_ids <- choices_up_kegg_reactive()
        updateSelectInput(
          session,
          "kegg_up_path_graphs",
          choices = kegg_up_ids)
  })
# 
#   # Downloading up kegg 
  output$downloading_up_kegg <- renderImage({
    req(input$kegg_up_path_graphsaction)
      diftable <- isolate(up_down_genes())
      # DEGs
      DEGs <- diftable[which(diftable$adj.P.Val < 0.05 & abs(diftable$logFC) > 1), ]
      up <- DEGs[which(DEGs$logFC > 0), ]
      down <- DEGs[which(DEGs$logFC < 0), ]
      up_down <- rbind(up, down)

      up_down$entrez <- mapIds(org.Hs.eg.db,
                               keys=up_down$ID,
                               column="ENTREZID",
                               keytype="SYMBOL",
                               multiVals="first")

      up_down$name <- mapIds(org.Hs.eg.db,
                             keys=up_down$ID,
                             column="GENENAME",
                             keytype="SYMBOL",
                             multiVals="first")
      # Load KEGG data
      data(kegg.sets.hs)
      data(sigmet.idx.hs)
      kegg.sets.hs <- kegg.sets.hs[sigmet.idx.hs]

      foldchanges <- up_down$logFC
      names(foldchanges) <- up_down$entrez
      kegg_up_down <- gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
      pathway_up <- data.frame(kegg_up_down$greater)
      pathway_up$ID <- str_extract(rownames(pathway_up), "hsa[0-9]+(?= )")
      pathway_up$Pathway <- gsub("hsa[0-9]+ ", "", rownames(pathway_up))
      pathway_up <- pathway_up[order(pathway_up$q.val, decreasing = FALSE), c(7, 8, 5, 4)]
      colnames(pathway_up) <- c("ID", "Pathway", "n_Gene", "FDR")
      kegg_top5_up_ids <- pathway_up[1:5, 1]

      plot_pathway <- function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)
      tmp <- sapply(kegg_top5_up_ids, function(pid) pathview(gene.data=foldchanges,
                                                             pathway.id=pid,
                                                             species="hsa",
                                                             same.layer = F))
      id <- input$kegg_up_path_graphs
      path <- paste0( id, ".pathview.png")
      expr = list(src = path, width = "100%", height = "180%")
    #} "plots/KEGG_UP/",
  }, deleteFile = FALSE)

  ### Reactive choices for down ides
  choices_down_kegg_reactive <- reactive({
      req(pathway())
      kegg_up_down <- isolate(pathway())
      pathway_down <- data.frame(kegg_up_down$less)
      pathway_down$ID <- str_extract(rownames(pathway_down), "hsa[0-9]+(?= )")
      pathway_down$Pathway <- gsub("hsa[0-9]+ ", "", rownames(pathway_down))
      pathway_down <- pathway_down[order(pathway_down$q.val, decreasing = FALSE), c(7, 8, 5, 4)]
      colnames(pathway_down) <- c("ID", "Pathway", "n_Gene", "FDR")
      rownames(pathway_down) <- NULL
      kegg_down_ids <- pathway_down[1:5,1]
  })

#   # Choices of top 5 for up kegg
  observeEvent(input$pathwayaction,{
    kegg_down_ids <- choices_down_kegg_reactive()
    updateSelectInput(
      session,
      "kegg_down_path_graphs",
      choices = kegg_down_ids)
  })
  # Downloading DOWN kegg
  output$downloading_down_kegg <- renderImage({
    req(input$kegg_down_path_graphsaction)
      diftable <- isolate(up_down_genes())
      # DEGs
      DEGs <- diftable[which(diftable$adj.P.Val < 0.05 & abs(diftable$logFC) > 1), ]
      up <- DEGs[which(DEGs$logFC > 0), ]
      down <- DEGs[which(DEGs$logFC < 0), ]
      up_down <- rbind(up, down)

      up_down$entrez <- mapIds(org.Hs.eg.db,
                               keys=up_down$ID,
                               column="ENTREZID",
                               keytype="SYMBOL",
                               multiVals="first")

      up_down$name <- mapIds(org.Hs.eg.db,
                             keys=up_down$ID,
                             column="GENENAME",
                             keytype="SYMBOL",
                             multiVals="first")
      # Load KEGG data
      data(kegg.sets.hs)
      data(sigmet.idx.hs)
      kegg.sets.hs <- kegg.sets.hs[sigmet.idx.hs]

      foldchanges <- up_down$logFC
      names(foldchanges) <- up_down$entrez
      kegg_up_down <- gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
      pathway_down <- data.frame(kegg_up_down$less)
      pathway_down$ID <- str_extract(rownames(pathway_down), "hsa[0-9]+(?= )")
      pathway_down$Pathway <- gsub("hsa[0-9]+ ", "", rownames(pathway_down))
      pathway_down <- pathway_down[order(pathway_down$q.val, decreasing = FALSE), c(7, 8, 5, 4)]
      colnames(pathway_down) <- c("ID", "Pathway", "n_Gene", "FDR")
      kegg_top5_down_ids <- pathway_down[1:5, 1]

      plot_pathway <- function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)
      tmp <- sapply(kegg_top5_down_ids, function(pid) pathview(gene.data=foldchanges,
                                                             pathway.id=pid,
                                                             species="hsa",
                                                             same.layer = F))
      id <- input$kegg_down_path_graphs
      path <- paste0("plots/KEGG_DOWN/", id, ".pathview.png")
      expr = list(src = path, width = "100%", height = "180%")
    
  }, deleteFile = FALSE)

  # 
  GSEA_analysis_reactive <- reactive({
    if(input$gsea_action){
      req(up_down_genes())
    diftable <- isolate(up_down_genes())
    DEGs <- diftable[which(diftable$adj.P.Val < 0.05 & abs(diftable$logFC) > 1), ]
    up <- DEGs[which(DEGs$logFC > 0), ]
    down <- DEGs[which(DEGs$logFC < 0), ]
    up_down <- rbind(up, down)
    H <- msigdbr(species = "Homo sapiens", category = input$gseacategory) #
    H.symbol.ls <- H %>% 
      dplyr::select(gs_name, gene_symbol) %>% 
      group_by(gs_name) %>% 
      summarise(all.genes = list(unique(gene_symbol))) %>% 
      deframe()
    # Just gene names and logFC
    FC <- up_down[,c(1,2)]
    ##format for gsea 
    FC.vec <- FC$logFC
    names(FC.vec) <- FC$ID
    #Run GSEA
    gsea.H <- fgseaSimple(pathways = H.symbol.ls,
                          stats = FC.vec,
                          scoreType = "std",
                          nperm=1000)
    
    my_pal <- c("#1B9E77","#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#9A7D0A")
    gsea.H$pathway <- gsub("HALLMARK_","", gsea.H$pathway)
    gsea.H$pathway <- gsub("_"," ", gsea.H$pathway)
    gsea.H <- gsea.H %>%
      mutate(Significance = "NotSignificance") %>%
      mutate(Significance = ifelse(padj <= input$gseasignificance & NES > 0, "Upregulated", Significance)) %>%
      mutate(Significance = ifelse(padj <= input$gseasignificance & NES < 0, "Downregulated", Significance))
    gsea.H <- gsea.H[which(gsea.H$Significance != "NotSignificance"), ]
    gsea.H <- na.omit(gsea.H)
    gsea.H <- gsea.H[1:input$gseanumber, ]
   p <- gsea.H %>%
      ggplot(aes(x=reorder(pathway, NES),
                 y=NES)) +
      geom_col(width = 0.3, aes(color= Significance, fill = Significance)) +
      labs(x= 'Gene set', y= "Normalized enrichment score (NES)") +
      theme_classic() +
      scale_color_manual(values=c("#1B9E77", "#7570B3")) +
      scale_fill_manual(values=c("#1B9E77", "#7570B3")) +
      coord_flip() +
      theme(axis.text = element_text(family = "Times",size = 24 , colour = "black"),
            axis.text.x = element_text(family = "Times",colour = "black", size = 11),
            axis.text.y = element_text(family = "Times",colour = "black", size = 10),
            plot.subtitle = element_text(family = "Times",size = 16, colour = "black", hjust = 0.5),
            axis.title.y = element_text(family = "Times", size = 16, angle = 90),
            axis.title.x = element_text(family = "Times", size = 16, angle = 00),
            legend.text = element_text(size = 10, family = "Times"),
            legend.title = element_text(size = 20, family = "Times")) +
      labs(subtitle = paste0("GSEA"))
    ggsave(filename = paste0("plots/GSEA.jpeg"),width = 10, height = 7)
    p
    }
  })

  output$GSEA_Plot <- renderPlot({
    if(input$gsea_action){
      gsea <- isolate(GSEA_analysis_reactive())
      gsea
    }
  })
  
  # Download GSEA plot
  output$download_GSEA <- renderUI({
    if(input$gsea_action) {
      downloadButton('GSEAPLOT', 'Save plot')
    }
  })
  output$GSEAPLOT <- downloadHandler(
    filename = function(){
      paste0(Sys.Date(),"_GSEA.jpeg")
    },
    content = function(file){
      file.copy("plots/GSEA.jpeg",
                file)
    }
  )

})