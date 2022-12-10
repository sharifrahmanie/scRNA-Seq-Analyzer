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
require(shinyWidgets)
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
dashheader <- dashboardHeader(title = "",titleWidth = 130,
                              tags$li(class = "dropdown",fluidRow(
                                
                                tags$div(style = "margin-top: 15px;margin-right: -160px;font-size: 15px",
                                         prettySwitch("togglemode", "Night mode", value = FALSE, fill = TRUE, status = "primary")
                                )
                              )),
                              
                              tags$li(class = "dropdown",  tags$a(href="mailto:rahami.biotech@gmail.com", icon("envelope"), "Contact", target= "_blank")),
                              
                              tags$li(class = "dropdown", tags$a(href="https://www.linkedin.com/in/sharifrahmanie/", icon("linkedin"), "Linkedin", target= "_blank"))
                              
                              
                              
)


dashSidebar <- dashboardSidebar(width = 130,collapsed = F,
  tags$head(tags$link(rel="shortcut icon", href="cell.png", height ="100%", width = "100%"),
            
            #################### Hide toggle
            tags$script(
              HTML(#code for hiding sidebar tabs 
                "Shiny.addCustomMessageHandler('manipulateMenuItem1', function(message)
        {
        var aNodeList = document.getElementsByTagName('a');

        for (var i = 0; i < aNodeList.length; i++) 
        {
        if(aNodeList[i].getAttribute('data-toggle') == message.toggle && aNodeList[i].getAttribute('role') == message.role) 
        {
        if(message.action == 'hide')
        {
        aNodeList[i].setAttribute('style', 'display: none;');
        } 
  
        };
        }
        });"
              ),
        #########################################
            
            )),
  shinyjs::useShinyjs(),
  use_darkmode(),
  sidebarMenu(id = "sidebar",
              menuItem(
                text = "scRNA-Seq",
                tabName = "scrna-seq",
                icon = icon("dna")
              ),
              uiOutput('style_tag')
  )
)

dashBody <- dashboardBody(
  useShinyjs(),
  ####################################################
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "app.css")),
  #########################################################
  ################################################## Home #######################################################
  tabItems(
  ################################################### SCRNA-Seq #######################################################
  tabItem(
    tabName = "scrna-seq",
    fluidRow(
      div(id = "inpupanelpos", column(3, 
             box(background = "navy", width = "100%", height = 80,
                 column(3,tags$div(style = "margin-top: -5px;margin-left: -16px;font-size: 20px",tags$img(src= "cell.png", style="text-align: right;", height ="70", width ="70"))),
                 br(),
                 column(9,tags$b(style="font-size: 20px;","Single-cell"))
                 
             ),
             box(status = "primary", collapsible = F, solidHeader = F, width = "100%", background = "navy",
                 tagList(
                  
                   includeScript("www/text.js"),
                              tags$div(class="form-group shiny-input-container", 
                                       tags$div(tags$label("Choose folder", class="btn btn-primary",
                                                           tags$input(id = "fileIn", webkitdirectory = TRUE, type = "file", style="display: none;", onchange="pressed()"))),
                                       tags$label("No folder choosen", id = "noFile"),
                                       tags$div(id="fileIn_progress", class="progress progress-striped active shiny-file-input-progress",
                                                tags$div(class="progress-bar")
                                       )     
                              ),
                 ),
                 HTML("<script type='text/javascript' src='getFolders.js'></script>"),
                 column(6,div(id = "barcodescolpos",numericInput(inputId = "barcodescol", label = "Barcode column", value = 1))),
                 column(6,div(id = "featurescolpos",numericInput(inputId = "featurescol", label = "Feauture column", value = 2))),
                 column(6, div(id = "projectnamepos", textInput(inputId = "projectname", label = "Project name", value = "Myproject"))),
                 column(6, div(id = "start_seuratpos", actionButton("start_seurat","Analyse", width = "100%"))),
             ))),
      column(6, box(status = "primary", collapsible = F, solidHeader = F, width = "100%", background = "navy",
                    tabsetPanel(type = "pills",
                                tabPanel(title = "QC",
                                         fluidRow(
                                           column(12,
                                                  tags$img(src="qc.png", height = 268, width = "100%")))),
                                tabPanel(title = "Clusters",
                                         fluidRow(
                                           column(12,
                                                  tags$img(src="clustering.png", height = 268, width = "100%")))),
                                tabPanel(title = "Volcano",
                                         fluidRow(
                                           column(12,
                                                  tags$img(src="volcano.png", height = 268, width = "100%")))),
                                tabPanel(title = "Heatmap",
                                         fluidRow(
                                           column(12,
                                                  tags$img(src="heatmap.png", height = 268, width = "100%")))),
                                tabPanel(title = "Network",
                                         fluidRow(
                                           column(12,
                                                  tags$img(src="network.png", height = 268, width = "100%")))),
                                tabPanel(title = "Pathway",
                                         fluidRow(
                                           column(12,
                                                  tags$img(src="pathway.png", height = 268, width = "100%")))),
                                tabPanel(title = "Help", 
                                         fluidRow(
                                           
                                           column(12,
                                                column(12, br()),
                                                  tags$b("Crucial steps:"),
                                                  tags$ul(
                                                    tags$li("Users should only provide ", tags$b("three"), "files in ", tags$b('a separate folder.')),
                                                    tags$li("Having more than three files will raise an error."),
                                                    tags$li("Files should be named", tags$b("matrix.mtx"),",", tags$b("barcodes.tsv"), ", and ",tags$b('genes.tsv'),"."),
                                                    tags$li("Compressed files with ", tags$b("gz"), "extension are also accepted."),
                                                    tags$li("Having files with other names will raise an error."),
                                                    tags$li("R is", tags$b("case sensitive"),"; therefore, all files should be named lower case"),
                                                    tags$li("When uploading your files, you should specify the " , tags$b("column of barcodes"), " and ", tags$b("genes.")),
                                                    tags$li("Unfortunately, using the free version of Shinyapps, this application only supports ", tags$b("5 Mb")," file size."),
                                                    tags$li("Example files can be downloaded from my ", tags$a(href="https://github.com/sharifrahmanie/scRNA-Seq-Analyzer", "github", target= "_blank"), "repository."),
                                                    
                                                  ),
                                                column(12, br()),
                                                  ))),
                                
                                
                    )
      )),
      column(3,
             box(background = "navy", width = "100%", height = 80,
                 column(3,tags$div(style = "margin-top: -5px;margin-left: -16px;font-size: 20px",tags$img(src= "rna.png", style="text-align: right;", height ="70", width ="70"))),
                 br(),
                 column(9,tags$b(style="font-size: 20px;","RNA-Seq"))
                 
             ),
             
             
             div(id = "describtionpos", box(status = "primary", collapsible = F, solidHeader = F, width = "100%", background = "navy", height = 234,
                 
                 p(style="text-align: justify;text-justify: inter-word;", tags$b("scRNA-SeqAnalyzer"), "is a free shiny web application that uses the Seurat package to automatically analyze single-cell RNA-Seq
                   data produced by 10X genomics.",tags$a(href="https://www.youtube.com/watch?v=9O3ejxw-81Q", "Tutorial.", target= "_blank")),
                  p(style="text-align: justify;text-justify: inter-word;","Contributor:
                   Edris Sharif Rahmani. He has developed another shiny web application,",
                   tags$a(href= "https://bioinformatics-ml.shinyapps.io/deepbc/", target= "_blank","DeepBC,"),
                   "to simplify performing bioinformatics and ML."),
                 # column(12, br()),
                 
             ))),
    ),
    column(12, br()),
    ## check errors
    div(id = "errorpos",conditionalPanel(condition = "input.start_seurat && output.files_sanity_checking == false",
                     infoBox(title = "Error", subtitle = "It must be files' name, or the number of files in the folder (> 3). Please check the help panel.",color = "olive",icon = icon("bug"),fill = TRUE))),
    conditionalPanel(condition = c("input.start_seurat && output.files_sanity_checking"),
                     tabsetPanel(type = "tabs",
                                 tabPanel(title = "Table", icon = icon("table"),
                                          fluidRow(
                                            column(12,box(width = "100%",status = "primary",solidHeader = T, collapsible = T,
                                                          column(12,DTOutput("metadata"))
                                            )
                                            ))),
                                 tabPanel(title = "QC", icon = icon("chart-line"),
                                          box(width = "100%",status = "primary",solidHeader = T, collapsible = T,
                                              column(5, column(6,uiOutput("download_box_plot"), column(6))),
                                              column(2, actionButton(inputId = "qcplotaction", label = "Display plots", width = "100%")),
                                              column(5),
                                              column(12,br()),
                                              column(12,plotOutput("qcplot0",width = "100%",height = 600)),
                                              column(12,br()),
                                              column(4,plotOutput("qcplot1",width = "100%",height = 600)),
                                              column(4,plotOutput("qcplot2",width = "100%",height = 600)),
                                              column(4,plotOutput("qcplot3",width = "100%",height = 600)),
                                              column(12,br()),
                                              column(4,plotOutput("qcplot4",width = "100%",height = 600)),
                                              column(4,plotOutput("qcplot5",width = "100%",height = 600)),
                                              column(4,plotOutput("qcplot6",width = "100%",height = 600)),
                                              )),
                                 tabPanel(title = "Filtering", icon = icon("check"),
                                          box(width = "100%",status = "primary",solidHeader = T, collapsible = T,
                                              column(4, sidebarPanel(width = "100%",
                                                                     
                                                                     sliderInput("filteringlessfeature", width = 300,"Discard cells with less feature than:",min =1,max = 500,value = 200,step = 1),
                                                                     sliderInput("filteringmorefeature", width = 300, "Discard cells with more feature than:",min = 1000, max = 20000, value = 2500,step = 1),
                                                                     sliderInput("filteringmtgenes", width = 300,"Keep cells with MT genes less than: (%)", min = 0,max = 100,value = 5,step = 1),
                                                                     sliderInput("filteringgeneswithzero", width = 300,"Discard a gene with more zeroes than:", min = 0,max = 100,value = 10,step = 1),
                                                                     actionButton(inputId = "filteringaction", label = "Apply filters", width = "100%"))),
                                              conditionalPanel(condition = "input.filteringaction",
                                              column(4, column(12, box(title = "Filtering stats", width = "100%",status = "primary",solidHeader = T, collapsible = T, DTOutput("stats_before_after"))))))),
                                 
                                 tabPanel(title = "Normalization", icon = icon("angular"),
                                          box(width = "100%",status = "primary",solidHeader = T, collapsible = T,
                                          column(3, box(title = "Parameters", width = "100%",status = "primary",solidHeader = T, collapsible = T, 
                                                        selectInput(inputId = "normalizationmethod",width = 300, label = "Method", choices = c("LogNormalize", "CLR", "RC"),multiple = F),
                                                        sliderInput(inputId = "normalization_scalefactor", width = 300, label = "Scale factor", min = 5000, max = 50000, value = 10000, step = 1000),
                                                        selectInput(inputId = "normalizationmargin", width = 300, label = "Margin", choices = c("Features", "Cells"),multiple = F),
                                                        actionButton(inputId = "normalizationaction",label = "Normalize", width = "100%"), column(12,br()))),
                                          conditionalPanel(condition = "input.normalizationaction",
                                          column(9, box(title = "Head of normalized data", width = "100%",status = "primary",solidHeader = T, collapsible = T, DTOutput("head_normalized_data")))))),
                                 tabPanel(title = "Feature selection", icon = icon("check-circle"),
                                          box(width = "100%",status = "primary",solidHeader = T, collapsible = T,
                                              br(),
                                              br(),
                                              column(12, box(title = "Parameters", width = "100%",status = "primary",solidHeader = T, collapsible = T, 
                                                             column(4,selectInput(inputId = "featureselectionmethod", label = "Method", choices = c("vst", "mvp", "disp"), multiple = F)),
                                                             column(4,sliderInput(inputId = "featureselectionnfeatures", label = "Number of features", min = 1000, max = 5000, value = 2000, step = 100)),
                                                             column(4,div(id = "featureselectionactionbt",actionButton(inputId = "featureselectionaction",label = "Start", width = "100%")),
                                                                    br(),
                                                                    br())
                                                            )    
                                              ),
                                              conditionalPanel(condition = "input.featureselectionaction",
                                              column(12,box(title = "Top 10 features", width = "100%",status = "primary",solidHeader = T, collapsible = T, plotOutput("feature_selection_plots")))))),
                                 tabPanel(title = "Clustering", icon = icon("sitemap"),
                                          box(width = "100%",status = "primary",solidHeader = T, collapsible = T,
                                              column(5), column(2,actionButton(inputId = "applydimreduction", label = "Apply PCA", width = "100%")), column(2),column(12, br()),
                                              column(12, br()),
                                              column(12, br()),
                                              column(12, br()),
                                              conditionalPanel(condition = "input.applydimreduction",
                                              column(12, box(title = "Elbow method", width = "100%",status = "primary",solidHeader = T, collapsible = T,
                                                             column(3, sidebarPanel(width = "100%",
                                                                                    sliderInput(inputId = "number_components", label = "Number of PCs", min = 2, max = 20, value = 10, step = 1),
                                                                                    sliderInput(inputId = "clusteringresolution", label = "Resolution", min = 0.4, max = 1.2, value = 0.4, step = 0.1),
                                                                                    selectInput(inputId = "clustringalgorithm", label = "Algorithm", choices = c("Louvain", "LouvainMR", "SLM", "Leiden"), multiple = F),
                                                                                    # pip install leidenalg
                                                                                    selectInput(inputId = "leidenmethod", label = "Leiden method", choices = c("matrix", "igraph"), multiple = F),
                                                                                    selectInput(inputId = "reductionmethod", label = "Reduction method", choices = c("umap", "pca"), multiple = F),
                                                                                    actionButton(inputId = "clusteringaction", label = "Find clusters", width = "100%")
                                                             )),
                                                             column(9,plotOutput("elbowplot", height = 500)))),
                                                    ),
                                              conditionalPanel(condition = "input.clusteringaction",
                                                     column(12, box(title = "Clusters", width = "100%", height = 700, status = "primary",solidHeader = T, collapsible = T,
                                                                    column(12, br()),
                                                                    uiOutput("download_Cluster_Plot"),
                                                                    column(12, br()),         
                                                         plotOutput("clusterplot", height = 500, width = "100%")))))
                                              ),
                                 tabPanel(title = "DGE", icon = icon("table"),
                                          box(width = "100%",status = "primary",solidHeader = T, collapsible = T,
                                              box(title = "DGE tables",width = "100%", status = "primary",solidHeader = T, collapsible = T,
                                              column(12,
                                                     box(title = "Parameters",width = "100%", status = "primary",solidHeader = T, collapsible = T,
                                                         column(3, sliderInput("findmarkerslogfcthreshold", width = 300,"Least logFC",min = 0.01,max = 5,value = 1,step = 0.01)),
                                                         column(3, pickerInput(inputId = "findmarkerstestuse", label = "Test", choices = c("DESeq2", "wilcox", "bimod", "roc", "t", "poisson", "LR", "MAST"), multiple = F)),
                                                         column(3,sliderInput("findmarkersminpct", width = 300,"Minimum apparance",min =0.01,max = 1,value = 0.25,step = 0.01)),
                                                         column(3,div(id = "findmarkersactionbt", actionButton(inputId = "findmarkersaction", label = "Find DEGs", width = "100%"))))),
                                              conditionalPanel(condition = "input.findmarkersaction",
                                              column(12, box(title = "Up regulated genes", width = "100%", height = 1450,status = "primary",solidHeader = T, collapsible = T,
                                                             column(12, br()),
                                                             uiOutput("download_DEGs_up"),
                                                             column(12, br()),
                                                             DTOutput("up_regulated"),
                                                             column(12, br()),
                                                             column(12,column(5), column(2,actionButton(inputId = "dotplot_up_action", label = "Dot plot", width = "100%")), column(5)),
                                                     column(12, 
                                                            column(12, br()),
                                                            uiOutput("download_dotplot_Up"),
                                                            column(12, br()),
                                                            plotOutput("dotplot_Up",  height = 400, width = "100%"))),
                                                     ),
                                              column(12, box(title = "Down regulated genes", width = "100%", height = 1350, status = "primary",solidHeader = T, collapsible = T,
                                                             column(12, br()),
                                                             uiOutput("download_DEGs_down"),
                                                             column(12, br()),
                                                             DTOutput("down_regulated"),
                                                             column(12, br()),
                                                             column(12,column(5), column(2,actionButton(inputId = "dotplot_down_action", label = "Dot plot", width = "100%")), column(5)),
                                                             column(12, 
                                                                    column(12, br()),
                                                                    uiOutput("download_dotplot_down"),
                                                                    column(12, br()),
                                                                    plotOutput("dotplot_Down", height = 400, width = "100%")),
                                                             )))),
                                              
                                              conditionalPanel(condition = "input.findmarkersaction",
                                              box(title = "Compare markers",width = "100%", height = 700, status = "primary",solidHeader = T, collapsible = T,
                                                  column(12, column(4), column(2, textInput(inputId = "marker1", label = "Marker 1", value = "MS4A1", width = "100%")),column(2, textInput(inputId = "marker2", label = "Marker 2", value = "CD79A", width = "100%")), column(4)),
                                                  column(12, column(5), column(2, actionButton(inputId = "markersvis_action", label = "Find markers", width = "100%")), column(5)),
                                                  column(12, 
                                                         column(12, br()),
                                                         uiOutput("download_Compare_Markers"),
                                                         column(12, br()),
                                                         plotOutput("markervis_up_down"))
                                              ),
                                            
                                              
                                              box(title = "Volcano plot",width = "100%",height =870, status = "primary",solidHeader = T, collapsible = T,
                                                column(12, column(5),column(2, actionButton(inputId = "volcanoaction", label = "Display plot", width = "100%")),column(5)), 
                                                
                                                column(3, sidebarPanel(width = "100%",selectInput("pval_adjpval", width = 300,"P-value or Adjusted P-value",choices = c("Pvalue", "AdjustedPvalue")),
                                                           sliderInput("significance", width = 300,"Change the significance level",min = 0.01,max = 0.1,value = 0.05,step = 0.01),
                                                           sliderInput("logfcup", width = 300, "Change the log2FC (Up)",min = 0, max = 10, value = 1,step = 0.1),
                                                           sliderInput("logfcdown", width = 300,"Change the log2FC (Down)", min = -10,max = 0,value = -1,step = 0.1)
                                                           )),
                                                
                                                
                                              column(9,
                                                     column(12, br()),
                                                     uiOutput("download_VolcanoPlots"),
                                                     column(12, br()),
                                                     plotOutput("volcano_plot", width = "100%", height = 650))
                                              ),
                                              box(title = "Heatmap", width = "100%", status = "primary",solidHeader = T, collapsible = T, 
                                                
                                                      column(3, sidebarPanel(width = "100%",
                                                                             sliderInput("heatmapsignificance", width = 300,"Change the significance level",min = 0.01,max = 0.1,value = 0.05,step = 0.01),
                                                                             sliderInput("heatmaplogfcup", width = 300, "Change the log2FC (Up)",min = 0, max = 10, value = 1,step = 0.1),
                                                                             sliderInput("heatmaplogfcdown", width = 300,"Change the log2FC (Down)", min = -10,max = 0,value = -1,step = 0.1),
                                                                             sliderInput("heatmap_degsnumber", width = 300,"Number of DEGs to include (each)", min = 1,max = 100,value = 50,step = 1),
                                                                             sliderInput("clusterchiecheatmap", width = 300,"Which cluster to plot", min = 0,max = 10,value = 2,step = 1),
                                                                             sliderInput(inputId = "heatmapgenelabesize", label = "Label size",min = 0, max = 20, value = 7, step = 1),
                                                                             actionButton(inputId = "heatmapaction", label = "Display map", width = "100%"))),
                                                      column(9,
                                                             downloadButton("downloadHeatmap",label = "Save plot"),
                                                             plotOutput("heatmap_analysis", width = "100%", height = 730)))),
                              ############# , width = "100%", height = 730
                                              )),
                              tabPanel(title = "Network", icon = icon("project-diagram"),
                                       box(width = "100%", status = "primary",solidHeader = T, collapsible = T, 
                                           box(title = "Network configuration", width = "100%", status = "primary",solidHeader = T, collapsible = T, 
                                               column(3,
                                                      sliderInput("genenetworksignificance", width = 300,"Change the significance level",min = 0.01,max = 0.1,value = 0.05,step = 0.01),
                                                      sliderInput("genenetworklogfcup", width = 300, "Change the log2FC (Up)",min = 0, max = 10, value = 1,step = 0.1),
                                                      sliderInput("genenetworklogfcdown", width = 300,"Change the log2FC (Down)", min = -10,max = 0,value = -1,step = 0.1)),
                                               column(3,
                                                      sliderInput("genenetwork_degsnumber", width = 300,"Number of DEGs to include (each)", min = 1,max = 1000,value = 150,step = 1),
                                                      sliderInput("clustergenenetwork", width = 300,"Which cluster to plot", min = 0,max = 10,value = 2,step = 1),
                                                      selectInput("genenetworkcormethod",label = "Correlation method",choices = c("pearson","spearman"), multiple = F, width = "100%")),
                                               
                                               column(3,
                                                      sliderInput(inputId = "genenetworkcordegree", label = "Remove edges below correlation",min = 0, max = 0.99, value = 0.6, step = 0.01),
                                                      selectInput("genenetworkmode", label = "Network mode",choices = c("directed", "undirected", "max", "min", "upper", "lower", "plus"), multiple = F, width = "100%"),
                                                      sliderInput(inputId = "genenetworklabledist", label = "Lables distance",min = 0.1, max = 1.5, value = 0.2, step = 0.1)),
                                               column(3,
                                                      selectInput(inputId = "networklayout", label = "Layout", choices = c("layout.auto", "layout.random" ,"layout.circle","layout.sphere", "layout.fruchterman.reingold", "layout.kamada.kawai", "layout.spring", "layout.reingold.tilford","layout.fruchterman.reingold.grid", "layout.lgl", "layout.graphopt", "layout.mds", "layout.svd", "layout.norm"), multiple = F, selected = "layout.random"),
                                                      sliderInput(inputId = "genenetworklablesize", label = "Size of lables",min = 0.1, max = 1.5, value = 1, step = 0.1),
                                                      br(),
                                                      br(),
                                                      actionButton(inputId = "genenetworkeaction", label = "Display network", width = "100%")),
                                           ),
                                           column(12,
                                                  downloadButton("downloadNetwork",label = "Save plot"),
                                                  plotOutput("genenetwork_analysis", width = "100%", height = 850)),
                                           
                                           
                                       )),
                                 
                                 tabPanel(title = "Pathway", icon = icon("table"),
                                          box(width = "100%", status = "primary",solidHeader = T, collapsible = T,
                                              column(3),
                                              column(3, sliderInput("pathwaysignificance", width = 300, "Significance level",min = 0.01, max = 0.1, value = 0.05,step = 0.01)),
                                              column(3, sliderInput("numberoftoppathways", width = 300, "Top pathways:",min = 1, max = 50, value = 20,step = 1)),
                                              column(3),
                                              column(5), column(3, actionButton(inputId = "pathwayaction", label = "Find pathways", width = "100%")),column(5),
                                              column(12, br()),
                                              column(12, br()),
                                              column(12, br()),
                                              conditionalPanel(condition = "input.pathwayaction",
                                              column(12,box(title = "KEGG Pathways", width = "100%", status = "primary",solidHeader = T, collapsible = T,
                                                            column(12, 
                                                                   uiOutput("download_KEGGPathway"),
                                                                   plotOutput("kegg_up_down_pathways"))),
                                                     box(title = "Up regulated pathways", width = "100%",height = 1200, status = "primary",solidHeader = T, collapsible = T,
                                                         column(12, br()),
                                                         column(12, column(5), column(2, selectInput(inputId ="kegg_up_path_graphs" , label = "Select a KEGG ID", choices = "", multiple = F),
                                                                                      actionButton(inputId = "kegg_up_path_graphsaction", label = "Find KEGG graphs", width = "100%")), column(5)),
                                                         column(12,plotOutput("downloading_up_kegg"))),
                                                     
                                              box(title = "Down regulated pathways", width = "100%", height = 1200, status = "primary",solidHeader = T, collapsible = T,
                                                            column(12, column(5), column(2, selectInput(inputId ="kegg_down_path_graphs" , label = "Select a KEGG ID", choices = NULL, multiple = F),
                                                                   actionButton(inputId = "kegg_down_path_graphsaction", label = "Find KEGG graphs", width = "100%")), column(5)),
                                                            column(12, plotOutput("downloading_down_kegg"))
                                                            ))
                                          ))),
                                 tabPanel(title = "GSEA", icon = icon("table"),
                                          box(width = "100%", height = 1000, status = "primary",solidHeader = T, collapsible = T,
                                              column(3, selectInput("gseacategory", label = "Category",choices = c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"), multiple = F, width = "100%")),
                                              column(3, sliderInput("gseasignificance", width = 300, "Significance level",min = 0.01, max = 0.1, value = 0.05,step = 0.01)),
                                              column(3, sliderInput("gseanumber", width = 300, "Top pathways:",min = 1, max = 50, value = 20,step = 1)),
                                              column(3, 
                                                     br(),
                                                     br(),
                                                     actionButton(inputId = "gsea_action", label = "GSEA", width = "100%")),
                                              column(12,     uiOutput("download_GSEA")), 
                                              br(),
                                                               plotOutput("GSEA_Plot", width = "100%", height = 700) #
                                                  
                                              ))
                     )),
  )
)
)











dashboardPage(
  header = dashheader,
  sidebar = dashSidebar,
  body = dashBody,
  title = "scRNA-Seq Analyzer",
  skin = "purple"
  
)