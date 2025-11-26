# Clean start :
# ----------- :
rm(list = ls())
graphics.off()
cat("\014")

# SETUP ENVIRONNEMENT :
# =================== :
source("00_setup.R")

# Shiny App maxsize file loading :
# ============================== :
options(shiny.maxRequestSize = 900*1024^2)

# UI :    ----------------------------------------------------------------------

ui <- fluidPage(
  title = "Shiny-App scRNAseq",
  includeCSS("www/tmp-style.css"),
  includeHTML("tmp-shiny-header.html"),
  
  h1("I- Importing the data", class = "section", id = "I"),
  
  fluidRow(
    column(
      width = 12, 
      class = "import-data-container",
      
      div(
        class = "import-file-container",
        p("Importing the data", class = "sm-title green"),
        fileInput(inputId = "loadseurat",
                  label = "Load your SeuratObject",
                  accept = ".rds"),
        verbatimTextOutput(outputId = "infoseurat", 
                           placeholder = F)
      ),
      
      div(
        class = "plt1-params",
        p("Exploring the data", class = "sm-title green"),
        selectInput(inputId = "plt1_select",
                    label = "Select meta to plot",
                    choices = c()),
        radioButtons(inputId = "plt1_split",
                     label = "Split the plot",
                     choices = c("Yes", "No"),
                     inline = T),
        selectInput(inputId = "plt1_split_by",
                    label = "Split by",
                    choices = c()),
        actionButton(inputId = "plt1_view", 
                     label = "View", 
                     class = "btn btn-1")
      ),
      
      div(
        class = "plt1-container",
        plotOutput("plt1", width = "auto", height = "auto")
      )
    )
  ),
  
  
  
  h1("II- Differentially expressed analysis (DEA)", class = "section", id = "II"),
  
  fluidRow(
    column(width = 12, 
           class = "c1")
  )
)



server <- function(session, input, output){
  # ---------------------------------------------------------------------------- I- Import Data ----
  
  ## --- Import and load seurat object:
  seurat <- eventReactive(input$loadseurat, {
    req(input$loadseurat)
    readRDS(input$loadseurat$datapath)
  })
  
  ## --- Print seurat object info:
  output$infoseurat <- renderPrint({
    req(seurat())
    print(seurat())
  })
  
  ## --- Update forms:
  observeEvent(seurat(), {
    req(seurat())
    seurat <- seurat()
    meta <- names(seurat@meta.data)[sapply(seurat@meta.data, function(x) !is.numeric(x))]
    
    updateSelectInput(session, "plt1_select", choices = meta)
    updateSelectInput(session, "plt1_split_by", choices = meta)
  })
  
  ## --- Create plot1:
  plot1 <- eventReactive(input$plt1_view, {
    req(seurat())
    seurat <- seurat()
    val <- length(unique(seurat@meta.data[[input$plt1_select]]))
    
    if(input$plt1_split == "Yes"){
      DimPlot(seurat,
              group.by = input$plt1_select,
              split.by = input$plt1_split_by,
              pt.size = 1.2,
              alpha = 0.6,
              cols = MyPalette)+
        My_umap_shiny()+
        theme(plot.title = element_text(margin = margin(b=0.1, unit = "in")))+
        guides(colour = guide_legend(override.aes = list(size = 4),
                                     ncol = ifelse(val <= 36, 1, 2)))
    } else {
      DimPlot(seurat,
              group.by = input$plt1_select,
              pt.size = 1.2,
              alpha = 0.6,
              cols = MyPalette)+
        My_umap_shiny()+
        guides(colour = guide_legend(override.aes = list(size = 4),
                                     ncol = ifelse(val <= 36, 1, 2)))
    }
  })
  
  ## --- Plot the output:
  output$plt1 <- renderPlot({
    req(plot1())
    print(plot1())
  }, width = "auto", height = 600)
}



shinyApp(ui, server)