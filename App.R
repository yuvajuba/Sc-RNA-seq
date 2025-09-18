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



# User Interface    ------------------------------------------------------------
ui <- fluidPage(
  includeCSS("www/ShinyStyle.css"),
  includeHTML("ShinyHeader.html"),
  
  ## --------------------------------------------------------------------------- I- Import Data  ----
  h1("I- Importing the data", class = "body", style = "margin-top:50px"),
  fluidRow(
    column(
      width = 5,
      fileInput(inputId = "LoadSeurat",
                label = "Load your SeuratObject",
                accept = ".rds",
                width = "90%"),
      verbatimTextOutput(outputId = "InfoSeurat", placeholder = F),
      p("Quick overview :", class = "sm-title1", style = "margin-top:50px"),
      selectInput(inputId = "explore_plot_select",
                  label = "Select meta to plot",
                  choices = c(),
                  width = "50%"),
      radioButtons(inputId = "split_plot",
                   label = "Split the plot",
                   choices = c("Yes", "No"),
                   inline = T),
      selectInput(inputId = "split_by",
                  label = "Split by",
                  choices = c(),
                  width = "50%"),
      actionButton("explore_meta", "View", class = "btn-perso2")
    ),
    
    column(
      width = 7,
      div(
        class = "box1-for-plot",
        plotOutput("plt1", width = "auto", height = "auto")
      )
    )
  ),
  
  
  ## --------------------------------------------------------------------------- II- DEA ----
  h1("II- Differential Expression Analysis (DEA)", class = "body"),
  
  ### -------------------------------------------------------------------------- 1- Select contrasts ----
  h1("1- Contrasts selection", class = "body2", style = "margin-top:0px"),
  p("In this part, you are going to select the 2 groups you want to have the differentially 
    expressed genes (DEGs) from.", style = "margin-bottom:30px"),
  fluidRow(
    column(
      width = 5,
      selectInput(inputId = "Select_meta",
                  label = "Select your meta",
                  choices = c(),
                  width = "50%"),
      selectInput(inputId = "Select_group1",
                  label = "Select group 1",
                  choices = c(),
                  width = "50%"),
      selectInput(inputId = "Select_group2",
                  label = "Select group 2",
                  choices = c(),
                  width = "50%"),
      div(
        style = "margin-top:40px ; margin-bottom:20px",
        actionLink(inputId = "summUp",
                   label = "View cell counts",
                   icon = icon("hand-point-right"),
                   style = "font-size:18px ; color:purple ; font-weight:600")
      ),
      tableOutput("SummaryTable"),
      actionButton("ViewContrasts", 
                   "View contrats", 
                   class = "btn-perso2",
                   style = "margin-top:50px")
    ),
    
    column(
      width = 7,
      div(
        class = "box1-for-plot",
        plotOutput("plt2", width = "auto", height = "auto")
      )
    )
  ),
  
  
  ### -------------------------------------------------------------------------- 2- DE Analysis ----
  hr(),
  h1("2- DE Analysis", class = "body2"),
  p("In this part, you are going to define parameters for the DEA such as the different thresholds 
    and then running the analysis.", 
    style = "margin-bottom:30px"),
  fluidRow(
    
    #### ----------------------------------------------------------------------- 2-1 Setup params ----
    column(
      width = 3,
      p("Setup DEA params", class = "sm-title1"),
      numericInput("log2fc",
                   "Log2FC threshold",
                   value = 0.5, min = 0, max = 3, step = 0.1, width = "50%"),
      numericInput("min.pct",
                   "Cell pct threshold",
                   value = 0.3, min = 0, max = 1, step = 0.05, width = "50%"),
      numericInput("min.diff.pct",
                   "Cell diff.pct",
                   value = 0, min = 0, max = 1, step = 0.05, width = "50%"),
      
      p("Additional filtering params", class = "sm-title1"),
      sliderInput("log2fc_up",
                  "Log2FC up.thresh",
                  value = c(0.5,10), min = 0, max = 10, step = 0.25, width = "100%"),
      sliderInput("log2fc_down",
                  "Log2FC down.thresh",
                  value = c(-10,-0.5), min = -10, max = 0, step = 0.25, width = "100%"),
      
      actionButton("runDEA", 
                   "RUN DEA", 
                   class = "btn-perso2",
                   style = "margin-top:50px")
    ),
    
    #### ----------------------------------------------------------------------- 2-1 RUN DEA ----
    column(
      width = 9
    )
  )
)






# Shiny Function    ------------------------------------------------------------
server <- function(session, input, output){
  
  # ---------------------------------------------------------------------------- I- Import Data ----
  seurat <- eventReactive(input$LoadSeurat, {
    req(input$LoadSeurat)
    readRDS(input$LoadSeurat$datapath)
  })
  
  output$InfoSeurat <- renderPrint({
    req(seurat())
    print(seurat())
  })
  
  observeEvent(seurat(), {
    req(seurat())
    seurat <- seurat()
    meta <- names(seurat@meta.data)[sapply(seurat@meta.data, function(x) !is.numeric(x))]
    
    updateSelectInput(session, "explore_plot_select", choices = meta)
    updateSelectInput(session, "split_by", choices = meta)
    updateSelectInput(session, "Select_meta", choices = meta)
  })
  
  plot1 <- eventReactive(input$explore_meta, {
    req(seurat())
    seurat <- seurat()
    val <- length(unique(seurat@meta.data[[input$explore_plot_select]]))
    
    if(input$split_plot == "Yes"){
      DimPlot(seurat,
              group.by = input$explore_plot_select,
              split.by = input$split_by,
              pt.size = 1.2,
              alpha = 0.6,
              cols = MyPalette)+
        My_umap_shiny()+
        theme(plot.title = element_text(margin = margin(b=0.1, unit = "in")))+
        guides(colour = guide_legend(override.aes = list(size = 4),
                                     ncol = ifelse(val <= 36, 1, 2)))
    } else {
      DimPlot(seurat,
              group.by = input$explore_plot_select,
              pt.size = 1.2,
              alpha = 0.6,
              cols = MyPalette)+
        My_umap_shiny()+
        guides(colour = guide_legend(override.aes = list(size = 4),
                                     ncol = ifelse(val <= 36, 1, 2)))
    }
    
    
  })
  
  output$plt1 <- renderPlot({
    req(plot1())
    print(plot1())
  }, width = "auto", height = 600, res = 85)
  
  
  
  # ---------------------------------------------------------------------------- II- DEA ----
  ## --------------------------------------------------------------------------- 1- Select contrasts ----
  observeEvent(input$Select_meta, {
    req(seurat(), input$Select_meta)
    seurat <- seurat()
    grp <- unique(seurat@meta.data[[input$Select_meta]])
    
    updateSelectInput(session, "Select_group1", choices = grp)
    updateSelectInput(session, "Select_group2", choices = grp)
  })
  
  summary_table <- eventReactive(input$summUp, {
    req(seurat(), input$Select_group1, input$Select_group2)
    seurat <- seurat()
    grp <- c(input$Select_group1, input$Select_group2)
    
    seurat@meta.data %>% 
      count(!!sym(input$Select_meta)) %>% 
      filter(!!sym(input$Select_meta) %in% grp)
  })
  
  output$SummaryTable <- renderTable({
    req(summary_table())
    summary_table()
  })
  
  plot2 <- eventReactive(input$ViewContrasts, {
    req(seurat(), input$Select_group1, input$Select_group2)
    seurat <- seurat()
    grp <- c(input$Select_group1, input$Select_group2)
    meta_column <- input$Select_meta
    
    labeled_umap(seuratobject = seurat,
                 metadata = meta_column,
                 cluster = grp,
                 palette = MyPalette)
  })
  
  output$plt2 <- renderPlot({
    req(plot2())
    print(plot2())
  }, width = "auto", height = 600, res = 85)
  
  
  
}







# Run Shiny     ----------------------------------------------------------------
shinyApp(ui, server)


































