## ShinyMeth

library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(FlowSorted.Blood.450k)
library(FlowSorted.CordBlood.450k)
library(FlowSorted.DLPFC.450k)
library(RColorBrewer)
library(wateRmelon)
library(car)
library(sva)
library(limma)

library(shiny)
library(shinyFiles)
library(shinysky)
library(shinythemes)

ui <- tagList(tags$header(
  tags$style(HTML("
                  @import url('//fonts.googleapis.com/css?family=Cherry+Swash');
                  
                  header {
                  font-family: 'Cherry Swash', cursive;
                  font-weight: 500;
                  line-height: 1.1;
                  }
                  
                  "))),
  navbarPage(
    title = tags$header(
      HTML('<a style="color: 	#FFA500" target  = "_blank" href="https://jangehlen.github.io/shinyOmics/">shinyOmics - Methylation</a>')),
    
    windowTitle = "ShinyOmics - Meth",
    id = "meth",
    theme = shinytheme("paper"),
    
    
    tabPanel("Input Section", icon = icon("cloud-upload", lib = "font-awesome"),
      sidebarLayout(
        sidebarPanel(
          fileInput("idat", "Select Idat File in Project directory"),
          fileInput("samplesheet", "Upload Samplesheet"),
          radioButtons("array", "Select Methylation Array", choices = c("Illumina 450K", "Illumina EPIC")),
          htmlOutput("sample_id_slc"),
          br(),
          br(),
          actionButton("idat_upload", "Start data import !", "warning")
        ),
      mainPanel(
        dataTableOutput("samplesheet"),
        verbatimTextOutput("test")
      ) 
            
      )
    )
  )
)

server <- function(input, output, session) {
  
  annotation <- reactive({
   array <- switch(input$array,
      "Illumina 450K" = IlluminaHumanMethylation450kanno.ilmn12.hg19,
      "Illumina EPIC" = IlluminaHumanMethylationEPICanno.ilm10b2.hg19
    )
   getAnnotation(array)
  })
  
  targets <- reactive({
    read.metharray.sheet(dirname(input$idat), basename(input$samplesheet))
  })

  output$sample_id_slc <- renderUI({
    req(targets())
    selectInput("s_ids", "Select Column with Sample IDs", choices = colnames(targets()))
  })
  
  pre_rgSet <- eventReactive(
    {
      input$idat_upload
      },
    {
      read.metharray.exp(targets=targets())
      })
  
  rgSet <- reactive({
    req(pre_rgSet())
    sampleNames(pre_rgSet) <- targets[[input$s_ids]]
  })
  
  norm_objects <- reactive({
    norm_raw  <- preprocessRaw(rgSet())
    norm_illumina <- preprocessIllumina(rgSet())
    norm_swan <- preprocessSWAN(rgSet())
    norm_noob <- preprocessNoob(rgSet())
    norm_quantile <- preprocessQuantile(rgSet())
    norm_funnorm  <- preprocessFunnorm(rgSet())
    
    norm_objects <- list(raw=norm_raw, 
                         illumina=norm_illumina,
                         swan=norm_swan, 
                         noob=norm_noob,
                         quantile=norm_quantile, 
                         funnorm=norm_funnorm)
    return(norm_objects)
  })
  
}
options(shiny.trace = T, shiny.launch.browser = T)
shinyApp(ui, server)
