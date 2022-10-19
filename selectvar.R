#SELECTVAR


#load libraries

library(shiny)
library(dplyr)
library(ggplot2)
library(shinydashboard)
library(DT)
library(fontawesome)
library(lattice)
library(ggcorrplot)
library(plotly)
library(remotes)
library(elasticnet)

#load datas

miRNA<-readRDS("Normal_miRNA.Rds")[2:1882]

#RNA<-readRDS("Normal_RNA.Rds")[2:60661]

#ui

selectvar<-function(id, label, choices) {
  ns <- NS(id)
  tagList(selectInput( inputId =  ns("mirinput"),
                       label = "Choose the gene...",
                       choices = colnames(miRNA)
                       ),
          # selectInput(inputId =  ns("mlmodel"),
          #             label = "Choose the machine learning algorithm",
          #             choices = c("LASSO","RIDGE REGRESSION","ELASTICNET")
          # ),
    submitButton("SUBMIT",icon = icon("redo")),
    
    # inputId =  ns("orp"),
    #                   label = "Choose RNA encoding genes or miRNA encoding genes...",
    #                   choices = "miRNA","RNA"),
    #      
        textOutput(outputId = ns("miroutput"))
    )}
           

#server

selectvarSERVER<-function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      
      output$miroutput<-renderText (
      
      paste("Chosen gene ",input$mirinput, sep = ":")
      
      )})}


        
#       output$mir<-renderDataTable({
#       # data<-switch(input$orp,
#       #              "miRNA"=miRNA,
#       #              "RNA"=RNA)
#       #  head(data) 
#       miRNA %>% select(input$mir)
#         # x<-(miRNA[,colnames(miRNA)==hsa.let.7a.1])
#         # print(x)
#         })})}