#ElasticNET


#load libraries

library(shiny)
library(dplyr)
library(ggplot2)
library(shinydashboard)
library(DT)
library(glmnet)
library(caret)
library(glmnetUtils)
#load datas

miRNA<-readRDS("Normal_miRNA.Rds")[2:1882]

#RNA<-readRDS("Normal_RNA.Rds")[2:60661]

#ui

elasticnetUI<-function(id, label, choices) {
  ns <- NS(id)
  tagList(h3("ElasticNET"),
          h4("The Most Important Predictor Genes"),
          #plotOutput(outputId = ns("rmse")),
          dataTableOutput(outputId = ns("elasticDToutput")),
          plotlyOutput(outputId = ns("lollipop3")),
          plotlyOutput (outputId = ns("scatterall3")),
          plotlyOutput(outputId = ns("top203"))
)}


#server

elasticnetSERVER<-function(id) {moduleServer(
  id,function(input, output, session) {
    
    # output$rmse <- renderPlot({
    #   
    #   p<-miRNA %>% select(input$mirinput)
    #   p<-as.numeric(unlist(p))
    #   response <- p
    #   rm(p)
    #   a<-miRNA$response
    #   
    #   miRNA_elnet <- train(
    #     hsa.mir.142 ~., data = miRNA, method = "glmnet",
    #     trControl = trainControl("cv", number = 10),
    #     tuneLength = 10)
    #   
    # })
    
    output$elasticDToutput <- renderDataTable({
      
      p<-miRNA %>% select(input$mirinput)
      p<-as.numeric(unlist(p))
      response <- p
      rm(p)
      varmtx<-miRNA %>% select(!(input$mirinput))
      varmtx<-as.matrix(varmtx)
      
      miRNA_elnet <- train(
        as.formula(paste(input$mirinput,"~.")) , data = miRNA, method = "glmnet",
        trControl = trainControl("cv", number = 10),
        tuneLength = 10)
      
      # miRNA_elnet <- train(
      #   hsa.mir.142 ~., data = miRNA, method = "glmnet",
      #   trControl = trainControl("cv", number = 10),
      #   tuneLength = 10)
      
      
      alpha2<-miRNA_elnet$bestTune$alpha
      lambda2<-miRNA_elnet$bestTune$lambda
      elnet<-glmnet(scale(varmtx),response,alpha = alpha2,lambda = lambda2)
  
      el_coef<-coef(elnet)
      el_coef <- as.matrix(el_coef)
      el_coef<-as.data.frame(el_coef)
      el_coef <- subset(el_coef, rownames(el_coef) != "(Intercept)")
      
      el_coef <- el_coef %>% 
        filter(!is.na(el_coef)) %>% 
        filter_at(vars(starts_with("s")), all_vars(. != 0))
      
      colnames(el_coef)[1]<-"Coefficient_Value"
      
      el_coef<-el_coef%>%arrange(desc(Coefficient_Value))
      
      print(el_coef)
  
       }) 
    output$lollipop3 <- renderPlotly({
      p<-miRNA %>% select(input$mirinput)
      p<-as.numeric(unlist(p))
      response <- p
      rm(p)
      varmtx<-miRNA %>% select(!(input$mirinput))
      varmtx<-as.matrix(varmtx)
      
      miRNA_elnet <- train(
        as.formula(paste(input$mirinput,"~.")), data = miRNA, method = "glmnet",
        trControl = trainControl("cv", number = 10),
        tuneLength = 10)
      
      alpha2<-miRNA_elnet$bestTune$alpha
      lambda2<-miRNA_elnet$bestTune$lambda
      elnet<-glmnet(scale(varmtx),response,alpha = alpha2,lambda = lambda2)
      
      el_coef<-coef(elnet)
      el_coef <- as.matrix(el_coef)
      el_coef<-as.data.frame(el_coef)
      el_coef <- subset(el_coef, rownames(el_coef) != "(Intercept)")
      
      el_coef <- el_coef %>% 
        filter(!is.na(el_coef)) %>% 
        filter_at(vars(starts_with("s")), all_vars(. != 0))
      
      el_coef<-el_coef%>%arrange(desc(s0))
      
      el_coef <- el_coef %>%
        arrange(s0) %>%
        slice(c(head(row_number(), 20), tail(row_number(), 20)))
      
      el_coef <- data.frame(h = rownames(el_coef), v = el_coef[,1])
      el_coef$h<-as.factor(el_coef$h)
      
      lol<-el_coef%>%
        arrange(v)%>%
        ggplot(aes(x=reorder(h,-v), y=v)) +
        geom_segment( aes(xend=h, yend=0), color="grey") +
        geom_point( color="orange", size=4) +
        theme_light() +
        theme(
          axis.text.x = element_text(angle =90 ),
          panel.grid.major.x = element_blank(),
          panel.border = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
          axis.title.x = element_text(vjust=-2),
          axis.title.y = element_text(angle=90, vjust=3)
        ) +
        xlab("Explanatory variable") +
        ylab("Regression coefficient")
      
      ggplotly(lol)
      
    })
    
    output$scatterall3 <- renderPlotly({
      p<-miRNA %>% select(input$mirinput)
      p<-as.numeric(unlist(p))
      response <- p
      rm(p)
      varmtx<-miRNA %>% select(!(input$mirinput))
      varmtx<-as.matrix(varmtx)
      
      miRNA_elnet <- train(
        as.formula(paste(input$mirinput,"~.")), data = miRNA, method = "glmnet",
        trControl = trainControl("cv", number = 10),
        tuneLength = 10)
      
      alpha2<-miRNA_elnet$bestTune$alpha
      lambda2<-miRNA_elnet$bestTune$lambda
      elnet<-glmnet(scale(varmtx),response,alpha = alpha2,lambda = lambda2)
      
      el_coef<-coef(elnet)
      el_coef <- as.matrix(el_coef)
      el_coef<-as.data.frame(el_coef)
      el_coef <- subset(el_coef, rownames(el_coef) != "(Intercept)")
      
      el_coef <- el_coef %>% 
        filter(!is.na(el_coef)) %>% 
        filter_at(vars(starts_with("s")), all_vars(. != 0))
      
      el_coef<-el_coef%>%arrange(desc(s0))
      
      colnames(el_coef)[1]<-"value"
      el_coef$miRNA<-rownames(el_coef)

      plot_ly(el_coef, x = ~miRNA, y = ~value, name = 'Predictor miRNAs', mode = 'markers',type ='scatter',color = "Orange")
    })
    
    output$top203 <- renderPlotly({
      p<-miRNA %>% select(input$mirinput)
      p<-as.numeric(unlist(p))
      response <- p
      rm(p)
      varmtx<-miRNA %>% select(!(input$mirinput))
      varmtx<-as.matrix(varmtx)
      
      miRNA_elnet <- train(
        as.formula(paste(input$mirinput,"~.")), data = miRNA, method = "glmnet",
        trControl = trainControl("cv", number = 10),
        tuneLength = 10)
      
      alpha2<-miRNA_elnet$bestTune$alpha
      lambda2<-miRNA_elnet$bestTune$lambda
      elnet<-glmnet(scale(varmtx),response,alpha = alpha2,lambda = lambda2)
      
      el_coef<-coef(elnet)
      el_coef <- as.matrix(el_coef)
      el_coef<-as.data.frame(el_coef)
      el_coef <- subset(el_coef, rownames(el_coef) != "(Intercept)")
      
      el_coef <- el_coef %>% 
        filter(!is.na(el_coef)) %>% 
        filter_at(vars(starts_with("s")), all_vars(. != 0))
      
      el_coef<-el_coef%>%arrange(desc(s0))
      
      el_coef <- el_coef %>%
        arrange(s0) %>%
        slice(c(head(row_number(), 20), tail(row_number(), 20)))
      
      el_coef <- data.frame(h = rownames(el_coef), v = el_coef[,1])
      el_coef$h<-as.factor(el_coef$h)
      
      x<-ggplot(el_coef,aes(v, reorder(h, v), color = v > 0)) +
        geom_point(show.legend = FALSE) +
        ggtitle("Influential variables") +
        xlab("Coefficient") +
        ylab(NULL)
      ggplotly(x)
    })
  
    })}