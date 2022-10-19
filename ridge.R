#RIDGE


#load libraries

library(shiny)
library(dplyr)
library(ggplot2)
library(shinydashboard)
library(DT)
library(glmnet)
library(caret)
library(MXM)
library(ggcorrplot)
library(broom)

#load datas

miRNA<-readRDS("Normal_miRNA.Rds")[2:1882]

#RNA<-readRDS("Normal_RNA.Rds")[2:60661]

#ui

ridgeUI<-function(id, label, choices) {
  ns <- NS(id)
  tagList(plotOutput(outputId = ns("ridgeoutput")),
          h4("The Most Important Predictor Genes"),
          dataTableOutput(outputId = ns("ridgeDToutput")),
          plotlyOutput(outputId = ns("lollipop")),
          plotlyOutput (outputId = ns("scatterall")),
          plotlyOutput(outputId = ns("top20")),
          downloadButton(NS(id, "download"), label = "Download Graph"))
          
          
}

#server

ridgeSERVER<- function(id) {moduleServer(
  id,function(input, output, session) {
    ns <- session$ns
    
      output$ridgeoutput <- renderPlot({
        p<-miRNA %>% select(input$mirinput)
        p<-as.numeric(unlist(p))
        varmtx<-miRNA %>% select(!(input$mirinput))
        varmtx<-as.matrix(varmtx)
        response <- p
        
        ridge <- glmnet(scale(varmtx), response, alpha=0)
        # Cross validation to find the optimal lambda penalization
        cv.ridge <- cv.glmnet(varmtx, response, alpha=0)
        
        lbs_fun <- function(fit, offset_x=1, ...) {
          L <- length(fit$lambda)
          x <- log(fit$lambda[L])+ offset_x
          y <- fit$beta[, L]
          labs <- names(y)
          text(x, y, labels=labs, ...)
          abline(v=cv.ridge$lambda.min, col = "red", lty=2)
          abline(v=cv.ridge$lambda.1se, col="blue", lty=2)
        }
        plot(ridge, xvar = "lambda", label=T)
        lbs_fun(ridge)
      })
      
      output$ridgeDToutput <- renderDataTable({
        p<-miRNA %>% select(input$mirinput)
        p<-as.numeric(unlist(p))
        varmtx<-miRNA %>% select(!(input$mirinput))
        varmtx<-as.matrix(varmtx)
        response <- p
        
        ridge <- glmnet(scale(varmtx), response, alpha=0)
        # Cross validation to find the optimal lambda penalization
        cv.ridge <- cv.glmnet(varmtx, response, alpha=0)
        
          best_model <- glmnet(scale(varmtx), response, alpha = 0, lambda = cv.ridge$lambda.min)
          best_model<-coef(best_model)

          best_model <- as.matrix(best_model)
          best_model<-as.data.frame(best_model)

          best_model <- subset(best_model, rownames(best_model) != "(Intercept)")

          best_model <- best_model %>%
            arrange(s0) %>%
            slice(c(head(row_number(), 50), tail(row_number(), 50)))

          colnames(best_model)[1]<-"Coefficient_Value"

          best_model<-best_model%>%arrange(desc(Coefficient_Value))

          best_model
        })    
      
      output$lollipop <- renderPlotly({
        p<-miRNA %>% select(input$mirinput)
        p<-as.numeric(unlist(p))
        varmtx<-miRNA %>% select(!(input$mirinput))
        varmtx<-as.matrix(varmtx)
        response <- p
        
        ridge <- glmnet(scale(varmtx), response, alpha=0)
        # Cross validation to find the optimal lambda penalization
        cv.ridge <- cv.glmnet(varmtx, response, alpha=0)
        
        lbs_fun <- function(fit, offset_x=1, ...) {
          L <- length(fit$lambda)
          x <- log(fit$lambda[L])+ offset_x
          y <- fit$beta[, L]
          labs <- names(y)
          text(x, y, labels=labs, ...)
          abline(v=cv.ridge$lambda.min, col = "red", lty=2)
          abline(v=cv.ridge$lambda.1se, col="blue", lty=2)
        }
        
        best_model <- glmnet(scale(varmtx), response, alpha = 0, lambda = cv.ridge$lambda.min)
        best_model<-coef(best_model)
        
        best_model <- as.matrix(best_model)
        best_model<-as.data.frame(best_model)
        
        best_model <- subset(best_model, rownames(best_model) != "(Intercept)")
        
        best_model <- best_model %>%
          arrange(s0) %>%
          slice(c(head(row_number(), 20), tail(row_number(), 20)))
        
        best_model <- data.frame(h = rownames(best_model), v = best_model[,1])
        best_model$h<-as.factor(best_model$h)
        
        lol<-best_model%>%
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
      
      output$scatterall <- renderPlotly({
        p<-miRNA %>% select(input$mirinput)
        p<-as.numeric(unlist(p))
        varmtx<-miRNA %>% select(!(input$mirinput))
        varmtx<-as.matrix(varmtx)
        response <- p
        
        ridge <- glmnet(scale(varmtx), response, alpha=0)
        # Cross validation to find the optimal lambda penalization
        cv.ridge <- cv.glmnet(varmtx, response, alpha=0)
        
        lbs_fun <- function(fit, offset_x=1, ...) {
          L <- length(fit$lambda)
          x <- log(fit$lambda[L])+ offset_x
          y <- fit$beta[, L]
          labs <- names(y)
          text(x, y, labels=labs, ...)
          abline(v=cv.ridge$lambda.min, col = "red", lty=2)
          abline(v=cv.ridge$lambda.1se, col="blue", lty=2)
        }
        
        coef_ridge<- glmnet(scale(varmtx), response, alpha = 0, lambda = cv.ridge$lambda.min)
        coef_ridge<-coef(coef_ridge)
        coef_ridge<-as.data.frame(as.matrix(coef_ridge))
        coef_ridge <- subset(coef_ridge, rownames(coef_ridge) != "(Intercept)")
        colnames(coef_ridge)[1]<-"value"
        coef_ridge$miRNA<-rownames(coef_ridge)
        coef_ridge<-coef_ridge%>%arrange(desc(value))
        
        plot_ly(coef_ridge, x = ~miRNA, y = ~value, name = 'Predictor miRNAs', mode = 'markers',type ='scatter',color = "Orange")
      })
      
      output$top20 <- renderPlotly({
        p<-miRNA %>% select(input$mirinput)
        p<-as.numeric(unlist(p))
        varmtx<-miRNA %>% select(!(input$mirinput))
        varmtx<-as.matrix(varmtx)
        response <- p
        
        ridge <- glmnet(scale(varmtx), response, alpha=0)
        # Cross validation to find the optimal lambda penalization
        cv.ridge <- cv.glmnet(varmtx, response, alpha=0)
        
        best_model <- glmnet(scale(varmtx), response, alpha = 0, lambda = cv.ridge$lambda.min)
        best_model<-coef(best_model)
        
        best_model <- as.matrix(best_model)
        best_model<-as.data.frame(best_model)
        
        best_model <- subset(best_model, rownames(best_model) != "(Intercept)")
        
        best_model <- best_model %>%
          arrange(s0) %>%
          slice(c(head(row_number(), 20), tail(row_number(), 20)))
        
        best_model <- data.frame(h = rownames(best_model), v = best_model[,1])
        best_model$h<-as.factor(best_model$h)
        x<-ggplot(best_model,aes(v, reorder(h, v), color = v > 0)) +
          geom_point(show.legend = FALSE) +
          ggtitle("Influential variables") +
          xlab("Coefficient") +
          ylab(NULL)
        ggplotly(x)
      })
      
      output$download<-downloadHandler(filename = function(){
        paste0('zort.png')},
          content = function(file) {ggsave(file, output$top20())})
      
      
      
      
        # varimp$variableMiRNAS<-rownames(varimp)
        # varimp<-varimp[order(varimp$Overall, rev(varimp$variableMiRNAS), decreasing = TRUE), ]
        # varimp$variableMiRNAS<-NULL
        # varimp<-varimp%>%select() 
        
        })}
