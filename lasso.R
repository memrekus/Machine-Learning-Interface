#LASSO


#load libraries

library(shiny)
library(dplyr)
library(ggplot2)
library(shinydashboard)
library(DT)
library(glmnet)

#load datas

miRNA<-readRDS("Normal_miRNA.Rds")[2:1882]

#RNA<-readRDS("Normal_RNA.Rds")[2:60661]

# ui

lassoUI<-function(id, label, choices) {
  ns <- NS(id)
  tagList(plotOutput(outputId = ns("lassooutput")),
          h3("The Most Important Predictor Genes"),
          dataTableOutput(outputId = ns("lassoDToutput")),
          plotlyOutput(outputId = ns("lollipop2")),
          plotlyOutput (outputId = ns("scatterall2")),
          plotlyOutput(outputId = ns("top202"))
          
  )}

#server

lassoSERVER<- function(id) {moduleServer(
  id,function(input, output, session) {
    ns <- session$ns
      output$lassooutput <- renderPlot({
        
        p<-miRNA %>% select(input$mirinput)
        pp<-as.numeric(unlist(p))
        varmtx<-miRNA %>% select(!(input$mirinput))
        varmtx<-as.matrix(varmtx)
        response <- pp
        
        lasso <- glmnet(scale(varmtx), response, alpha=1)
        
        ridge <- glmnet(scale(varmtx), response, alpha=0)
        
        # Cross validation to find the optimal lambda penalization
        cv.lasso <- cv.glmnet(varmtx, response, alpha=1)
        
        lbs_fun <- function(fit, offset_x=1, ...) {
          L <- length(fit$lambda)
          x <- log(fit$lambda[L])+ offset_x
          y <- fit$beta[, L]
          labs <- names(y)
          text(x, y, labels=labs, ...)
          abline(v=cv.lasso$lambda.min, col = "red", lty=2)
          abline(v=cv.lasso$lambda.1se, col="blue", lty=2)
        }
        
        plot(lasso, xvar = "lambda", label=T)
        lbs_fun(ridge, offset_x = -2)
     
      })
      
      output$lassoDToutput <- renderDataTable({
        
        p<-miRNA %>% select(input$mirinput)
        pp<-as.numeric(unlist(p))
        varmts<-miRNA %>% select(!(input$mirinput))
        varmtx<-as.matrix(varmts)
        response <- pp

        lasso <- glmnet(scale(varmtx), response, alpha=1)
        cv.lasso <- cv.glmnet(varmtx, response, alpha=1)
      
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
      
      output$lollipop2 <- renderPlotly({
        p<-miRNA %>% select(input$mirinput)
        p<-as.numeric(unlist(p))
        varmtx<-miRNA %>% select(!(input$mirinput))
        varmtx<-as.matrix(varmtx)
        response <- p
        
        lasso <- glmnet(scale(varmtx), response, alpha=1)
        cv.lasso <- cv.glmnet(varmtx, response, alpha=1)
        
        best_model <- glmnet(scale(varmtx), response, alpha = 0, lambda = cv.lasso$lambda.min)
        best_model<-coef(best_model)
        
        best_model <- as.matrix(best_model)
        best_model<-as.data.frame(best_model)
        
        best_model <- subset(best_model, rownames(best_model) != "(Intercept)")
        
        best_model <- best_model %>%
          arrange(s0) %>%
          slice(c(head(row_number(), 20), tail(row_number(), 20)))
        
        best_model <- data.frame(h = rownames(best_model), v = best_model[,1])
        best_model$h<-as.factor(best_model$h)
        
        lol2<-best_model%>%
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
        
        ggplotly(lol2)
        
      })
      
      output$scatterall2 <- renderPlotly({
        p<-miRNA %>% select(input$mirinput)
        p<-as.numeric(unlist(p))
        varmtx<-miRNA %>% select(!(input$mirinput))
        varmtx<-as.matrix(varmtx)
        response <- p
        
        lasso <- glmnet(scale(varmtx), response, alpha=1)
        cv.lasso <- cv.glmnet(varmtx, response, alpha=1)
        
        coef_lasso<- glmnet(scale(varmtx), response, alpha = 0, lambda = cv.lasso$lambda.min)
        coef_lasso<-coef(coef_lasso)
        coef_lasso<-as.data.frame(as.matrix(coef_lasso))
        coef_lasso <- subset(coef_lasso, rownames(coef_lasso) != "(Intercept)")
        colnames(coef_lasso)[1]<-"value"
        coef_lasso$miRNA<-rownames(coef_lasso)
        coef_lasso<-coef_lasso%>%arrange(desc(value))
        
        plot_ly(coef_lasso, x = ~miRNA, y = ~value, name = 'Predictor miRNAs', mode = 'markers',type ='scatter',color = "Orange")
      })
      
      output$top202 <- renderPlotly({
        
        p<-miRNA %>% select(input$mirinput)
        p<-as.numeric(unlist(p))
        varmtx<-miRNA %>% select(!(input$mirinput))
        varmtx<-as.matrix(varmtx)
        response <- p
        
        lasso <- glmnet(scale(varmtx), response, alpha=1)
        cv.lasso <- cv.glmnet(varmtx, response, alpha=1)
        
        best_model <- glmnet(scale(varmtx), response, alpha = 0, lambda = cv.lasso$lambda.min)
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
      
      })}