#https://github.com/MrFuguDataScience/MLE_Basics_N_Stats/blob/master/Lasso_Ridge_exp_Rstudio.ipynb

#https://glmnet.stanford.edu/articles/glmnet.html

# http://www.sthda.com/english/articles/37-model-selection-essentials-in-r/153-penalized-regression-essentials-ridge-lasso-elastic-net/

# https://www.anycodings.com/1questions/5640334/how-to-supress-code-from-showing-up-in-shiny-ui

# https://www.datacamp.com/tutorial/tutorial-ridge-lasso-elastic-net

p<-miRNA %>% select(hsa.mir.142)
p<-as.numeric(unlist(p))
varmtx<-miRNA %>% select(!(hsa.mir.142))
varmtx<-as.matrix(varmtx)
response <- as.vector(p)

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

######

##top predictor genes LOLLIPOP

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

###############SCATTERALL
coef_ridge<- glmnet(scale(varmtx), response, alpha = 0, lambda = cv.ridge$lambda.min)
coef_ridge<-coef(coef_ridge)
coef_ridge<-as.data.frame(as.matrix(coef_ridge))
coef_ridge <- subset(coef_ridge, rownames(coef_ridge) != "(Intercept)")
colnames(coef_ridge)[1]<-"value"
coef_ridge$miRNA<-rownames(coef_ridge)
coef_ridge<-coef_ridge%>%arrange(desc(value))

plot_ly(coef_ridge, x = ~miRNA, y = ~value, name = 'Predictor miRNAs', mode = 'markers',type ='scatter',color = "Orange")
###############

################TOP20-+
x<-ggplot(best_model,aes(v, reorder(h, v), color = v > 0)) +
  geom_point(show.legend = FALSE) +
  ggtitle("Influential variables") +
  xlab("Coefficient") +
  ylab(NULL)
ggplotly(x)
######



set.seed(123)

n <- 1000
p <- 2000

pred <- matrix(rnorm(n*p), nrow = n, ncol = p)



dv <- (rowSums(pred[,1:5]) + .8*rowSums(pred[,6:10]) + 
         .6 * rowSums(pred[,11:15]) + .4*rowSums(pred[,16:20]) +
         .2 * rowSums(pred[,21:25]) + rnorm(n))


pred <- scale(pred)

train_rows <- sample(1:n, .7*n, replace = F)

pred.train <- pred[train_rows,]
dv.train <- dv[train_rows]

pred.test <- pred[-train_rows,]
dv.test <- dv[-train_rows]

####

dv2<- miRNA$hsa.mir.142


train_rows2 <- miRNA$hsa.mir.142

pred.train2 <- pred[train_rows2,]
dv.train2 <- dv[train_rows2]

pred.test2 <- pred[-train_rows2,]
dv.test2 <- dv[-train_rows2]


####################################

# Get coefficients of all 100 models
ridge_coef <- coef(ridge)

# Display coefficients for 6 models. 
# Ridge Models are displayed in decreasing order of lambda.
ridge_coef<-as.matrix(ridge_coef)

round(ridge_coef[c(1:3, 98:100),], 6)
plot(ridge_coef,xvar = "lambda", label=T)

#############

varimp<-varImp(ridge,lambda = cv.ridge$lambda.min)
varimp<-as.data.frame(varimp)
xd<-varimp %>% slice_max(Overall, n = 20)


model <- train(
  hsa.mir.142 ~ .,
  data = miRNA,
  method = 'ridge'
)
model

axd<-Vi[1:20,1:20]

plot(varImp(ridge,lambda = cv.ridge$lambda.min))
ggplot(axd)
plot(Vi, xvar = "lambda", label=T)

# dataset <- matrix(runif(300 * 20, 1, 20), nrow = 300 ) 
# #the target feature is the last column of the dataset as a vector
# target <- dataset[, 20]
# dataset <- dataset[, -20]
# ridge.plot(target, dataset)
##################

#######################################
data<-DirtyData

cc<-data[,!duplicated(col(data))]

ee <- data[!duplicated(as.list(data)) ]

bb<-data[!duplicated(unclass(data))]

df <- DirtyData %>%
  mutate(meta.ajcc_pathologic_m = coalesce(meta.ajcc_pathologic_m.y,meta.ajcc_pathologic_m.x)) %>%
  select(-c(meta.ajcc_pathologic_m.y,meta.ajcc_pathologic_m.x))

#######################################

# set.seed(42)
#   cv_5 = trainControl(method = "cv", number = 5)
#   
#   miRNA_elnet = train(
#     hsa.mir.142 ~ ., data = miRNA,
#     method = "glmnet",
#     trControl = cv_5)
#   
#   miRNA_elnet_int = train(
#     hsa.mir.142 ~ ., data = miRNA,
#     method = "glmnet",
#     trControl = cv_5,
#     tuneLength = 10)
#   
#   # miRNA_elnet_int = train(
#   #   hsa.mir.142 ~ . ^ 2, data = miRNA,
#   #   method = "glmnet",
#   #   trControl = cv_5,
#   #   tuneLength = 10)
#   # Error: cannot allocate vector of size 12.4 Gb
#   
#   get_best_result = function(caret_fit) {
#     best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
#     best_result = caret_fit$results[best, ]
#     rownames(best_result) = NULL
#     best_result
#   }
#   
#   get_best_result(miRNA_elnet)
#   get_best_result(miRNA_elnet_int)
# 
# X <- miRNA %>% 
#   select(hsa.let.7a.1) %>% 
#   scale(center = TRUE, scale = FALSE) %>% 
#   as.matrix()
# 
# Y <- miRNA %>% 
#   select(-hsa.let.7a.1) %>% 
#   as.matrix()
# 
# control <- trainControl(method = "repeatedcv",
#                         number = 5,
#                         repeats = 5,
#                         search = "random",
#                         verboseIter = TRUE)
# elastic_model <- train(hsa.let.7a.1 ~ .,
#                        data =cbind(X, Y),
#                        method = "glmnet",
#                        preProcess = c("center", "scale"),
#                        tuneLength = 25,
#                        trControl = control)
# # Model Prediction
# premodel <- predict(elastic_model,X)
# 
# # Multiple R-squared
# rsq <- cor(X, premodel)^2
# 
# # Plot
# plot(elastic_model)

##########################

p<-miRNA %>% select(hsa.mir.142)
p<-as.numeric(unlist(p))
varmtx<-miRNA %>% select(!(hsa.mir.142))
varmtx<-as.matrix(varmtx)
response <- p
rm(p)

miRNA_elnet <- train(
  hsa.mir.142 ~., data = miRNA, method = "glmnet",
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

el_coef



# #set a seed so you can reproduce the results
# set.seed(1212)
# #split the data into training and test data
# sample_size <- floor(0.75 * nrow(Boston))
# training_index <- sample(nrow(Boston), size = sample_size)
# train <- Boston[training_index, ]
# test <- Boston[-training_index, ]
# 
# 
# #set a seed so you can reproduce the results
# set.seed(1212)
# #split the data into training and test data
# sample_size <- floor(0.75 * nrow(miRNA))
# response1<-response[1:sample_size]
# sample_size2<-floor(0.25 * nrow(miRNA))
# train <- miRNA[response1, ]
# test <- miRNA[-response1, ]
# test<-test[1:sample_size2,]
# 
# model.net <- train(
#   hsa.mir.142 ~., data = miRNA, method = "glmnet",
#   trControl = trainControl("cv", number = 10),
#   tuneLength = 10)
# 
# model.net$bestTune
# 
# coef(model.net$finalModel, model.net$bestTune$lambda)
# x.test.net <- model.matrix(hsa.mir.142 ~., test)[,-1]
# predictions.net <- model.net %>% predict(x.test.net)
# data.frame(
#   RMSE.net = RMSE(predictions.net, test$hsa.mir.142),
#   Rsquare.net = R2(predictions.net, test$hsa.mir.142))
# 


miRNA_elnet <- train(
  hsa.mir.142 ~., data = miRNA, method = "glmnet",
  trControl = trainControl("cv", number = 10),
  tuneLength = 10)

ggplot(miRNA_elnet) +
  labs(title = "Elastic Net Regression Parameter Tuning", x = "lambda")



lll<-function(x){
  
  p<-miRNA %>% select(x)
  pp<-as.numeric(unlist(p))

  varmtx<-miRNA %>% select(!(p))
  varmtx<-as.matrix(varmtx)
  
  miRNA_elnet <- train(
    p ~., data = miRNA, method = "glmnet",
    trControl = trainControl("cv", number = 10),
    tuneLength = 10)
  
  alpha2<-miRNA_elnet$bestTune$alpha
  lambda2<-miRNA_elnet$bestTune$lambda
  elnet<-glmnet(scale(varmtx),p,alpha = alpha2,lambda = lambda2)
  
  el_coef<-coef(elnet)
  el_coef <- as.matrix(el_coef)
  el_coef<-as.data.frame(el_coef)
  el_coef <- subset(el_coef, rownames(el_coef) != "(Intercept)")
  
  el_coef <- el_coef %>% 
    filter(!is.na(el_coef)) %>% 
    filter_at(vars(starts_with("s")), all_vars(. != 0))
  
  el_coef<-el_coef%>%arrange(desc(s0))
  
  el_coef
}

lll(hsa.mir.142)


# 

p<-miRNA %>% select(hsa.mir.142)
p<-as.numeric(unlist(p))
response <- p
rm(p)
varmtx<-miRNA %>% select(!(hsa.mir.142))
varmtx<-as.matrix(varmtx)

miRNA_elnet <- train(
  hsa.mir.142 ~., data = miRNA, method = "glmnet",
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

el_coef



#############

elasticnet(hsa.mir.142~.,miRNA,lambda = "lambda.min")

plot.elasticnet(result)

summary.elasticnet(result)












