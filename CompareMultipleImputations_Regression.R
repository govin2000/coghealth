

#Load necessary libraries in R
library(plyr)

library(caret)
library(doParallel)
library(ModelMetrics)
library(tidyLPA)
library(dplyr)
library(mclust)
library(tidyr)
library(tidyverse)
library(pROC)
library(ggpubr)
library(foreach)


# Clean theme for publication quality figures
theme_set(theme_pubclean())


#Set working directory data filenames 
setwd("/Users/gopoudel/Documents/Research/AusDiab/data")
filename <- 'Imputations2.RData'
varname <- 'imp.data'



#Function to estimate linear equation model
lm_eqn <- function(y,x){
  m <- lm(y ~ x);
  eq <- substitute(~~italic(r)^2~"="~r2, 
                   list(r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


#function to load the data
myLoad <- function(filename, variable) {
  tmp.env <- new.env()
  varnames <- load(filename, envir=tmp.env)
  if (missing(variable)) {
    ## Choose the variable name as the only variable in the file or
    ## give an error.
    if (length(varnames) == 1) {
      return(get(varnames, envir=tmp.env))
    } else {
      stop("More than one variable present in the file, please specify which variable to extract.  Choices are: ",
           paste(varnames, sep=", "))
    }
  } else {
    return(get(variable, envir=tmp.env))
  }
}



#Load the data

data <- myLoad(filename, varname)


#Define machine learing models formula
models_sdmtscore_12 <- c("sdmtscore_12 ~ drage_12 + drsex_12 + Living_arrange + ESB + educatio + income_cat + irsad",
                         
                         "sdmtscore_12 ~ SIDN3_N10 + PPDN_N10 + DWDN_N10 + RELU_PR_N10 + COLU_PR_N10 + PALU_PR_N10 + BSLU_PR_N10 + PALU_DT_A 
              + BSLU_DR_A + AP_NO2  + AP_PM25 + AP_TL_DT_A  + AP_RLMA_DN_C2 + AP_RLMI_DN_C2",
                         
                         
                         "sdmtscore_12 ~ gq22_time_12 + gq28_time_12  + gq30_wtran_12 + hkaq17c_12 + hkaq17d_12 + sit_trans_day
              + sit_screen_day + sit_comput_day + sit_work_day + sit_other_day",
                         
                         
                         "sdmtscore_12 ~ drage_12 + drsex_12 + Living_arrange + ESB + educatio + income_cat + irsad 
              + SIDN3_N10 + PPDN_N10 + DWDN_N10 + RELU_PR_N10 + COLU_PR_N10 + PALU_PR_N10 + BSLU_PR_N10 + PALU_DT_A 
              + BSLU_DR_A + AP_NO2  + AP_PM25 + AP_TL_DT_A  + AP_RLMA_DN_C2 + AP_RLMI_DN_C2", 
                         
                         "sdmtscore_12 ~ drage_12 + drsex_12 + Living_arrange + ESB + educatio + income_cat + irsad + gq28_time_12  
                + gq30_wtran_12 + hkaq17c_12 + hkaq17d_12 + sit_trans_day
              + sit_screen_day + sit_comput_day + sit_work_day + sit_other_day",
                         
                         "sdmtscore_12 ~ SIDN3_N10 + PPDN_N10 + DWDN_N10 + RELU_PR_N10 + COLU_PR_N10 + PALU_PR_N10 + BSLU_PR_N10 + PALU_DT_A 
              + BSLU_DR_A + AP_NO2  + AP_PM25 + AP_TL_DT_A  + AP_RLMA_DN_C2 + AP_RLMI_DN_C2+ gq22_time_12 + gq28_time_12  + gq30_wtran_12 + hkaq17c_12 + hkaq17d_12 +
                         sit_trans_day + sit_screen_day + sit_comput_day + sit_work_day + sit_other_day",
                         
                         
                         "sdmtscore_12 ~ drage_12 + drsex_12 + Living_arrange + ESB + educatio + income_cat + irsad 
              + SIDN3_N10 + PPDN_N10 + DWDN_N10 + RELU_PR_N10 + COLU_PR_N10 + PALU_PR_N10 + BSLU_PR_N10 + PALU_DT_A 
              + BSLU_DR_A + AP_NO2  + AP_PM25 + AP_TL_DT_A  + AP_RLMA_DN_C2 + AP_RLMI_DN_C2 + gq22_time_12
              + gq28_time_12  + gq30_wtran_12 + hkaq17c_12 + hkaq17d_12 + sit_trans_day  
              + sit_screen_day + sit_comput_day + sit_work_day + sit_other_day"
                         
)

models_cvltscore_12 <- c(
  
  "cvltscore_12 ~ drage_12 + drsex_12 + Living_arrange + ESB + educatio + income_cat + irsad",
  
  "cvltscore_12 ~ SIDN3_N10 + PPDN_N10 + DWDN_N10 + RELU_PR_N10 + COLU_PR_N10 + PALU_PR_N10 + BSLU_PR_N10 + PALU_DT_A 
              + BSLU_DR_A + AP_NO2  + AP_PM25 + AP_TL_DT_A  + AP_RLMA_DN_C2 + AP_RLMI_DN_C2",
  
  
  "cvltscore_12 ~ gq22_time_12 + gq28_time_12  + gq30_wtran_12 + hkaq17c_12 + hkaq17d_12 + sit_trans_day
              + sit_screen_day + sit_comput_day + sit_work_day + sit_other_day",
  
  
  "cvltscore_12 ~ drage_12 + drsex_12 + Living_arrange + ESB + educatio + income_cat + irsad 
              + SIDN3_N10 + PPDN_N10 + DWDN_N10 + RELU_PR_N10 + COLU_PR_N10 + PALU_PR_N10 + BSLU_PR_N10 + PALU_DT_A 
              + BSLU_DR_A + AP_NO2  + AP_PM25 + AP_TL_DT_A  + AP_RLMA_DN_C2 + AP_RLMI_DN_C2", 
  
  
  "cvltscore_12 ~ drage_12 + drsex_12 + Living_arrange + ESB + educatio + income_cat + irsad 
                + gq22_time_12 + gq28_time_12  + gq30_wtran_12 + hkaq17c_12 + hkaq17d_12 + sit_trans_day
              + sit_screen_day + sit_comput_day + sit_work_day + sit_other_day",
  
  "cvltscore_12 ~ SIDN3_N10 + PPDN_N10 + DWDN_N10 + RELU_PR_N10 + COLU_PR_N10 + PALU_PR_N10 + BSLU_PR_N10 + PALU_DT_A 
              + BSLU_DR_A + AP_NO2  + AP_PM25 + AP_TL_DT_A  + AP_RLMA_DN_C2 + AP_RLMI_DN_C2+ gq22_time_12 + gq28_time_12  + 
               gq30_wtran_12 + hkaq17c_12 + hkaq17d_12 + sit_trans_day
              + sit_screen_day + sit_comput_day + sit_work_day + sit_other_day",
  
  
  "cvltscore_12 ~ drage_12 + drsex_12 + Living_arrange + ESB + educatio + income_cat + irsad 
              + SIDN3_N10 + PPDN_N10 + DWDN_N10 + RELU_PR_N10 + COLU_PR_N10 + PALU_PR_N10 + BSLU_PR_N10 + PALU_DT_A 
              + BSLU_DR_A + AP_NO2  + AP_PM25 + AP_TL_DT_A  + AP_RLMA_DN_C2 + AP_RLMI_DN_C2 + gq22_time_12
              + gq28_time_12  + gq30_wtran_12 + hkaq17c_12 + hkaq17d_12 + sit_trans_day  
              + sit_screen_day + sit_comput_day + sit_work_day + sit_other_day"
  
)
# Define grid for models defined above

nnGrid <- expand.grid(.decay = c(0.01, 0.05,0.1), .size = c(5, 6, 7,8,9,10))
gbmGrid <-  expand.grid(interaction.depth = c(1, 5, 9), n.trees = (1:30)*50, shrinkage = c(0.01, 0.05, 0.1), n.minobsinnode = c(15))
ctrl <- fitControl <- trainControl(method = "repeatedcv", number = 10, repeats = 10, allowParallel = T)



cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)


r2 <- array(NA, dim = c(10, 7, 4))
se <- array(NA, dim = c(10, 7, 4))

for (vl in 1:10) {
  
  ad_data <-data.frame(data[vl])
  names(ad_data) <- sub(".*\\.", "", names(ad_data))
  
  set.seed(124)
  trainIndex <- createDataPartition(ad_data$sdmtscore_12,p=.8, list=FALSE)
  trainData <- ad_data[trainIndex,]
  testData  <- ad_data[-trainIndex,]
  
  imp <- c()
  for (i in 1:7) {
    
  
    
    set.seed(124) #for reproducability
    i_new <- gsub("[\r\n]", "", models_sdmtscore_12[i]) 
    
    
    # gbm
    set.seed(124) #for reproducability
    gbm.model <- train(formula(i_new),  data =trainData, method = 'gbm', tuneGrid=gbmGrid, metric = "RMSE", trControl =  ctrl)    
    # Predict using linear model
    gbm.pred <- predict(gbm.model, newdata = testData)
    est_r2se<-r2se(testData$sdmtscore_12,  gbm.pred)
    
    r2[vl,i,1] <- est_r2se$R2
    se[vl,i,1] <- est_r2se$SE
    
    
    # SVM
    set.seed(124) #for reproducability
    svm.model <- train(formula(i_new),  data =trainData, method = 'svmLinear',  preProcess = c('center', 'scale'), metric = "RMSE", trControl =  ctrl)    
    svm.pred <- predict(svm.model, newdata = testData)
    est_r2se<-r2se(testData$sdmtscore_12,  svm.pred)
    
    r2[vl,i,2] <- est_r2se$R2
    se[vl,i,2] <- est_r2se$SE
    
    
    
    
    # nnet
    set.seed(124) #for reproducability
    nnet.model <- train(formula(i_new),  data =trainData, method = 'nnet', tuneGrid = nnGrid, preProcess = c('center', 'scale'),
                        metric = "RMSE", trControl =  ctrl, linout=TRUE)    
    # Predict using linear model
    nnet.pred <- predict(nnet.model, newdata = testData)
    est_r2se<-r2se(testData$sdmtscore_12,  nnet.pred)
    
    r2[vl,i,3] <- est_r2se$R2
    se[vl,i,3] <- est_r2se$SE
    
    
    
    
    # Linear model
    lm.model <- train(formula(i_new),  data =trainData, method = 'glm', metric = "RMSE", trControl =  ctrl)    
    # Predict using linear model
    lm.pred <- predict(lm.model, newdata = testData)
    est_r2se<-r2se(testData$sdmtscore_12,  lm.pred)
    
    r2[vl,i,4] <- est_r2se$R2
    se[vl,i,4] <- est_r2se$SE

    
  }
  
}

saveRDS(r2, 'R2_SDMT.rds')
saveRDS(se, 'SE_SDMT.rds')

r2_sdmt<-readRDS('R2_SDMT.rds')
se_sdmt<-readRDS('SE_SDMT.rds')


stopCluster(cl)



cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)
r2 <- array(NA, dim = c(10, 7, 4))
se <- array(NA, dim = c(10, 7, 4))

for (vl in 1:10) {
  
  ad_data <-data.frame(data[vl])
  names(ad_data) <- sub(".*\\.", "", names(ad_data))
  
  set.seed(124)
  trainIndex <- createDataPartition(ad_data$cvltscore_12,p=.8, list=FALSE)
  trainData <- ad_data[trainIndex,]
  testData  <- ad_data[-trainIndex,]
  
  imp <- c()
  for (i in 1:7) {
    
    
    
    set.seed(124) #for reproducability
    i_new <- gsub("[\r\n]", "", models_cvltscore_12[i]) 
    
    
    # gbm
    set.seed(124) #for reproducability
    gbm.model <- train(formula(i_new),  data =trainData, method = 'gbm', tuneGrid=gbmGrid, metric = "RMSE", trControl =  ctrl)    
    # Predict using linear model
    gbm.pred <- predict(gbm.model, newdata = testData)
    est_r2se<-r2se(testData$cvltscore_12,  gbm.pred)
    
    r2[vl,i,1] <- est_r2se$R2
    se[vl,i,1] <- est_r2se$SE
    
    gbmImp=varImp(gbm.model)
    
    
    
    
    # SVM
    set.seed(124) #for reproducability
    svm.model <- train(formula(i_new),  data =trainData, method = 'svmLinear',  preProcess = c('center', 'scale'), metric = "RMSE", trControl =  ctrl)    
    svm.pred <- predict(svm.model, newdata = testData)
    est_r2se<-r2se(testData$cvltscore_12,  svm.pred)
    
    r2[vl,i,2] <- est_r2se$R2
    se[vl,i,2] <- est_r2se$SE
    
    
    
    
    # nnet
    set.seed(124) #for reproducability
    nnet.model <- train(formula(i_new),  data =trainData, method = 'nnet', tuneGrid = nnGrid, preProcess = c('center', 'scale'),
                        metric = "RMSE", trControl =  ctrl, linout=TRUE)    
    # Predict using linear model
    nnet.pred <- predict(nnet.model, newdata = testData)
    est_r2se<-r2se(testData$cvltscore_12,  nnet.pred)
    
    r2[vl,i,3] <- est_r2se$R2
    se[vl,i,3] <- est_r2se$SE
    
    
    
    
    # Linear model
    lm.model <- train(formula(i_new),  data =trainData, method = 'glm', metric = "RMSE", trControl =  ctrl)    
    # Predict using linear model
    lm.pred <- predict(lm.model, newdata = testData)
    est_r2se<-r2se(testData$cvltscore_12,  lm.pred)
    
    r2[vl,i,4] <- est_r2se$R2
    se[vl,i,4] <- est_r2se$SE
    
    
  }
  
}

saveRDS(r2, 'R2_cvlt.rds')
saveRDS(se, 'SE_cvlt.rds')

r2_cvlt<-readRDS('R2_cvlt.rds')
se_cvlt<-readRDS('SE_cvlt.rds')

stopCluster(cl)



r2se <- function (measured, predicted){
  
  fit <-lm(predicted ~ measured)
  sfit <-summary(fit)
  r2<-data.frame(sfit$r.squared, sfit$coefficients[2,2])
  colnames(r2) = c("R2", "SE")
 
  return(r2)
}


