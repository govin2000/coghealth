

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

library(gbm)

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


#function to load the multiple imputation data
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
models <- c("class ~ drage_12 + drsex_12 + Living_arrange + ESB + educatio + income_cat + irsad",
            
            "class ~ SIDN3_N10 + PPDN_N10 + DWDN_N10 + RELU_PR_N10 + COLU_PR_N10 + PALU_PR_N10 + BSLU_PR_N10 + PALU_DT_A 
          + BSLU_DR_A + AP_NO2  + AP_PM25 + AP_TL_DT_A  + AP_RLMA_DN_C2 + AP_RLMI_DN_C2",
            
            
            "class ~ gq22_time_12 + gq28_time_12  + gq30_wtran_12 + hkaq17c_12 + hkaq17d_12 + sit_trans_day
          + sit_screen_day + sit_comput_day + sit_work_day + sit_other_day",
            
            
            "class ~ drage_12 + drsex_12 + Living_arrange + ESB + educatio + income_cat + irsad 
          + SIDN3_N10 + PPDN_N10 + DWDN_N10 + RELU_PR_N10 + COLU_PR_N10 + PALU_PR_N10 + BSLU_PR_N10 + PALU_DT_A 
          + BSLU_DR_A + AP_NO2  + AP_PM25 + AP_TL_DT_A  + AP_RLMA_DN_C2 + AP_RLMI_DN_C2", 
            
            "class ~ drage_12 + drsex_12 + Living_arrange + ESB + educatio + income_cat + irsad 
            +gq22_time_12 + gq28_time_12  + gq30_wtran_12 + hkaq17c_12 + hkaq17d_12 + sit_trans_day
          + sit_screen_day + sit_comput_day + sit_work_day + sit_other_day",
            
            "class ~ SIDN3_N10 + PPDN_N10 + DWDN_N10 + RELU_PR_N10 + COLU_PR_N10 + PALU_PR_N10 + BSLU_PR_N10 + PALU_DT_A 
          + BSLU_DR_A + AP_NO2  + AP_PM25 + AP_TL_DT_A  + AP_RLMA_DN_C2 + AP_RLMI_DN_C2+ gq22_time_12 + gq28_time_12  + gq30_wtran_12 + hkaq17c_12 + hkaq17d_12 + sit_trans_day
          + sit_screen_day + sit_comput_day + sit_work_day + sit_other_day",
            
            
            "class ~ drage_12 + drsex_12 + Living_arrange + ESB + educatio + income_cat + irsad 
          + SIDN3_N10 + PPDN_N10 + DWDN_N10 + RELU_PR_N10 + COLU_PR_N10 + PALU_PR_N10 + BSLU_PR_N10 + PALU_DT_A 
          + BSLU_DR_A + AP_NO2  + AP_PM25 + AP_TL_DT_A  + AP_RLMA_DN_C2 + AP_RLMI_DN_C2 + gq22_time_12
          + gq28_time_12  + gq30_wtran_12 + hkaq17c_12 + hkaq17d_12 + sit_trans_day  
          + sit_screen_day + sit_comput_day + sit_work_day + sit_other_day"
            
)




# Define grid for models defined above

nnGrid <- expand.grid(.decay = c(0.01, 0.05,0.1), .size = c(5, 6, 7,8,9,10))
gbmGrid <-  expand.grid(interaction.depth = c(1, 5, 9), n.trees = (1:30)*50, shrinkage = c(0.01, 0.05, 0.1), n.minobsinnode = c(15))
ctrl <- fitControl <- trainControl(method = "repeatedcv", number = 10, repeats = 10, classProbs = TRUE, 
                                   summaryFunction = twoClassSummary,  allowParallel = T)





cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)


auc <- array(NA, dim = c(10, 7, 4))
se <- array(NA, dim = c(10, 7, 4))

df_bic <- data.frame()
dfn <-data.frame()
df_allad<- data.frame()
bic_mod <- c ('EEV', 'EEV', 'EEV', 'EEV', 'EEV', 'EEV',  'EEV', 'EEV', 'EEV', 'EEV')

All_Data <-list()
for (vl in 1:10) {
  
  ad_data <-data.frame(data[vl])
  
  names(ad_data) <- sub(".*\\.", "", names(ad_data))
  cog_cluster <- ad_data %>%
    select(cvltscore_12, sdmtscore_12) %>%
    mutate_all(list(scale))
  
   set.seed(124)
   xv <- mclustBIC(cog_cluster)
   y <- xv %>%
     as.data.frame.matrix() %>%
     rownames_to_column("n_profiles") 
   
   df_y<-data.frame(y, imputation=vl)
   df_bic<-rbind(df_bic, df_y)
  
   set.seed(124)
   mod1 <- Mclust(cog_cluster, modelNames = bic_mod[vl], G = 2, x = xv)
   
   dfx <- data.frame(mod1$classification) 
   dft<-dfx %>% group_by(mod1.classification) %>% tally()
   
   dfn[vl,1]<-dft$n[1]
   dfn[vl,2]<-dft$n[2]
   
   
   ad_data$class <- as.factor( mod1$classification)
   
   ad_data$class<- factor(ad_data$class, 
                          levels = c(1, 2), 
                          labels = c("One", "Two"))
   All_Data[[vl]] <- ad_data
   
   df_addata<-ad_data %>% select(sdmtscore_12, cvltscore_12, class)
   df_allad<-rbind(df_allad, data.frame(df_addata, imputation = vl))
   
}

write.csv(df_bic,'ImputeData_BIC.csv')
write.csv(dfn,'classnumber_BIC.csv')
write.csv(df_allad,'df_allad.csv')
write.csv(dfn,'class_numbers.csv')






#Script for plotting the BIC values
bic<-read.csv('ImputeData_BIC.csv')
bic_sel<- bic %>%
  select("n_profiles", "EEV", "EEE", "VEE", "imputation") %>%
  group_by(n_profiles) %>%
  abs()%>% 
  mutate(
    maxEEV = max(EEV),
    minEEV = min(EEV),
    maxEEE = max(EEE),
    minEEE = min(EEE),
    maxVEE = max(VEE),
    minVEE = min(VEE)
    
  ) 

plt<- ggplot(bic_sel, aes(x=n_profiles, group=imputation)) +
  geom_ribbon(aes(ymin=minEEV, ymax=maxEEV), alpha=0.1) +
  ylab("BIC") +
  xlab ("Number of profiles")+
  scale_x_continuous(breaks=c(1:9))+
  theme_pubclean()+
  theme(text = element_text(size = 20))

ggsave('BIC_EEV.pdf',plt)


testd<-All_Data[[5]]


ggplot(testd, aes(y = sdmtscore_12, x = class, group = class, color= class)) + 
  geom_boxplot(fill="lightgray", width =0.2)+ 
  theme_pubclean()+
  theme(text = element_text(size = 30)) +
  ylab("SDMT")+
  xlab("Class")

df_allad %>% select(class, imputation) %>% group_by(class, imputation) %>%tally


plt <- ggplot(df_allad, aes(y = sdmtscore_12, x = class, fill=imputation)) + 
  geom_boxplot(fill="lightgray", width =0.2)+ 
  facet_wrap(~imputation)+
  theme_pubclean()+
  theme(text = element_text(size = 30)) +
  ylab("SDMT")+
  xlab("Class")
ggsave('SDMT_ClassSeparation.pdf',plt)

plt<-ggplot(df_allad, aes(y = cvltscore_12, x = class, fill=imputation)) + 
  geom_boxplot(fill="lightgray", width =0.2)+ 
  facet_wrap(~imputation)+
  theme_pubclean()+
  theme(text = element_text(size = 30)) +
  ylab("CVLT")+
  xlab("Class")
ggsave('CVLT_ClassSeparation.pdf',plt)




cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)
auc <- array(NA, dim = c(1, 7, 4))
se <- array(NA, dim = c(1, 7, 4))

df_total = data.frame()
df_vimp = data.frame()

vals<- c(1,2,3,4,5,6,7,8,9,10)  
#vals <-c(7)
  for (vl in vals) {
    
    ad_data <-c()
    
    ad_data <-All_Data[[vl]]
    
    set.seed(110)
    trainIndex <- createDataPartition(ad_data$class,p=.8, list=FALSE)
    trainData <- ad_data[trainIndex,]
    testData  <- ad_data[-trainIndex,]
    
    model_weights <- ifelse(trainData$class == "One",
                            (1/table(trainData$class)[1]) * 0.5,
                            (1/table(trainData$class)[2]) * 0.5)
    
    
      for (i in 1:7) {
        
        set.seed(110) #for reproducability
        i_new <- gsub("[\r\n]", "", models[i]) 
        
        
        
        # GBM
        set.seed(110) #for reproducability
        gbm.model <- train(formula(i_new),  data =trainData, method = 'gbm', weights=model_weights, tuneGrid=gbmGrid, metric = "ROC", trControl =  ctrl)    
        roc <-test_roc(gbm.model, testData)
        auc[vl,i,1] <- roc$auc
        se[vl,i,1] <- var(roc)
        
        df <- data_frame(tpr = roc$sensitivities, fpr = 1 - roc$specificities, model = 'GBM', Type =i, imputation = vl)
        df_total <- rbind(df_total,df)
        
        
        vimp<-summary(gbm.model, method = permutation.test.gbm, las = 2 )
        dfv<-data.frame(vimp, model='GBM', Type =i, imputation = vl)
        df_vimp <-rbind(df_vimp, dfv)
       
       
        # SVM
        set.seed(110) #for reproducability
        svm.model <- train(formula(i_new),  data =trainData, method = 'svmLinear', weights=model_weights, preProcess = c('center', 'scale'), metric = "ROC", trControl =  ctrl)    
        
        roc <-test_roc(svm.model, testData)
        auc[vl,i,2] <- roc$auc
        se[vl,i,2] <- var(roc)
        df <- data_frame(tpr = roc$sensitivities, fpr = 1 - roc$specificities, model = 'SVM', Type =i, imputation = vl)
        df_total <- rbind(df_total,df)
        
        
        # nnet
        set.seed(110) #for reproducability
        nnet.model <- train(formula(i_new),  data =trainData, method = 'nnet', tuneGrid = nnGrid, weights=model_weights, preProcess = c('center', 'scale'), metric = "ROC", trControl =  ctrl)    
        roc <-test_roc(nnet.model, testData)
        auc[vl,i,3] <- roc$auc
        se[vl,i,3] <- var(roc)
        
        df <- data_frame(tpr = roc$sensitivities, fpr = 1 - roc$specificities, model = 'nnet', Type =i, imputation = vl)
        df_total <- rbind(df_total,df)
        
        # Linear model
        set.seed(110) #for reproducability
        lm.model <- train(formula(i_new),  data =trainData, method = 'glm',family = 'quasibinomial',  weights=model_weights, metric = "ROC", trControl =  ctrl)  
        roc <-test_roc(lm.model, testData)
        auc[vl,i,4] <- roc$auc
        se[vl,i,4] <- var(roc)
        
        df <- data_frame(tpr = roc$sensitivities, fpr = 1 - roc$specificities, model = 'LM', Type =i, imputation = vl)
        df_total <- rbind(df_total,df)
        
      
      }
      
  }
  

saveRDS(auc, 'AUC3.rds')
saveRDS(se, 'SE3.rds')
saveRDS(df_total,'roc_df3.rds')
saveRDS(df_vimp,'vimp3.rds')



stopCluster(cl)




test_roc <- function(model, data) {
  
  roc(data$class,
      predict(model, data, type='prob')[, "Two"])
  
}



r2se <- function (measured, predicted){
  
  fit <-lm(predicted ~ measured)
  sfit <-summary(fit)
  r2<-data.frame(sfit$r.squared, sfit$coefficients[2,2])
  colnames(r2) = c("R2", "SE")
 
  return(r2)
}


