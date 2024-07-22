# Script name:  combat_sites
#  Created on:  26/02/2024
#      Author:  Clara Riegis 
#     Version:  1.0
#       Notes:  /

#_______________________________________________________________________________
#install.packages("lmerTest",dependencies=TRUE)
#install.packages('patchwork', dependencies = TRUE)

library(lme4)
library(lmerTest)
library(dplyr)
library(tidyr)
library(lattice)
library(patchwork)
library(ggplot2)
library(officer)
library(MuMIn)
library(openxlsx)

#_______________________________________________________________________________

setwd("~/Desktop/Cambridge/mphil_project/ABCD/abcd_format_data")  

#  Demographics. 
demographics <- read.csv('demographics.csv')

#  Covariates.
df_covar <- read.csv("covar.csv", header = TRUE)
# ROI names.
roinames <- readLines("roi_names.txt") # paste(, collapse = " ")

# alpha = 0.000058055152395

setwd("~/Desktop/Cambridge/mphil_project")  

thresholds = c(0.00031, 6e-05)
scales = c('scaled') #, 'no_scaled'
hamonizations = c('combat') #, 'no_harmo'
models = c('pgs_age') # 'interview_age',
dx_on_off = c('dx', 'no_dx')


      
 

# Iterate through the condition-measure-ROI 
# combination selected during cross-sectional 
# analysis. 
wb <- createWorkbook()

for (alpha in thresholds){
  
  for (norm in hamonizations) {
    
    for (scale in scales){
    
      for (dx in dx_on_off){
      
        # create files if they don't exist. 
        dir_path <- paste("Figures/mixed_model/", alpha, '/', norm, '/', scale, '/', dx, sep = '')
        if (!dir.exists(dir_path)) {
          dir.create(dir_path, recursive = TRUE)
          message("Directory created")
        } else {
          message("Directory already exists")
        }
        
        # Initialize storage variables.
        pred <- list()
        est_int <- list()
        est_age <- list()
        pval_int <- list()
        pval_age <- list()
        r_squared <- list()
        intercept <- list()
        data_list <- list()
        df_cross_res <- read.csv(paste('Figures/cross_sect/',norm, '/', scale, '/',
                                       alpha, '/result_table_noduplicates.csv', sep = ''))
        
          for (i in seq_len(nrow(df_cross_res))) { 
            
            xfile_name = df_cross_res$measure[[i]]
            print(xfile_name)
            
            roi = df_cross_res$roi[[i]]
            print(roi)
            
            condition = df_cross_res$condition[[i]]
            df_pgs = read.csv(paste("ABCD/PGS/all_final_scores/", 
                                    condition, "_final_scores.csv", sep = ''))
          
              
            # Load the data for the corresponding harmo, scaled, and measures. 
            path <- paste("ABCD/abcd_format_data/", norm, "/", scale, "/",  
                          xfile_name, ".csv", sep = '') 
            data_list[[xfile_name]] <- read.csv(path)  
            
            
            # Merge the demographics, struct. measures, pgs and covariates together.
            data_merged1 = merge(demographics, data_list[[xfile_name]], 
                                 by = "ID_event", 
                                suffixes = c("", "")) 
            # Make sure to delete the duplicates as follow:
            data_merged1 <- data_merged1[, !duplicated(names(data_merged1))]
            # Merging the remaining df (pgs and covariates). 
            data_merged2 <- merge(merge(data_merged1, df_pgs, by = "ID", 
                                        suffixes = c("", "")), df_covar, 
                                  by = "ID", 
                                  suffixes = c("", ""))
            # Getting rid of the duplicated once again.
            data_merged2 <- data_merged2[, !duplicated(names(data_merged2))]
          
            data_merged3 <- data_merged2
              
            # Rename the column of the measure we are looking into. 
            colnames(data_merged3)[which(colnames(data_merged2) == 
                                             roi)] <- 'measure' 
              
            # print(data_merged3[1, 'measure'])
            # print(list_measures$V1[i])
              
            data_merged3 <- data_merged3 %>% drop_na()
            # time affects the pheno?
            # use age rather than eventn num -> scale the age / can try with age^2
              
            data_merged3$eventnum <- factor(data_merged3$eventnum, 
                                            levels = unique(data_merged3$eventnum))
            
            
            
            # Define common terms
            common_terms <- c("interview_age", "SCORE", "gender", "acq", 
                              "fd_1", "fd_max_1", "leftplusright", "C1", 
                              "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", 
                              "C10", "C11", "C12", "C13", "C14", "C15", "C16", 
                              "C17", "C18", "C19", "C20",
                              "SCORE * interview_age", "(1 | ID)")
            
            # Define the formula based on the conditions
            if (norm == 'no_harmo') {
              if (dx == 'dx') {
                formula <- paste("measure ~", paste(common_terms, collapse = " + "), 
                                 "+ site_id_l + dx")
                
              } else if (dx == 'no_dx') {
                formula <- paste("measure ~", paste(common_terms, collapse = " + ")
                                 , "+ site_id_l")
              }
            } else {
              if (dx == 'dx') {
                formula <- paste("measure ~", paste(common_terms, collapse = " + "), 
                                 "+ dx")
                
              } else if (dx == 'no_dx') {
                formula <- paste("measure ~", paste(common_terms, collapse = " + "))
              }
            }
            
            
            mlr = lmer(formula, data = data_merged3)     # Dataframe.
          
                
            summary = summary(mlr) 
          
            # Store the predicted values to later create a calibration plot.
            pred[[i]] <- predict(mlr, re.form = ~ 0)
            
            # Store the variance explained by the model (conditional r-squared).
            r_squared[[i]] <- r.squaredGLMM(mlr, conditional = TRUE)[2]
              
            # Select the coefficients for later plotting 
            # (interaction and age only). 
            coeff_int <- as.data.frame(summary$coefficients) %>% select('Estimate') 
            est_int[[i]] <- coeff_int[rownames(coeff_int) %in% c("interview_age:SCORE"),]
            est_age[[i]] <- coeff_int[rownames(coeff_int) %in% c("interview_age"),]
            
            # Select the intercepts. 
            intercept[[i]] <- coeff_int[rownames(coeff_int) %in% c("(intercept)"),]
          
            # Now, select the corresponding p-values. 
            pval <- as.data.frame(summary$coefficients) %>% select('Pr(>|t|)')
            pval_int[[i]] <- pval[rownames(pval) %in% c("interview_age:SCORE"),]
            pval_age[[i]] <- pval[rownames(pval) %in% c("interview_age"),]
            
          }
        # Storing everything with the cross-sectional results. 
        df_cross_res$mlm_est_int = est_int
        df_cross_res$mlm_pval_int = pval_int
        
        df_cross_res$mlm_est_age = est_age
        df_cross_res$mlm_pval_age = pval_age
        
        df_cross_res$mlm_r2 = r_squared
        
        df_cross_res$mlm_intercept = intercept
        
  
        # Save this cross-section + longitudinal results df for later plotting. 
        file_path <- paste(dir_path, "/mlm_results.csv", sep = '') 
        df1 <- apply(df_cross_res,2,as.character)
        write.csv(df1, file = file_path, row.names = FALSE)
        
        sheet_name <- paste(norm, scale, dx, alpha, sep = '-') 
        addWorksheet(wb, sheet_name)
        writeData(wb, sheet = sheet_name, df_cross_res)
      }
      
    }
    
  }
  
}

saveWorkbook(wb, "Figures/mixed_model/results_tables.xlsx", overwrite = TRUE)



