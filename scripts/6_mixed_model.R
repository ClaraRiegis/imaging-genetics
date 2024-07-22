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
#_______________________________________________________________________________

setwd("~/Desktop/Cambridge/mphil_project/ABCD/abcd_format_data")  

#  Demographics. 
demographics <- read.csv('demographics.csv')
#  Polygenic scores. 
df_pgs <- read.csv("autism_unstratified_finalscore.csv", header = TRUE)
#  Covariates.
df_covar <- read.csv("covar.csv", header = TRUE)
#  Phenotypes. 
str_meas <- c("area", "sulc", "thk", "vol")
# ROI names.
roinames <- readLines("roi_names.txt") # paste(, collapse = " ")

alpha = 0.01

setwd("~/Desktop/Cambridge/mphil_project")  

data_list <- list()
norms = c('combat')
scales = c('scaled')
models = c('interview_age','pgs_age')


  
  for (norm in norms){
    print(norm)
    
    for (scale in scales){
      print(scale)
      
      for (model in models) {
      
        all_est = list()
        all_pval = list()
        
          for (xfile_name in str_meas) { # Iterate over each file name
          
            # Load the data for the corresponding harmo, scaled, and measures. 
            path <- paste("ABCD/abcd_format_data/", norm, "/", scale, "/",  
                          xfile_name, ".csv", sep = '') 
            data_list[[xfile_name]] <- read.csv(path)  
            
            
            # Merge the demographics, struct. measures, pgs and covariates together.
            data_merged1 = merge(demographics, data_list[[xfile_name]], by = "ID_event", 
                                suffixes = c("", "")) 
            # Make sure to delete the duplicates as follow:
            data_merged1 <- data_merged1[, !duplicated(names(data_merged1))]
            # Merging the remaining df (pgs and covariates). 
            data_merged2 <- merge(merge(data_merged1, df_pgs, by = "ID", 
                                        suffixes = c("", "")), df_covar, by = "ID", 
                                  suffixes = c("", ""))
            # Getting rid of the duplicated once again.
            data_merged2 <- data_merged2[, !duplicated(names(data_merged2))]
            
            
            coeff_pval = list()
            coeff_est = list()
            mlr <- list()

            for (i in seq_along(roinames)) {  # Looping through the various ROIs. 
              print(roinames[i])
              
              
              data_merged3 <- data_merged2
              
              # Rename the column of the measure we are looking into. 
              colnames(data_merged3)[which(colnames(data_merged2) == 
                                             roinames[i])] <- 'measure' 
              
              # print(data_merged3[1, 'measure'])
              # print(list_measures$V1[i])
              
              data_merged3 <- data_merged3 %>% drop_na()
              # time affects the pheno?
              # use age rather than eventn num -> scale the age / can try with age^2
              
              data_merged3$eventnum <- factor(data_merged3$eventnum, levels = unique(data_merged3$eventnum))
              
              if (model == 'interview_age'){
                print(model)
                mlr[[xfile_name]] =   lmer(measure ~                # Phenotype.
                                           interview_age+
                                           dx +                     # Diagnosis
                                           sex +                    # Sex / gender?
                                           rsfmri_meanmotion +      # Framewise disp.
                                           interview_age * sex +
                                           (1 | ID),                # Random intercepts.
                                           data = data_merged3)     # Dataframe.
              
                
                fit_mlm <- data_merged3
                fit_mlm$pred <- predict(mlr[[xfile_name]], re.form = ~ 0)
                ggplot(fit_mlm, aes(x = interview_age, y = pred), alpha = 10) +
                  geom_point(aes(color = sex), size = 0.1) + 
                  geom_smooth(method = "lm", se = FALSE, aes(color = sex)) + 
                  #geom_smooth(method = "lm", se = FALSE, aes(color = 'all'))+
                  scale_color_manual(values = c( 'indianred2', 'steelblue3'))+
                  labs(x = "Age", y = paste("Predicted", xfile_name), color = "Sex") +
                  theme_minimal()
                
                summary = summary(mlr[[xfile_name]]) 
                  
                # Select the coefficients for later plotting. 
                coeff <- as.data.frame(summary$coefficients) %>% select('Estimate') 
                coeff_est[[i]] <- coeff[rownames(coeff) %in% c("interview_age:sexMale"),]
                
                pval <- as.data.frame(summary$coefficients) %>% select('Pr(>|t|)')
                coeff_pval[[i]] <- pval[rownames(pval) %in% c("interview_age:sexMale"),]
                
                }
            
            
              else if (model == 'pgs_age') {
                print(model)
                mlr[[xfile_name]] = lmer(measure ~                # Phenotype.
                                           interview_age +
                                           SCORE +                  # Polygenic Score. 
                                           dx +                     # Diagnosis
                                           sex +                    # Sex / gender?
                                           rsfmri_meanmotion +      # Framewise disp.
                                           C1 + C2 + C3 + C4 +      
                                           C5 + C6 + C7 + C8 +      # Principal
                                           C9 + C10 + C11 +         # Components.
                                           C12 + C13 + C14 + 
                                           C15 + C16 + C17 +
                                           C18 + C19 + C20 +
                                           SCORE * interview_age +
                                           (1 | ID),                # Random intercepts.
                                           data = data_merged3)     # Dataframe.
                
                
                summary = summary(mlr[[xfile_name]]) 
                
                
                # Select the coefficients for later plotting. 
                coeff <- as.data.frame(summary$coefficients) %>% select('Estimate') 
                coeff_est[[i]] <- coeff[rownames(coeff) %in% c("interview_age:SCORE"),]
                
                pval <- as.data.frame(summary$coefficients) %>% select('Pr(>|t|)')
                coeff_pval[[i]] <- pval[rownames(pval) %in% c("interview_age:SCORE"),]
                
                if (coeff_pval[[i]] <= alpha){
                  
                  fit_mlm <- data_merged3
                  fit_mlm$pred <- predict(mlr[[xfile_name]], re.form = ~ 0)
                  ggplot(fit_mlm, aes(x = interview_age, y = pred), alpha = 10) +
                    geom_point(aes(color = SCORE), size = 0.5) + 
                    scale_colour_continuous( low = "blue", high = "red", 
                                             space = "Lab", name = "PGS") +
                    #geom_line(aes(group = factor(ID), col = sex), linewidth=0.1) +
                    geom_smooth(method = "lm", se = FALSE) + # aes(color = sex)
                    #geom_smooth(method = "lm", se = FALSE, aes(color = 'all'))+
                    #scale_color_manual(values = c('purple4', 'indianred2', 'steelblue3'))+
                    labs(x = "Age", y = paste("predicted", xfile_name), color = "Sex") +
                    ggtitle(roinames[[i]])
                    theme_minimal()
                  
                  ggsave(paste("Figures/outputs/mixed_model/", norm, "/", 
                               scale, "/", model, "_",  roinames[i], "_line.png", sep =""))
                  
                  sink(paste("Figures/outputs/mixed_model/", norm, "/", 
                                                   scale, "/", model, "_",  roinames[i], 
                                                   "_summary.txt", sep =""))
                  print(summary)
                  sink()  # returns output to the console
                  
                  
                }
                
              }
              
            
            
             
              #significant_results <- summary$coefficients[summary$coefficients[, "Pr(>|t|)"] < alpha, ]
              #if (NCOL(significant_results) >= 2) {
              #  p_values <- list()
              #  t_values <- list()
              #  print(paste(" ___________ MODEL:", model, "| MEASURE:",
              #              xfile_name, "| ROI:",roinames[i] ," ___________"))
              #  # Iterate through significant predictors and provide interpretations
              #  for (i in seq_len(nrow(significant_results))) {
              #     predictor <- rownames(significant_results)[i]
              #    p_value <- significant_results[i, "Pr(>|t|)"]
              #    t_value <- significant_results[i, "t value"]
              #    t_values[[i]] <- t_value
              #    p_values[[i]] <- p_value
              #    interpretation <- paste("The predictor", predictor, "is statistically significant (t =", round(t_value, 5), ", p =", round(p_value, 5), ")")
              #    # Print interpretation
              #    print(interpretation)
              #  }
              #}
            
              #print("_____________________________________________________________")
            
            
            
              
            }
          
            
            all_est[[xfile_name]] = coeff_est
            all_pval[[xfile_name]] = coeff_pval
          
           
          
          }
        
        all_est[["rois"]] = as.list(roinames)
        all_pval[["rois"]] = as.list(roinames)
        
        df_est <- do.call(rbind, lapply(all_est, function(x) as.data.frame(t(x))))
        df_pval <- do.call(rbind, lapply(all_pval, function(x) as.data.frame(t(x))))
        
        df_est <- as.data.frame(df_est)
        df_pval <- as.data.frame(df_pval)
        
        df_est <- cbind(measures = rownames(df_est), df_est)
        rownames(df_est) <- 1:nrow(df_est)
        
        df_pval <- cbind(measures = rownames(df_pval), df_pval)
        rownames(df_pval) <- 1:nrow(df_pval)
        
        # Set the path where the results will be stored. 
        file_path <- paste("Figures/outputs/mixed_model/", norm, "/", scale, 
                           "/", model, "_estimates.csv", sep = '')
        df1 <- apply(df_est,2,as.character)
        write.csv(df1, file = file_path, row.names = FALSE)
        
        file_path <- paste("Figures/outputs/mixed_model/", norm, "/", scale, 
                           "/", model, "_p_values.csv", sep = '')
        df2 <- apply(df_pval,2,as.character)
        write.csv(df2, file = file_path, row.names = TRUE)
        
      }
    }
  }



if (is.data.frame(df)) {
  print("Your variable is a dataframe.")
} else {
  print("Your variable is not a dataframe.")
}
