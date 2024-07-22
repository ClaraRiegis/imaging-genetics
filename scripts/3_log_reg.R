# Script name:  x_log_reg
#  Created on:  14/05/2024
#      Author:  Clara Riegis 
#     Version:  1.0
#       Notes:  /


#_______________________________________________________________________________
#_______________________________________________________________________________



library(openxlsx)
library(dplyr)
library(lme4)
library(lmerTest)
library(tidyr)
library(lattice)
library(patchwork)
library(ggplot2)
library(officer)
library(MuMIn)
library(pROC)
if (!requireNamespace("scales", quietly = TRUE)) {
  install.packages("scales")
}

# Load necessary libraries
library(ggplot2)
library(scales)  # For scientific notation labels

# Set working directory
setwd("~/Desktop/Cambridge/mphil_project/ABCD/abcd_format_data")  


#_______________________________________________________________________________
#_______________________________________________________________________________


# Load data
demographics <- read.csv('demographics.csv')
df_covar <- read.csv("covar.csv", header = TRUE)
roinames <- readLines("roi_names.txt")

str_pgs <- c("ADHD", "ASD", "BIP", "MDD", "SCZ", "SUD")
str_meas <- c("area", "thk", "sulc", "vol")
alpha = 0.01

setwd("~/Desktop/Cambridge/mphil_project")  

pgs_list <- list()
data_list <- list()
log_reg <- list()
all_est <- data.frame(column1 = rep(NA, (length(str_meas) * length(str_pgs))))
all_pval <- data.frame(column1 = rep(NA, (length(str_meas) * length(str_pgs))))

# Load the measures. 
path <- paste("ABCD/abcd_format_data/combat/scaled/area.csv", sep = '') 
data_list[[xfile_name]] <- read.csv(path)  

#_______________________________________________________________________________
#_______________________________________________________________________________


# Function to preprocess data
preprocess_data <- function(demographics, data_list, df_covar, pgs_list, pgs, xfile_name, num_unique_timepoints) {
  # Merge demographics, structural measures, PGS, and covariates
  data_merged1 <- merge(demographics, data_list[[xfile_name]], by = "ID_event", suffixes = c("", ""), all = TRUE)
  data_merged1 <- data_merged1[, !duplicated(names(data_merged1))]
  data_merged2 <- merge(merge(data_merged1, pgs_list[[pgs]], by = "ID", suffixes = c("", ""), all = TRUE), df_covar, by = "ID", suffixes = c("", ""), all = TRUE)
  data_merged2 <- data_merged2[, !duplicated(names(data_merged2))]
  data_merged2 <- data_merged2 %>%
    arrange(desc(ID)) %>%
    na.omit() %>%
    group_by(ID) %>%
    mutate(num_time_points = case_when(
      length(unique(eventnum)) == num_unique_timepoints ~ 1,
      length(unique(eventnum)) == 1 ~ 0,
      TRUE ~ NA_real_
    )) %>%
    ungroup() %>%
    filter(num_time_points %in% c(0, 1))
  
  # Add a column with age at baseline for all participants
  baseline_ages <- data_merged2 %>%
    filter(eventname == 'baseline_year_1_arm_1') %>%
    select(ID, interview_age) %>%
    rename(baseline_age = interview_age)
  
  # Merge the baseline age column back to the original dataframe
  data_merged2 <- data_merged2 %>%
    left_join(baseline_ages, by = "ID") %>%
    distinct(ID, .keep_all = TRUE)
  
  return(data_merged2)
}

#_______________________________________________________________________________
#_______________________________________________________________________________


coeff_pval = list()
coeff_est = list()
mlr <- list()
i = 1
num_tps = c(2, 3)

wb <- createWorkbook()

for (pgs in str_pgs) {
  
  # Load the polygenic scores. 
  path <- paste("ABCD/PGS/all_final_scores/", pgs, "_final_scores.csv", sep = '') 
  pgs_list[[pgs]] <- read.csv(path)  
  
  for (x_tp in num_tps) {
    print("___________")
    print(x_tp)
    print("_______________________")
    
    # Preprocess the data
    data_merged2 <- preprocess_data(demographics, 
                                    data_list, 
                                    df_covar, 
                                    pgs_list, 
                                    pgs, 
                                    xfile_name, 
                                    x_tp)
    
    log_dat <- data_merged2 %>% filter(num_time_points %in% c(0, 1))
    
    # Fit the logistic regression model
    log_reg <- glm(num_time_points ~ 
                     SCORE + baseline_age + sex + 
                     C1 + C2 + C3 + C4 + C5 + C6 + 
                     C7 + C8 + C9 + C10 + C11 + 
                     C12 + C13 + C14 + C15 + C16 + 
                     C17 + C18 + C19 + C20,
                   data = log_dat, 
                   family = binomial)
    
    # Summary stats
    sum_stats <- coef(summary(log_reg))
    
    # Create a new sheet for each iteration
    sheet_name <- paste(pgs, "-", x_tp)
    addWorksheet(wb, sheet_name)
    
    # Save relevant results to the sheet
    writeData(wb, sheet = sheet_name, x = sum_stats, rowNames = TRUE)  # Example: saving coefficients
    
    
  }
  
}

# Save the Excel file
saveWorkbook(wb, file = "logistic_regression_results.xlsx", overwrite = TRUE)



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Apply symmetric log transformation
coeff_df$Estimate_transformed <- sym_log(coeff_df$Estimate)
coeff_df$Error_min <- sym_log(coeff_df$Estimate - 1.96 * coeff_df$`Std. Error`)
coeff_df$Error_max <- sym_log(coeff_df$Estimate + 1.96 * coeff_df$`Std. Error`)

# Plot the coefficients with symmetric log-transformed y-axis
ggplot(coeff_df, aes(x = Variable, y = Estimate_transformed)) +
  geom_point() +
  geom_errorbar(aes(ymin = Error_min, ymax = Error_max), width = 0.2) +
  scale_y_continuous(labels = function(x) scales::scientific(10^abs(x) - 1)) +  # Convert back to scientific notation
  geom_hline(yintercept = sym_log(0), linetype = "dashed", color = "darkgrey", size = 0.5) +  # Thinner, dark grey horizontal line at transformed zero
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Logistic Regression Coefficients", y = "Coefficient Estimate (Symmetric Log Scale)", x = "Predictor Variables")







all_est$estimate <- coeff_est
all_est$pval <- coeff_pval

# Set the path where the results will be stored. 
file_path <- paste("Figures/outputs/logist_reg/pval_n_estimates.csv", sep = '')
df1 <- apply(all_est, 2, as.character)
write.csv(df1, file = file_path, row.names = FALSE)







