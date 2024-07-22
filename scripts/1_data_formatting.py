#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 12:38:56 2024

Author: clarariegis

Purpose: Load, clean and format the phenotypic data 
        (measures of brain structure), and the demographic, PGS and covariates
        files that were made available. 
        The phenotypes were fetched directly from the ABCD release
        file available in the HPC storage. 
        The covariates and the ASD PGS were shared by Varun W.
        The demographics file was share by Richard B.
        Save the formatted data. At this stage all data is still raw 
        (it has not been scaled or anything else).

Description: 

    


"""



#%% - 0 - Libraries 

import pandas as pd
import numpy as np
import os
import glob
from sklearn.preprocessing import StandardScaler


os.chdir('/Users/clarariegis/Desktop/Cambridge/mphil_project/code/scripts')
from custom_functions import create_directory_if_not_exists

os.chdir('/Users/clarariegis/Desktop/Cambridge/mphil_project')

#%% - 1 - Load structural measures



# Polygenic Scores:
# Load the polygenic risk scores to keep participant retained after 
# the genetic quality control. 
path_autism = 'ABCD/PGS/autism_unstratified_finalscore.csv'
pgs_autism = pd.read_csv(path_autism, delim_whitespace=True)
# Change the particpant ID column name. 
pgs_autism = pgs_autism.rename(columns={'IID': 'ID'})
output = "ABCD/abcd_format_data/autism_unstratified_finalscore.csv"
pgs_autism.to_csv(output, index=False)


# Covariates:
# The covariates were separated into two different variables:
df_covar1 = pd.read_csv("ABCD/covar_demos/GWAS_ABCD_corticalcovar.txt", 
                        sep=' ')
df_covar2 = pd.read_csv("ABCD/covar_demos/GWAS_ABCD_corticalqcovar.txt", 
                        sep=' ')
# Replace the Xs in column names by Cs (components?).
df_covar2.rename(columns=lambda x: x.replace('X', 'C'), inplace=True)
# Putting them together:
df_covar = df_covar1.merge(df_covar2, how='inner')
sum(np.where(df_covar.isna() == True))  # Check for NaNs.
# Chance the participant ID column name. 
df_covar = df_covar.rename(columns={'IID': 'ID'})
# Save covariate file:
output = "ABCD/abcd_format_data/covar.csv"
df_covar.to_csv(output, index=False)



# Structural measures:
path_roi = 'ABCD/abcd_raw_data/mri_y_smr_*.csv'
all_measures = {}
str_meas = []

# Get the names of the ROI such that they can be used to be plotted.
data_dk = {'bankssts_left': 2.48,
           'caudalanteriorcingulate_left': 2.228,
           'caudalmiddlefrontal_left': 2.618,
           'cuneus_left': 2.137,
           'entorhinal_left': 3.332,
           'fusiform_left': 2.605,
           'inferiorparietal_left': 2.436,
           'inferiortemporal_left': 2.782,
           'isthmuscingulate_left': 2.138,
           'lateraloccipital_left': 2.258,
           'lateralorbitofrontal_left': 2.582,
           'lingual_left': 2.333,
           'medialorbitofrontal_left': 2.321,
           'middletemporal_left': 2.949,
           'parahippocampal_left': 2.731,
           'paracentral_left': 2.433,
           'parsopercularis_left': 2.561,
           'parsorbitalis_left': 2.715,
           'parstriangularis_left': 2.463,
           'pericalcarine_left': 1.978,
           'postcentral_left': 2.213,
           'posteriorcingulate_left': 2.206,
           'precentral_left': 2.73,
           'precuneus_left': 2.238,
           'rostralanteriorcingulate_left': 2.632,
           'rostralmiddlefrontal_left': 2.406,
           'superiorfrontal_left': 2.638,
           'superiorparietal_left': 2.244,
           'superiortemporal_left': 2.656,
           'supramarginal_left': 2.507,
           'frontalpole_left': 2.579,
           'temporalpole_left': 3.652,
           'transversetemporal_left': 2.567,
           'insula_left': 2.869,
           'bankssts_right': 2.579,
           'caudalanteriorcingulate_right': 2.501,
           'caudalmiddlefrontal_right': 2.649,
           'cuneus_right': 2.265,
           'entorhinal_right': 2.448,
           'fusiform_right': 2.602,
           'inferiorparietal_right': 2.424,
           'inferiortemporal_right': 2.609,
           'isthmuscingulate_right': 2.127,
           'lateraloccipital_right': 2.381,
           'lateralorbitofrontal_right': 2.533,
           'lingual_right': 2.424,
           'medialorbitofrontal_right': 2.266,
           'middletemporal_right': 2.741,
           'parahippocampal_right': 2.741,
           'paracentral_right': 2.45,
           'parsopercularis_right': 2.521,
           'parsorbitalis_right': 2.505,
           'parstriangularis_right': 2.442,
           'pericalcarine_right': 2.02,
           'postcentral_right': 2.171,
           'posteriorcingulate_right': 2.4,
           'precentral_right': 2.654,
           'precuneus_right': 2.348,
           'rostralanteriorcingulate_right': 2.541,
           'rostralmiddlefrontal_right': 2.362,
           'superiorfrontal_right': 2.642,
           'superiorparietal_right': 2.265,
           'superiortemporal_right': 2.587,
           'supramarginal_right': 2.459,
           'frontalpole_right': 2.551,
           'temporalpole_right': 3.486,
           'transversetemporal_right': 2.714,
           'insula_right': 2.994}
new_rois = list(data_dk.keys())



for filename in glob.glob(path_roi):
    print(filename)
    
    # Extract the name of the loaded measure - instead of having the  
    # whole file name. 
    str_meas = filename.split('_')[5] 
    
    # Load the file. 
    df_measures = pd.read_csv(filename, sep = ',')
    
    # Rename the ID column. 
    df_measures = df_measures.rename(columns={'src_subject_id': 'ID'})
    
    # Rename the ROI columns: 
    # - Extract the current ROI names: 
    old_rois = df_measures.iloc[:, 2:].columns
    # - Map the old and new ROIs into a dictionary: 
    dict_rois = dict(zip(old_rois, new_rois))
    # - Replace the original ROI names by more readable ones:
    df_measures = df_measures.rename(columns=dict_rois)
    
    # Create a column combining the participant ID and time point. 
    df_measures["ID_event"] = (df_measures["ID"] + "_" + 
                               df_measures["eventname"])
    
    df_measures = df_measures.sort_values(by='ID_event')
    
    # Only keep the participants that are in the PGS file. 
    df_measures = df_measures[df_measures['ID'].isin(pgs_autism['ID'])]
    
    # Only keep the participants that are in the covariate file. 
    df_measures = df_measures[df_measures['ID'].isin(df_covar['ID'])]
    
    # Numerise the event names. 
    df_measures['eventnum'] = (df_measures['eventname']
                             .replace(['baseline_year_1_arm_1', 
                                       '2_year_follow_up_y_arm_1', 
                                       '4_year_follow_up_y_arm_1'],
                                      [0, 1, 2], inplace=False))

    # Store all the dataframes toget into a dictionary. 
    all_measures[str_meas] = df_measures 
    


# Load the demographic file: 
demographics = pd.read_csv('ABCD/covar_demos/demographics_form.csv')
# Change ID column name. 
demographics = demographics.rename(columns={'src_subject_id': 'ID'})
# Create a column combining participant IDs and time point. 
demographics["ID_event"] = (demographics["ID"] + "_" + 
                            demographics["eventname"])
# Add demographic file to the dict with the measure dfs + order the rows. 
all_measures["demographics"] = demographics.sort_values(by='ID_event')  



# Keep only the participants that are present in all three dataframes and for
# and for the same time points: 
pcp_id = 'ID_event' # Target column. 
pcp_roi= 'thk'     # Reference df - doesn't matter which one.

# Find common participants
common_participants = set(all_measures[pcp_roi][pcp_id])
for key in all_measures:
    common_participants &= set(all_measures[key][pcp_id])

# Filter DataFrames
all_measures = {key: df[df[pcp_id].isin(common_participants)] for key, 
                      df in all_measures.items()}


# Save the files with the measures and demographics:
for key in all_measures:
    if key == 'demographics': 
        # create path to save output: 
        dir_path = f"ABCD/abcd_format_data/"
        create_directory_if_not_exists(dir_path)
        
        output = "ABCD/abcd_format_data/" + key + ".csv"
    else: 
        dir_path = "ABCD/abcd_format_data/no_harmo/no_scaled/"
        create_directory_if_not_exists(dir_path)
        
        output = "ABCD/abcd_format_data/no_harmo/no_scaled/" + key + ".csv"
    all_measures[key].to_csv(output, index=False)
    
all_measures["demographics"]["interview_age"].min()
all_measures["demographics"]["interview_age"].max()


# Scale the polygenic scores for every condition: 
str_pgs_list = []
df_all_pgs = {}
scaler = StandardScaler()

# Loading the PGS for each GWAS condition. 
path_pgs = "ABCD/PGS/all_final_scores/*.profile"
for filename in glob.glob(path_pgs): 

    # Extract the name of the condition.
    condition = filename.split('/')[3].split('_')[0]
    str_pgs_list += [condition]

    # Load the file. 
    df_pgs = pd.read_csv(filename, delim_whitespace=True)
    df_pgs = df_pgs.rename(columns={'IID': 'ID'})
    plt.hist(df_pgs["SCORE"], bins=60)
    # Store all the dfs in a dictionary.
    df_all_pgs[condition] = df_pgs

    df_pgs['SCORE'] = scaler.fit_transform(df_pgs['SCORE'].values.reshape(-1, 1))
    plt.hist(df_pgs["SCORE"], bins=60)
    # Change the particpant ID column name and save to csv. 
    output = f"ABCD/PGS/all_final_scores/{condition}_final_scores.csv"
    df_pgs.to_csv(output, index=False)














