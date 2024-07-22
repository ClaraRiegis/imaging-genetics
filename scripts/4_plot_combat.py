#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 15:21:21 2024

Author: clarariegis

Purpose: Visualise the association between age and brain structures. 

Description: - Plot of the lines of best fit for each brain region 
               (averaged across participants).
             - Show the slopes that are significant on the desikan killiany 
               atlas. 
             - Can adapt the code to show the unit of the measures if needed
               (not needed when using measures that have been scaled?????????)
             - Can plot Combat, no combat, scaled or not scaled data. 
             
"""

#%% - 0 - Libraries 

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import seaborn as sns
from pptx.util import Inches



os.chdir('/Users/clarariegis/Desktop/Cambridge/mphil_project/code/scripts')
from raincloud2 import raincloud2
from custom_functions import add_figures_to_presentation, create_directory_if_not_exists


os.chdir('/Users/clarariegis/Desktop/Cambridge/mphil_project')






#%%


# Get the regiong names as a list. 
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
roi_list = list(data_dk.keys())
xmin = {}
xmax = {}
fig_paths = []


normals = ['combat', 'no_harmo']
scales = ['scaled', 'no_scaled']

for harmo in normals:
    
    for scale in scales:
        # create path to save output: 
        dir_path = f'Figures/combat_plot/{harmo}/{scale}'
        create_directory_if_not_exists(dir_path)
        
        scatt_col = {}
        scatt_mark = {}
        print(harmo)
        
        path = f"ABCD/abcd_format_data/{harmo}/{scale}/*.csv"
        
        for filename in glob.glob(path):
            print(filename)
            
            # Extract the name of the loaded measure.
            str_meas = filename.split('/')[4].split('.')[0]
            
            # Load the file. 
            df_measures = pd.read_csv(filename, sep = ',')
            
    
            # Store all the dfs in a dictionary.
            # df_all_meas[str_meas] = df_measures
            
            
            
            sites = np.unique(df_measures['site'])
            
            random_colors = sns.color_palette(palette="light:#5A9", 
                                              n_colors = len(sites))
            
            
            # Average the regions of each participant within each site. 
            roi_avr = {}
            for site in sites:
                
                df_site = df_measures.loc[df_measures["site"] == site]
                roi_avr[site] = df_site[roi_list].mean(axis = 1).dropna().tolist()
                
                
                
                # Now prepare the colors and markers for plotting. 
                scatt_col[site] = ["black"] * len(roi_avr[site])
                scatt_mark[site] = ['o'] * len(roi_avr[site])
                
            xmin[str_meas] = min(min(values) for values in roi_avr.values())
            xmax[str_meas] = max(max(values) for values in roi_avr.values())
            
            fig, ax = plt.subplots(figsize=(10, 20))
            
            
            
            rain_colors = {'box_col': random_colors, 
                           'violin_col': random_colors, 
                           'scatt_col': scatt_col, 
                           'scatt_marker': scatt_mark,
                           'mark_size': 8}
    
    
    
            raincloud2(data_x = roi_avr, 
                       x_label = f'{str_meas}', 
                       title = ('Distribution of structural measures' +
                                f' ({str_meas}, {harmo}, {scale})'),
                       colors = True, 
                       ax = ax, 
                       size = 12, 
                       style = rain_colors)
        
        
            
    
            
            fig_path = f'{dir_path}/{harmo}_{scale}_{str_meas}.png'
            plt.savefig(fig_path)  # Save as PNG format
            fig_paths += [fig_path]
            
            
            
            plt.close()  # Close the figure

#%%


input_ppt_path = 'Figures/Slides/mphil_figures_clean2.pptx'
scale_ppt_path = 'Figures/Slides/scale.pptx'
output_ppt_path = 'Figures/Slides/mphil_figures_clean2.pptx'



# every path contains a different type of figure, for the same data. 
# we want them to be plotted in the same slide. 
list_fig_paths = [fig_paths]

# Desired width and height for the images (in inches)
desired_width = [Inches(10)]
desired_height = [Inches(10)]
left = [Inches(0)] # Adjust position as needed
top = [Inches(0)]



add_figures_to_presentation(input_ppt_path = input_ppt_path, 
                            scale_ppt_path = scale_ppt_path, 
                            fig_paths = list_fig_paths, 
                            desired_widths = desired_width, 
                            desired_heights = desired_height, 
                            left_positions = left, 
                            top_positions = top, 
                            output_ppt_path = output_ppt_path)











