#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 21:59:47 2024

@author: clarariegis
"""

#%% - 0 - Libraries 

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from matplotlib.pyplot import cm
from pptx import Presentation
from pptx.util import Inches
import seaborn as sns



os.chdir('/Users/clarariegis/Desktop/Cambridge/mphil_project/code/scripts')
from raincloud2 import raincloud2
# import func_plot_desikan2 
from func_plot_desikan2 import plot_dk
from func_plot_desikan2 import _get_cmap_
from custom_functions import add_figures_to_presentation, create_directory_if_not_exists, create_subplot

os.chdir('/Users/clarariegis/Desktop/Cambridge/mphil_project')



#%% Average across participants visualisation


harmonizations = ['combat'] #, 'no_harmo'
scales = ['scaled'] # , 'no_scaled'
dx_on_off = ['dx', 'no_dx']
list_cond = ['AD', 'ADHD', 'ASD', 'BIP', 'MDD', 'SCZ', 'SUD']
alpha = 0.05
#alpha = 0.002


demographics = pd.read_csv("ABCD/abcd_format_data/demographics.csv")


for harmo in harmonizations: 
    
    for scale in scales: 
        
        for dx in dx_on_off:
        
            dir_path = f"Figures/mixed_model/{harmo}/{scale}/{dx}"
            
            # Results of the longitudinal analysis.
            mlm_res0 = pd.read_csv(dir_path + "/mlm_results.csv")
            
            
            
            mlm_res = mlm_res0.loc[mlm_res0['mlm_pval_int'] < alpha]
            max_plots = mlm_res['condition'].value_counts().max()
            pheno = (np.unique(mlm_res['measure']))
            n_pheno = len(pheno)
            
            # Get the regiong names as a list. 
            roi_names = open('ABCD/abcd_format_data/roi_names.txt').read().split("\n")[:-1] #+['global']

            roi_list = pd.DataFrame(roi_names, columns=['roi'])
            
            
            # True measures.
            path_meas = 'ABCD/abcd_format_data/{hamor}/{scale}/*.csv'
            measures = {}
            for filename in glob.glob(path_meas):
                meas = filename.split('/')[-1].split('.')[0]
                measures[meas] = pd.read_csv(filename)
                
           
            
            # Listing the coefficients that will be used to evaluate the fit of the model.
            list_coeff = ['r2', 'slope', 'p-value']
            # Listing the x-limits for the raincloud plots. 
            list_lims = [4000, 0.05]
            
            
            # Create figures per PGS condition. 
            for condition in list_cond:
                
                cmap  = 'coolwarm'
                                
                # Select the results that correspond to one GWAS condition. 
                cond_res = mlm_res.loc[mlm_res['condition'] == condition]
                
                if len(cond_res) > 0:
                                        
                   
                    # +1 to create a place for the colorbar
                    fig, axs = create_subplot((n_pheno)+1, max_col = 5)
                    axs = axs.flat
                    
                    i = 0
                    cmap = 'Spectral'
                    
                    # Plot the slopes onto the desikan atlas. 
                    for ax in axs:
                
                        if i < len(axs)-1:
                            
                            # The measure currently looping through. 
                            t_pheno = pheno[i]
                            
                            # From the results df, we extract the rows that correspond to 
                            # the current measure. 
                            match_rows = cond_res.loc[cond_res['measure'] == t_pheno]
                            
                            # We merge the dataframes so that we have a result for all ROIs, 
                            # (the non-existent ones will be marked as nans).
                            dsk_df = pd.merge(roi_list, match_rows, how = 'outer', on = 'roi')
                            
                            # Save how many ROIs were significant to add onto the plot.
                            n_sign = sum(cond_res['measure'] == t_pheno)
                            
                            # That's a way to check whether the significant 
                            # is measure the global measure or one of the regions. 
                            if len(dsk_df) > len(roi_list): 
                                print('hello')
                                title = f'''N = {n_sign}, global = {dsk_df.iloc[-1]['estimate']}'''
                                
                                # drop last row:
                                dsk_df = dsk_df.drop(dsk_df.index[-1])
                                
                            else:
                            
                                title = f'N = {n_sign}'
                                
                            
                            # Make it into a directory so that it can be plotted. 
                            dictionary_slopes = dsk_df.set_index('roi')['mlm_est_age'].to_dict()
                            print(dictionary_slopes)
                            
                            
                            # Determine this min and max values of the colorbar. 
                            if abs(min(cond_res['mlm_est_age'])) > max(cond_res['mlm_est_age']):
                                cbar_lim = abs(min(cond_res['mlm_est_age']))
                            else: cbar_lim = max(cond_res['mlm_est_age'])
                            
                            plot_dk(dictionary_slopes, ax = ax, figsize=(20,20),
                                               cmap=cmap, 
                                               background='w', edgecolor='k', 
                                               bordercolor='k', ylabel='', 
                                               vminmax=[-cbar_lim, cbar_lim],
                                               fontsize = 6,
                                               title= title)
                            
                            
                            
                            ax.set_xlabel(f'{t_pheno}', size = 9)
                            # Remove borders or frame of the axis
                            ax.spines['top'].set_visible(False)
                            ax.spines['bottom'].set_visible(False)
                            ax.spines['left'].set_visible(False)
                            ax.spines['right'].set_visible(False)
                            
                            i += 1
                            print(i)
                                
                                
                    # Figure title.
                    fig.suptitle(f'{condition}, {scale}, {harmo}', fontsize=10, y=1.1) 
                    
                    # Adjust the position, colors, and add the colorbar.
                    cmap, norm = _get_cmap_(cmap, cond_res['mlm_est_age'], 
                                            vminmax=[-cbar_lim, cbar_lim])
                    
                    
                    left, bottom, width, height = 0.8, 0.1, 0.015, 0.7  # Example dimensions
                    axs[-1].set_position([left, bottom, width, height])
                    
                    cbar = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                                                    cax = axs[-1],
                                                    orientation='vertical',
                                                    ticklocation='left')
                    
                    t = cbar.ax.yaxis.get_offset_text() 
                    t.set_size(5)
                    cbar.ax.tick_params(labelsize = 10, size=1)
                    cbar.set_label('significant slopes', size = 10)
                    cbar.ax.yaxis.set_label_position('right')
                    
                    # path2 = f'{dir_path}/mlm_dsk_{condition}.png'
                    # plt.savefig(path2, bbox_inches='tight')
                
                
                    
                    
                    # _________________________________________________________________________
                    
                    
                    # # this could be completely out of the loop but just to 
                    # # make sure there is no clash/overwritting of the vars
                    # input_ppt_path = 'Figures/Slides/mphil_figures_clean2.pptx'
                    # scale_ppt_path = 'Figures/Slides/scale.pptx'
                    # output_ppt_path = 'Figures/Slides/mphil_figures_clean2.pptx'
                    
                    # # Desired width and height for the images (in inches)
                    # desired_width = [Inches(7)]
                    # desired_height = [Inches(7)]
                    # left = [Inches(0)] # Adjust position as needed
                    # top = [Inches(0)]
        
         
                   
                    # # every path contains a different type of figure, for the same data. 
                    # # we want them to be plotted in the same slide. 
                    # list_fig_paths = [[path2]]
                    
                    
                    
                    
                    # add_figures_to_presentation(input_ppt_path = input_ppt_path, 
                    #                             scale_ppt_path = scale_ppt_path, 
                    #                             fig_paths = list_fig_paths, 
                    #                             desired_widths = desired_width, 
                    #                             desired_heights = desired_height, 
                    #                             left_positions = left, 
                    #                             top_positions = top, 
                    #                             output_ppt_path = output_ppt_path)
                
                    
            
            
        
        


#%%



#__________________________________________________________________
            
# Raincloud plot 

z = 0
list_res = ['mlm_est_int', 'mlm_pval_int']
labels = ['estimate', 'p-value']
list_lims = [0.0001, 0.05]
pvals = {}

for result in list_res: 
        
    df_dict = {}
    for filename1 in glob.glob("Figures/mixed_model/*/scaled/***/mlm_results.csv"):
        print(filename1)
        label = f"{filename1.split('/')[2]}, {filename1.split('/')[4]}"
        df_dict[label] = pd.read_csv(filename1)[result]
    
        pvals[label] = pd.read_csv(filename1)['mlm_pval_int']
    
    random_colors = sns.cubehelix_palette(start=.5, rot=-.75, as_cmap=False,
                                          n_colors = 4)
    
    
    
    
    rain_font = 25
    alpha = 0.05
    
    # Prepare the color dictionaries. 
    scatt_dict = df_dict.copy()
    mark_dict = df_dict.copy()
    # Now switch the values to a maker or a color, depending on whether the results
    # are significant or not. 
    for key, values in df_dict.items():
        mark_dict[key] = ['o' if item>alpha else '*' for item in pvals[key]]
        scatt_dict[key] = ['black' if item>alpha else 'red' for item in pvals[key]]
    
        
        
    rain_colors = {'box_col': random_colors, 
                    'violin_col': random_colors, 
                    'scatt_col': scatt_dict, 
                    'scatt_marker': mark_dict,
                    'mark_size': 13}
    
    
        
        
    data_x = df_dict.copy()
    
    
    
    fig, ax = plt.subplots(figsize=(11, 10), dpi=800)
    
    
    raincloud2(data_x = df_dict, 
               x_label = labels[z], 
               title = (' '),
               colors = True, 
               ax = ax, 
               size = rain_font, 
               style = rain_colors)
    
    
    
    plt.xlim(min(value for values in df_dict.values() for value in values) - list_lims[z], 
             max(value for values in df_dict.values() for value in values) + list_lims[z])
    
   
    plt.xticks(rotation=45) 
    
    fig_path = f'Figures/mixed_model/sensitivity_{result}_rain.png'
    plt.savefig(fig_path, bbox_inches='tight')  # Save as PNG format
    

    
    # # this could be completely out of the loop but just to 
    # # make sure there is no clash/overwritting of the vars
    # input_ppt_path = 'Figures/Slides/mphil_figures_clean2.pptx'
    # scale_ppt_path = 'Figures/Slides/scale.pptx'
    # output_ppt_path = 'Figures/Slides/mphil_figures_clean2.pptx'
    
    # # Desired width and height for the images (in inches)
    # desired_width = [Inches(7)]
    # desired_height = [Inches(7)]
    # left = [Inches(0)] # Adjust position as needed
    # top = [Inches(0)]
    
    
    
    # # every path contains a different type of figure, for the same data. 
    # # we want them to be plotted in the same slide. 
    # list_fig_paths = [[fig_path]]
    
    
    
    
    # add_figures_to_presentation(input_ppt_path = input_ppt_path, 
    #                             scale_ppt_path = scale_ppt_path, 
    #                             fig_paths = list_fig_paths, 
    #                             desired_widths = desired_width, 
    #                             desired_heights = desired_height, 
    #                             left_positions = left, 
    #                             top_positions = top, 
    #                             output_ppt_path = output_ppt_path)
    
    z += 1
    
    


