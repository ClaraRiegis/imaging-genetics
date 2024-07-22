#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 14:24:07 2024

@author: clarariegis
"""

#%%


import pandas as pd
import matplotlib.pyplot as plt
import os
from pptx import Presentation
from pptx.util import Inches


os.chdir('/Users/clarariegis/Desktop/Cambridge/mphil_project/code/scripts')
from custom_functions import add_figures_to_presentation, create_directory_if_not_exists



#%%
# The genetic correlation matrix was calculated in Bash using LDSC

os.chdir('/Users/clarariegis/Desktop/Cambridge/mphil_project')
cmap = 'coolwarm'

dir_path = 'Figures/correlations'
create_directory_if_not_exists(dir_path)

df = pd.read_csv(f'{dir_path}/genetic_correlation_results.csv',
                 index_col=None)[["p1","p2","rg"]]
matrix = df.pivot(index='p1', columns='p2', values='rg')
matrix = matrix.combine_first(matrix.T)



fig, ax = plt.subplots(dpi=700)
plt.rcParams["figure.figsize"] = (15,15)
plt.rcParams.update({'font.size': 20})

im = ax.imshow(matrix, cmap=cmap, vmin = -1, vmax = 1)
ax.set_title('Genetic correlation')
fig.colorbar(im, ticks=[-1,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0])

# Set the tick labels
ax.set_xticks(range(len(matrix.columns)))
ax.set_xticklabels(matrix.columns, rotation=45, ha='left')
ax.set_yticks(range(len(matrix.index)))
ax.set_yticklabels(matrix.index)

fig_path = (f'{dir_path}/genetic_correlation_matrix.png')
plt.savefig(fig_path)  # Save as PNG format

matrix = matrix.fillna(1)
matrix.to_csv(f'{dir_path}/genetic_corr_mat.csv', index=False)



#%%


input_ppt_path = 'Figures/Slides/mphil_figures_clean2.pptx'
scale_ppt_path = 'Figures/Slides/scale.pptx'
output_ppt_path = 'Figures/Slides/mphil_figures_clean2.pptx'



# every path contains a different type of figure, for the same data. 
# we want them to be plotted in the same slide. 
list_fig_paths = [[fig_path]]

# Desired width and height for the images (in inches)
desired_width = [Inches(7)]
desired_height = [Inches(7)]
left = [Inches(0)] # Adjust position as needed
top = [Inches(1)]



add_figures_to_presentation(input_ppt_path = input_ppt_path, 
                            scale_ppt_path = scale_ppt_path, 
                            fig_paths = list_fig_paths, 
                            desired_widths = desired_width, 
                            desired_heights = desired_height, 
                            left_positions = left, 
                            top_positions = top, 
                            output_ppt_path = output_ppt_path)





