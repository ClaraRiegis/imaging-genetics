#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  8 14:32:27 2024

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
import matplotlib.colors as colors
from pptx.util import Inches



os.chdir('/Users/clarariegis/Desktop/Cambridge/mphil_project/code/scripts')
from raincloud2 import raincloud2
# import func_plot_desikan2 
from func_plot_desikan2 import plot_dk
from func_plot_desikan2 import _get_cmap_
from custom_functions import add_figures_to_presentation, create_directory_if_not_exists, create_subplot


os.chdir('/Users/clarariegis/Desktop/Cambridge/mphil_project')



#%% 


demographics = pd.read_csv("ABCD/abcd_format_data/demographics.csv")

    

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
           
roi_list = list(data_dk.keys()) + ['global']

# True measures.
path_meas = 'ABCD/abcd_format_data/combat/scaled/*.csv'
measures = {}
self_corr = {}

roi_avr = pd.DataFrame()
mean_corr = {}
measures2 = {}
meas_list = []
cmap = 'coolwarm'

avoid = ['vol-sulc', 'thk-sulc', 'thk-vol']

i = 0
fig, axs = create_subplot(12, max_col = 3)
for filename1 in glob.glob(path_meas):
    meas = filename1.split('/')[-1].split('.')[0]
    meas_list += [meas]
    print(meas)
    measures[meas] = pd.read_csv(filename1).sort_values(by='ID', ascending=False)[roi_list]
    self_corr[meas] = np.corrcoef(measures[meas])
    
    roi_avr[meas] = measures[meas].mean(axis=1)
    mean_corr = np.corrcoef(roi_avr)
    
    
        
    for filename2 in glob.glob(path_meas):
        
        meas2 = filename2.split('/')[-1].split('.')[0]
        measures2[meas2] = pd.read_csv(filename2)[roi_list]
        print(meas2)
        roi_corr = {}
        comb = f'{meas}-{meas2}'
        if filename1 != filename2 and meas != 'area':
            
            for roi in roi_list[:-1]:
                roi_corr[roi] = measures[meas][roi].corr(measures2[meas2][roi])
            
            
            if comb not in avoid:
                plot_dk(roi_corr, ax = axs[i], figsize=(20,20),
                                   cmap=cmap, 
                                   background='w', edgecolor='k', 
                                   bordercolor='black', ylabel='', 
                                   vminmax=[-1, 1],
                                   fontsize = 5,
                                   title= f'{comb}')
                
                
            # Remove borders or frame of the axis
            axs[i].spines['top'].set_visible(False)
            axs[i].spines['bottom'].set_visible(False)
            axs[i].spines['left'].set_visible(False)
            axs[i].spines['right'].set_visible(False)
            
            axs[i].axis('off')
            
            i += 1
            
cmap = 'coolwarm'
# Adjust the position, colors, and add the colorbar.
cmap, norm = _get_cmap_(cmap, values = 0 , vminmax=[-1, 1])


left, bottom, width, height = 0.8, 0.1, 0.015, 0.2  # Example dimensions
axs[-1].set_position([left, bottom, width, height])

cbar = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                                cax = axs[-1],
                                orientation='vertical',
                                ticklocation='right')

t = cbar.ax.yaxis.get_offset_text() 
t.set_size(5)
cbar.ax.tick_params(labelsize = 5, size=1)
cbar.set_label('pearson correlation', size = 5)
cbar.ax.yaxis.set_label_position('left')


fig_path = ('Figures/correlations/phenotypic_correlation_dsk.png')
plt.savefig(fig_path)  # Save as PNG format



#%% 

merged = pd.DataFrame()
list_measures = []
for filename in glob.glob(path_meas):
    meas = filename.split('/')[-1].split('.')[0]
    print(meas)
    data = pd.read_csv(filename).sort_values(by='ID', ascending=False)[roi_list]
    
    data.columns = [f'{meas}_' + col for col in data.columns]
    
    merged = pd.concat([data, merged], axis = 1)
    list_measures += [meas]
    
from mpl_toolkits.axes_grid1 import make_axes_locatable

# the order of the measure names needs to be reversed to match the plot:
list_measures.reverse()

# Correlation matrix
corr_mat = np.corrcoef(merged, rowvar = False)

    

start = 0
step = 69
num_lists = len(meas_list)
labels = np.array([[start + i*step, start + (i+1)*step - 1] for i in range(num_lists)])





def generate_colors(num_colors):
    # Define a colormap with a range of colors
    cmap = plt.get_cmap('tab10')  # You can choose other colormaps like 'viridis', 'jet', etc.

    # Generate a list of colors from the colormap
    col = [cmap(i) for i in range(num_colors)]

    return col


col = generate_colors(len(meas_list))


    
cmap = 'coolwarm'

fig, ax = plt.subplots()
plt.rcParams.update({'font.size': 30})

plt.rcParams["figure.figsize"] = (35,35)

im = ax.imshow(corr_mat, cmap=cmap, vmin = -1, vmax = 1)
ax.set_title('Phenotypic correlation')
fig.colorbar(im, ticks=[-1,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0])



# create axes next to plot
divider = make_axes_locatable(ax)
axb = divider.append_axes("bottom", "5%", pad=0.06, sharex=ax)
axl = divider.append_axes("left", "5%", pad=0.06, sharey=ax)
axb.invert_yaxis()
axl.invert_xaxis()
axb.axis("off")
axl.axis("off")

ax.tick_params(axis='x', colors='w')  # Changes the tick color for the x-axis
ax.tick_params(axis='y', colors='w')  # Changes the tick color for the y-axis


# plot colored bar plots to the axes
barkw = dict(color=['white', 'white', 'white', 'white'], linewidth=0.9, ec="k", clip_on=False, align='edge',) # col
axb.bar(labels[:,0],np.ones(len(labels)), 
        width=np.diff(labels, axis=1).flatten(), **barkw)
axl.barh(labels[:,0],np.ones(len(labels)), 
         height=np.diff(labels, axis=1).flatten(), **barkw)


# Add text labels to the bottom and left bar plots
for idx, label in enumerate(labels):
    axb.text((label[0] + label[1]) / 2, 0.5, f'{list_measures[idx]}', ha='center', va='center', fontsize=20)
    axl.text(0.5, (label[0] + label[1]) / 2, f'{list_measures[idx]}', ha='center', va='center', rotation='vertical', fontsize=20)

plt.show()
fig_path = ('Figures/correlations/phenotypic_correlation_matrix.png')
plt.savefig(fig_path)  # Save as PNG format



pd.DataFrame(corr_mat).to_csv('Figures/correlations/pheno_corr_mat.csv', index=False)
#merged.to_csv('ABCD/abcd_format_data/pheno_corr_mat.csv', index=False)



#%%




input_ppt_path = 'Figures/Slides/mphil_figures_clean2.pptx'
scale_ppt_path = 'Figures/Slides/scale.pptx'
output_ppt_path = 'Figures/Slides/mphil_figures_clean2.pptx'



# every path contains a different type of figure, for the same data. 
# we want them to be plotted in the same slide. 
list_fig_paths = [['Figures/correlations/phenotypic_correlation_matrix.png'],
                  ['Figures/correlations/phenotypic_correlation_dsk.png']]

# Desired width and height for the images (in inches)
desired_width = [Inches(7), Inches(7)]
desired_height = [Inches(7), Inches(7)]
left = [Inches(0), Inches(0)] # Adjust position as needed
top = [Inches(0),Inches(7)]



add_figures_to_presentation(input_ppt_path = input_ppt_path, 
                            scale_ppt_path = scale_ppt_path, 
                            fig_paths = list_fig_paths, 
                            desired_widths = desired_width, 
                            desired_heights = desired_height, 
                            left_positions = left, 
                            top_positions = top, 
                            output_ppt_path = output_ppt_path)














































