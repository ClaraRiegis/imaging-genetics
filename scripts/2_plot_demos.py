#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 12:36:12 2024

Author: clarariegis

Purpose: Visualise the demographic data of the 5th release of the ABCD dataset.

Description: Plot demographics
             - create a power point presentation
             - Plot of the gender, age and diagnosis distributions
             - Table counts to see how many pcp are in all follow ups. 
    
"""

#%%
import pandas as pd
import numpy as np
import os
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from docx import Document
from pptx.util import Inches

os.chdir('/Users/clarariegis/Desktop/Cambridge/mphil_project/code/scripts')
from custom_functions import add_figures_to_presentation, create_directory_if_not_exists


os.chdir('/Users/clarariegis/Desktop/Cambridge/mphil_project')

#%% Demographics plots


# Load the demographic data:
demographics = pd.read_csv("ABCD/abcd_format_data/demographics.csv")
# Check the different time points availables.
t_points = np.unique(demographics['eventname'])
# Reorder the time points such that they are in chronological order.
time_points = [t_points[2], t_points[0], t_points[1]]
# Covariates:
covars = pd.read_csv("ABCD/abcd_format_data/covar.csv")


# Create a figure
size_plt = 15
fig, ax = plt.subplot_mosaic("""ABC""", figsize=(25, 7))   
# Chosen RBG colors (they need to be scaled so that 0<color<1)
colors = [(37/255,75/255,124/255), 
          (59/255,104/255,161/255), 
          (105/255,145/255,195/255)]
plt.subplots_adjust(hspace=0.6, bottom=0.1, top=0.7, right=0.8)
patches = []
demo_table = pd.DataFrame(columns = ["Female", "Male", "Min age", "Max age", 
                                     "Mean age", "ASD", "CN", "IDD", "Dev.", 
                                     "Psych.", "SCZ"])
i=0

for time in time_points:
    
    time_demos = demographics.loc[demographics['eventname'] == time]
    n_pcp = len(time_demos)
    
    # - Sex 
    pcp_gender = time_demos['sex']
    n_fem = sum(pcp_gender == "Female")
    n_mal = sum(pcp_gender == "Male")
    
    # - Age 
    pcp_ages = time_demos['interview_age']
    min_age = min(pcp_ages)
    max_age = max(pcp_ages)
    mean_age = np.mean(pcp_ages)
    
    # - Diagnosis
    pcp_diag = time_demos['dx']
    unique_diag = np.unique(pcp_diag)
    n_asd = sum(pcp_diag == "ASD")
    n_cn = sum(pcp_diag == "CN")
    n_idd = sum(pcp_diag == "IDD")
    n_dev = sum(pcp_diag == "Other Developmental")
    n_psy = sum(pcp_diag == "Other Psychiatric")
    n_scz = sum(pcp_diag == "SCZ")
    
    # Store these descriptive statistics into a dictionary. 
    dic = {"Female": n_fem, "Male": n_mal, "Min age": min_age, 
            "Max age": max_age, "Mean age": mean_age, "ASD": n_asd, "CN": n_cn, 
            "IDD": n_idd, "Dev.": n_dev, "Psych.": n_psy, "SCZ": n_scz}
    
    # Convert the dictionary into a dataframe and append it for each t-point. 
    demo_table = pd.concat([demo_table, pd.DataFrame(dic, index=[0])],
                           ignore_index=True)
    
    # Plot the sex, age and diagnosis distributions:
    ax["A"].bar(x = [-0.9+i/2.5, 0.9+i/2.5], height = [n_fem, n_mal], 
                    color= colors[i], width = 0.4)
    ax["B"].hist(pcp_ages, bins=50, color = colors[i], alpha=0.7)
    ax['C'].bar(x = [-9+i/2, -6+i/2, -3+i/2, 0+i/2, 3+i/2, 6+i/2],
                height = [n_cn, n_dev, n_asd, n_psy, n_idd, n_scz], 
                color= colors[i], width = 0.5)
        
    # Store the patches and mumber of participants for the legend. 
    patches += [mpatches.Patch(color = colors[i], 
                               label= f't{i} (n = {n_pcp})')]

    
    
    # Add titles and name axes.
    axes = ["A", "B", "C"]
    x_labels = ('Gender', 'Age', 'Diagnosis')
    if x_labels[i] == 'Age': unit = " (mo.) "
    else: unit = " "
    ax[axes[i]].set_xlabel((x_labels[i] + unit ), size = size_plt)
    ax[axes[i]].set_ylabel('Frequency', size = size_plt)
    ax[axes[i]].set_title((x_labels[i] + ' distribution') , 
                     size = (size_plt + 2))
    
        
    i+=1
    
# Show the plot and legends.
fig.legend(handles=patches, loc = (0.79,0.81), fontsize = size_plt) 

# Define custom ticks and labels.
# - Sex
custom_ticks = [-0.9+1/2.5, 0.9+1/2.5]  # tick locations.
custom_labels = ['Females', 'Males']    # tick labels.
ax['A'].set_xticks(custom_ticks)        # Set custom ticks.
ax['A'].set_xticklabels(custom_labels)  # Set custom labels.

# - Diagnosis
custom_ticks = [-9+1/2, -6+1/2, -3+1/2, 0+1/2, 3+1/2, 6+1/2]
custom_labels = ["ASD", "CN", "IDD", "Dev.", "Psych.", "SCZ"]
ax['C'].set_xticks(custom_ticks)        # Set custom ticks.
ax['C'].set_xticklabels(custom_labels)  # Set custom labels.


# If the file we want to save the output in does not exist, create it: 
file_path = 'Figures/demographics'
create_directory_if_not_exists(file_path)
fig_path1 = [file_path + '/demographic_plot.png']
plt.savefig(fig_path1[0])  # Save as PNG format
    

# Save the table with the descriptive statistics.
fig_path1 += [file_path + '/demographics_table.csv']
demo_table.to_csv(fig_path1[1], index=False)


#%% Count table

# Create a dictionary to store counts for each time point
counts = {tp: len(demographics[demographics['eventname'] == tp]
                  ['ID'].unique()) for tp in time_points}

# Create a DataFrame from the counts
count_table = pd.DataFrame.from_dict(counts, orient='index', 
                                     columns=['Participant Count'])

# Add a column for the intersection counts
for tp1 in time_points:
    for tp2 in time_points:
        intersection_count = (len(set(demographics[demographics['eventname'] 
                              == tp1]['ID']) & set(demographics[demographics
                                            ['eventname'] == tp2]['ID'])))
        count_table.at[tp1, f'Intersection with {tp2}'] = intersection_count

print(count_table)



# Save the table with the participant count.
fig_path1 += [file_path + '/count_table.csv']
count_table.to_csv(fig_path1[2], index=False)



#%% save figures in slides

input_ppt_path = 'Figures/Slides/mphil_figures_clean.pptx'
scale_ppt_path = 'Figures/Slides/scale.pptx'
output_ppt_path = 'Figures/Slides/mphil_figures_clean2.pptx'



# every path contains a different type of figure, for the same data. 
# we want them to be plotted in the same slide. 
fig_paths = [fig_path1]

# Desired width and height for the images (in inches)
desired_width = [Inches(7)]
desired_height = [Inches(5)]
left = [Inches(0)] # Adjust position as needed
top = [Inches(0.5)]



add_figures_to_presentation(input_ppt_path = input_ppt_path, 
                            scale_ppt_path = scale_ppt_path, 
                            fig_paths = fig_paths, 
                            desired_widths = desired_width, 
                            desired_heights = desired_height, 
                            left_positions = left, 
                            top_positions = top, 
                            output_ppt_path = output_ppt_path, 
                            num_format='int')


