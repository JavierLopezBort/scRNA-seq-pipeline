# Attempt to determine the cluster IDs 
# Build the reference if possible

import pandas as pd
from logging_func import handle_exception, setup_logger
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import json
matplotlib.use('Agg')

#########################################################3
# LOAD RAW COUNT MATRIX
#########################################################3

adata = sc.read_h5ad(snakemake.input.adata)

with open(snakemake.input.data_list_perc, 'r') as json_file:
    data_list_perc_pre = json.load(json_file)

with open(snakemake.input.data_list_numb, 'r') as json_file:
    data_list_numb_pre = json.load(json_file)

"""
# CHANGE ORDER
new_order_indices = [0, 6, 1, 7, 2, 8, 3, 9, 4, 10, 5, 11]  # This defines the new order for the elements
# Using list comprehension to create a new list based on the order of indices
data_list_perc = [data_list_perc_pre[i] for i in new_order_indices]
data_list_numb = [data_list_numb_pre[i] for i in new_order_indices]   
"""

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

logger = setup_logger(snakemake.log[0])

dim_red_plot = snakemake.output.dim_red_plot
barcodes_counts_perc_plot = snakemake.output.barcodes_counts_perc_plot
barcodes_counts_numb_plot = snakemake.output.barcodes_counts_numb_plot

#########################################################3
# PLOTS
#########################################################3

# Embedding plots
fig, ax = plt.subplots(figsize=(16,10), dpi=300, ncols=1)
sc.pl.embedding(adata, basis = "X_emb_2", color = "sample", ax = ax[0], frameon=False, legend_loc='right margin', show = False)
plt.tight_layout()
fig.savefig(dim_red_plot)

# BARCODE PERC PLOT
nrows = 2
ncols = 3
fig, ax = plt.subplots(figsize=(20,10), dpi=300, ncols=ncols, nrows=nrows)

i = 0
for row in range(nrows):
    ax[row,0].set_ylabel('% transfection')
    for col in range(ncols):
        if i < len(data_list_perc): 
            sample = data_list_perc[i][0]
            sample_counts = pd.Series(data_list_perc[i][2], index=data_list_perc[i][1])
            bars = sample_counts.plot(kind='bar', ax=ax[row,col])
            ax[row,col].tick_params(axis='x', rotation=45)
            ax[row,col].set_title(sample)
            ax[row,col].set_ylim(0, 35)

            # Adding text above each bar
            for bar in bars.containers[0]:
                height = bar.get_height()
                ax[row, col].text(bar.get_x() + bar.get_width() / 2., height, 
                                    f'{height}', ha='center', va='bottom')
        
            i += 1 
        else:
            fig.delaxes(ax[row,col])

plt.tight_layout()
#fig.subplots_adjust(wspace=0.25, hspace=0)
fig.savefig(barcodes_counts_perc_plot)

# BARCODE NUMB PLOT
nrows = 2
ncols = 3
fig, ax = plt.subplots(figsize=(20,10), dpi=300, ncols=ncols, nrows=nrows)

i = 0
for row in range(nrows):
    ax[row,0].set_ylabel('N cells')
    for col in range(ncols):
        if i < len(data_list_numb): 
            sample = data_list_numb[i][0]
            sample_counts = pd.Series(data_list_numb[i][2], index=data_list_numb[i][1])
            bars = sample_counts.plot(kind='bar', ax=ax[row,col])
            ax[row,col].tick_params(axis='x', rotation=45)
            ax[row,col].set_title(sample)
            ax[row,col].set_ylim(0, 400)

            # Adding text above each bar
            for bar in bars.containers[0]:
                height = bar.get_height()
                ax[row, col].text(bar.get_x() + bar.get_width() / 2., height, 
                                    f'{height}', ha='center', va='bottom')
        
            i += 1 
        else:
            fig.delaxes(ax[row,col])

plt.tight_layout()
#fig.subplots_adjust(wspace=0.25, hspace=0)
fig.savefig(barcodes_counts_numb_plot)
