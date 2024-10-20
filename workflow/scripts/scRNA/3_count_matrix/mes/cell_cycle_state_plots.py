# Attempt to determine the cluster IDs 
# Build the reference if possible

import pandas as pd
from logging_func import handle_exception, setup_logger
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
matplotlib.use('Agg')

#########################################################3
# LOAD RAW COUNT MATRIX
#########################################################3

adata = sc.read_h5ad(snakemake.input.adata)

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

logger = setup_logger(snakemake.log[0])

dim_red_plots = snakemake.params.dim_red_plots
cluster_plots = snakemake.params.cluster_plots
cluster_an = snakemake.params.cluster_an

dim_red_plot = snakemake.output.dim_red_plot
violin_plot = snakemake.output.violin_plot

#########################################################3
# PLOTS
#########################################################3

# Embedding plots (preannot)
fig, ax = plt.subplots(figsize=(16,10), dpi=300, ncols=3, nrows=2)

sc.pl.embedding(adata, basis = dim_red_plots, color = cluster_plots, ax = ax[0,0], frameon=False, legend_loc='right margin')
ax[0,0].set_xticklabels(ax[0,0].get_xticklabels(), rotation=90)

sc.pl.embedding(adata, basis = dim_red_plots, color = cluster_an, ax = ax[0,1], frameon=False, legend_loc='right margin')
ax[0,1].set_xticklabels(ax[0,1].get_xticklabels(), rotation=90)

sc.pl.embedding(adata, basis = dim_red_plots, color='cell_type', ax = ax[0,2], frameon=False, legend_loc='right margin')
ax[0,2].set_xticklabels(ax[0,2].get_xticklabels(), rotation=90)

sc.pl.embedding(adata, basis = dim_red_plots, color='S_score', ax = ax[1,0], frameon=False, legend_loc='right margin')
ax[1,0].set_xticklabels(ax[1,0].get_xticklabels(), rotation=90)

sc.pl.embedding(adata, basis = dim_red_plots, color='G2M_score', ax = ax[1,1], frameon=False, legend_loc='right margin')
ax[1,1].set_xticklabels(ax[1,1].get_xticklabels(), rotation=90)

sc.pl.embedding(adata, basis = dim_red_plots, color = 'phase', ax=ax[1,2], frameon=False, legend_loc='right margin')
ax[1,2].set_xticklabels(ax[1,2].get_xticklabels(), rotation=90)

#fig.subplots_adjust(wspace=0.2, hspace=0.2)
plt.tight_layout()
fig.savefig(dim_red_plot)


ncols_violin = 1
nrows_violin = 2
# VIOLIN GENE PLOTS
fig, ax = plt.subplots(figsize=(16,10), dpi=300, ncols = ncols_violin, nrows = nrows_violin) 

for row in range(nrows_violin):
    for col in range(ncols_violin):
        if row == 0:
            sc.pl.violin(adata, keys = "S_score", groupby = "sample", size = 0, jitter = 0.4, ax = ax[row])
        else:
            sc.pl.violin(adata, keys = "G2M_score", groupby = "sample", size = 0, jitter = 0.4, ax = ax[row])      
        
plt.tight_layout()
fig.savefig(violin_plot)