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

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

logger = setup_logger(snakemake.log[0])

threshold = snakemake.params.threshold

dim_red_plot = snakemake.output.dim_red_plot

#########################################################3
# GET THE TICKS
#########################################################3

min_values = np.min(adata.obsm["X_emb_2"], axis=0)
max_values = np.max(adata.obsm["X_emb_2"], axis=0)

X_mde_1 = np.linspace(min_values[0], max_values[0], 20)
X_mde_2 = np.linspace(min_values[1], max_values[1], 20)

X_mde_1_str = list()
X_mde_2_str = list()

for i in X_mde_1: X_mde_1_str.append(str(round(i,2)))
for i in X_mde_2: X_mde_2_str.append(str(round(i,2)))

#########################################################3
# PLOTS
#########################################################3

# Embedding plots (preannot)
fig, ax = plt.subplots(figsize=(20,10), dpi=300, ncols=2, nrows=1)

sc.pl.embedding(adata, basis = "X_emb_2", color = "sample", ax = ax[0], frameon=True, legend_loc='right margin', show = False)
ax[0].set_xticks(X_mde_1, labels = X_mde_1_str, rotation=90)
ax[0].set_yticks(X_mde_2, labels = X_mde_2_str)
"""
sc.pl.embedding(adata, basis = "X_emb_2", color = "cell_annot", ax = ax[1], frameon=True, legend_loc='right margin', show = False)
ax[1].set_xticks(X_mde_1, labels = X_mde_1_str, rotation=90)
ax[1].set_yticks(X_mde_2, labels = X_mde_2_str)
ax[1].axvline(x=threshold, color='black', linestyle='--')
"""
sc.pl.embedding(adata, basis = "X_emb_2", color = "cell_annot_comb", ax = ax[1], frameon=True, legend_loc='right margin', show = False)
ax[1].set_xticks(X_mde_1, labels = X_mde_1_str, rotation=90)
ax[1].set_yticks(X_mde_2, labels = X_mde_2_str)
ax[1].axvline(x=threshold, color='black', linestyle='--')

#fig.subplots_adjust(wspace=0.2, hspace=0.2)
plt.tight_layout()
fig.savefig(dim_red_plot)
