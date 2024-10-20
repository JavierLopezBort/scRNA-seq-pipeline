#########################################################3
# IMPORT PACKAGES
#########################################################3

import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

#########################################################3
# IMPORT FROM SNAKEMAKE
#########################################################3

adata_preannot = sc.read_h5ad(snakemake.input.adata_preannot)
adata_annot = sc.read_h5ad(snakemake.input.adata_annot)

dim_red_annot_plot = snakemake.output.dim_red_annot_plot
dim_red_preannot_plot = snakemake.output.dim_red_preannot_plot

#########################################################3
# PLOTTING
#########################################################3

# Embedding plots (preannot)
fig, ax = plt.subplots(figsize=(25,10), dpi=300, ncols=4, nrows=2)

sc.pl.embedding(adata_preannot, basis = 'X_emb_2', color = "cluster", ax = ax[0,0], frameon=False, legend_loc='right margin')
ax[0,0].set_xticklabels(ax[0,0].get_xticklabels(), rotation=90)

sc.pl.embedding(adata_preannot, basis = 'X_emb_2', color = "study", ax = ax[0,1], frameon=False, legend_loc='right margin')
ax[0,1].set_xticklabels(ax[0,1].get_xticklabels(), rotation=90)

sc.pl.embedding(adata_preannot, basis = 'X_emb_2', color='gender', ax = ax[0,2], frameon=False, legend_loc='right margin')
ax[0,2].set_xticklabels(ax[0,2].get_xticklabels(), rotation=90)

sc.pl.embedding(adata_preannot, basis = 'X_emb_2', color='sample', ax = ax[0,3], frameon=False, legend_loc='right margin')
ax[0,3].set_xticklabels(ax[0,3].get_xticklabels(), rotation=90)

sc.pl.embedding(adata_preannot, basis = 'X_emb_2', color = 'specie', ax = ax[1,0], frameon=False, legend_loc='right margin')
ax[1,0].set_xticklabels(ax[1,0].get_xticklabels(), rotation=90)

sc.pl.embedding(adata_preannot, basis = 'X_emb_2', color = 'cell_type', ax=ax[1,1], frameon=False, legend_loc='right margin')
ax[1,1].set_xticklabels(ax[1,1].get_xticklabels(), rotation=90)

sc.pl.embedding(adata_preannot, basis = 'X_emb_2', color = 'cell_annot_comb', ax=ax[1,2], frameon=False, legend_loc='right margin')
ax[1,2].set_xticklabels(ax[1,2].get_xticklabels(), rotation=90)

sc.pl.embedding(adata_preannot, basis = 'X_emb_2', color = 'cell_annot', ax=ax[1,3], frameon=False, legend_loc='right margin')
ax[1,3].set_xticklabels(ax[1,3].get_xticklabels(), rotation=90)

#fig.subplots_adjust(wspace=0.2, hspace=0.2)
plt.tight_layout()
fig.savefig(dim_red_preannot_plot)

# Embedding plots (preannot)
fig, ax = plt.subplots(figsize=(25,10), dpi=300, ncols=4, nrows=2)

sc.pl.embedding(adata_annot, basis = 'X_emb_2', color = "cluster", ax = ax[0,0], frameon=False, legend_loc='right margin')
ax[0,0].set_xticklabels(ax[0,0].get_xticklabels(), rotation=90)

sc.pl.embedding(adata_annot, basis = 'X_emb_2', color = "study", ax = ax[0,1], frameon=False, legend_loc='right margin')
ax[0,1].set_xticklabels(ax[0,1].get_xticklabels(), rotation=90)

sc.pl.embedding(adata_annot, basis = 'X_emb_2', color='gender', ax = ax[0,2], frameon=False, legend_loc='right margin')
ax[0,2].set_xticklabels(ax[0,2].get_xticklabels(), rotation=90)

sc.pl.embedding(adata_annot, basis = 'X_emb_2', color='sample', ax = ax[0,3], frameon=False, legend_loc='right margin')
ax[0,3].set_xticklabels(ax[0,3].get_xticklabels(), rotation=90)

sc.pl.embedding(adata_annot, basis = 'X_emb_2', color = 'specie', ax = ax[1,0], frameon=False, legend_loc='right margin')
ax[1,0].set_xticklabels(ax[1,0].get_xticklabels(), rotation=90)

sc.pl.embedding(adata_annot, basis = 'X_emb_2', color = 'cell_type', ax=ax[1,1], frameon=False, legend_loc='right margin')
ax[1,1].set_xticklabels(ax[1,1].get_xticklabels(), rotation=90)

sc.pl.embedding(adata_annot, basis = 'X_emb_2', color = 'cell_annot_comb', ax=ax[1,2], frameon=False, legend_loc='right margin')
ax[1,2].set_xticklabels(ax[1,2].get_xticklabels(), rotation=90)

sc.pl.embedding(adata_annot, basis = 'X_emb_2', color = 'cell_annot', ax=ax[1,3], frameon=False, legend_loc='right margin')
ax[1,3].set_xticklabels(ax[1,3].get_xticklabels(), rotation=90)

#fig.subplots_adjust(wspace=0.2, hspace=0.2)
plt.tight_layout()
fig.savefig(dim_red_annot_plot)