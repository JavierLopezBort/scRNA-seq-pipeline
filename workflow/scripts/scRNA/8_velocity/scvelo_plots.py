import scanpy as sc
import matplotlib.pyplot as plt
import louvain
import scvelo as scv

#########################################################3
# LOAD RAW COUNT MATRIX
#########################################################3

adata = sc.read_h5ad(snakemake.input.adata)

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

# Output plots
velocity_plot = snakemake.output.velocity_plot

#########################################################3
# PLOTS
#########################################################3

# Velocity plot
if "cell_annot_comb" in adata.obs.columns:
    
    # Dimensionality reduction
    fig, ax = plt.subplots(figsize=(25,10), dpi=300, ncols=4, nrows=2)

    scv.pl.velocity_embedding_stream(adata, basis = "X_emb_2", color = "cluster", ax = ax[0,0], frameon=False, legend_loc='on data')
    ax[0,0].set_xticklabels(ax[0,0].get_xticklabels(), rotation=90)

    scv.pl.velocity_embedding_stream(adata, basis = "X_emb_2", color='study', ax = ax[0,1], frameon=False, legend_loc='right margin')
    ax[0,1].set_xticklabels(ax[0,1].get_xticklabels(), rotation=90)

    scv.pl.velocity_embedding_stream(adata, basis = "X_emb_2", color='gender', ax = ax[0,2], frameon=False, legend_loc='right margin')
    ax[0,2].set_xticklabels(ax[0,2].get_xticklabels(), rotation=90)

    scv.pl.velocity_embedding_stream(adata, basis = "X_emb_2", color='sample', ax = ax[0,3], frameon=False, legend_loc='right margin')
    ax[0,3].set_xticklabels(ax[0,3].get_xticklabels(), rotation=90)

    scv.pl.velocity_embedding_stream(adata, basis = "X_emb_2", color = 'specie', ax=ax[1,0], frameon=False, legend_loc='right margin')
    ax[1,0].set_xticklabels(ax[1,0].get_xticklabels(), rotation=90)

    scv.pl.velocity_embedding_stream(adata, basis = "X_emb_2", color = 'cell_type', ax=ax[1,1], frameon=False, legend_loc='right margin')
    ax[1,1].set_xticklabels(ax[1,1].get_xticklabels(), rotation=90)

    scv.pl.velocity_embedding_stream(adata, basis = "X_emb_2", color = 'cell_annot_comb', ax=ax[1,2], frameon=False, legend_loc='right margin')
    ax[1,2].set_xticklabels(ax[1,2].get_xticklabels(), rotation=90)

    scv.pl.velocity_embedding_stream(adata, basis = "X_emb_2", color = 'cell_annot', ax=ax[1,3], frameon=False, legend_loc='right margin')
    ax[1,3].set_xticklabels(ax[1,3].get_xticklabels(), rotation=90)

    #fig.subplots_adjust(wspace=0.2, hspace=0.2)
    plt.tight_layout()
    fig.savefig(velocity_plot)

else:
    
    # Dimensionality reduction
    fig, ax = plt.subplots(figsize=(20,10), dpi=300, ncols=3, nrows=2)

    scv.pl.velocity_embedding_stream(adata, basis = "X_emb_2", color = "cluster", ax = ax[0,0], frameon=False, legend_loc='on data')
    ax[0,0].set_xticklabels(ax[0,0].get_xticklabels(), rotation=90)

    scv.pl.velocity_embedding_stream(adata, basis = "X_emb_2", color='study', ax = ax[0,1], frameon=False, legend_loc='right margin')
    ax[0,1].set_xticklabels(ax[0,1].get_xticklabels(), rotation=90)

    scv.pl.velocity_embedding_stream(adata, basis = "X_emb_2", color='gender', ax = ax[0,2], frameon=False, legend_loc='right margin')
    ax[0,2].set_xticklabels(ax[0,2].get_xticklabels(), rotation=90)

    scv.pl.velocity_embedding_stream(adata, basis = "X_emb_2", color='sample', ax = ax[1,0], frameon=False, legend_loc='right margin')
    ax[1,0].set_xticklabels(ax[1,0].get_xticklabels(), rotation=90)

    scv.pl.velocity_embedding_stream(adata, basis = "X_emb_2", color = 'specie', ax=ax[1,1], frameon=False, legend_loc='right margin')
    ax[1,1].set_xticklabels(ax[1,1].get_xticklabels(), rotation=90)

    scv.pl.velocity_embedding_stream(adata, basis = "X_emb_2", color = 'cell_type', ax=ax[1,2], frameon=False, legend_loc='right margin')
    ax[1,2].set_xticklabels(ax[1,2].get_xticklabels(), rotation=90)

    #fig.subplots_adjust(wspace=0.2, hspace=0.2)
    plt.tight_layout()
    fig.savefig(velocity_plot)

