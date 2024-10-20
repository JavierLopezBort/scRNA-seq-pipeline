import scanpy as sc
import matplotlib.pyplot as plt
import louvain

#########################################################3
# LOAD RAW COUNT MATRIX
#########################################################3

adata = sc.read_h5ad(snakemake.input.adata)

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

# Output plots
pca_variance_plot = snakemake.output.pca_variance_plot
dim_red_plot = snakemake.output.dim_red_plot
dendrogram_plot = snakemake.output.dendrogram_plot
corr_mat_plot = snakemake.output.corr_mat_plot

#########################################################3
# PLOTS
#########################################################3

# Proportions
sc.pl.pca_variance_ratio(
    adata,
    log=True,
    show = False
)
fig = plt.gcf()
fig.savefig(pca_variance_plot)

if "cell_annot_comb" in adata.obs.columns:
    
    # Dimensionality reduction
    fig, ax = plt.subplots(figsize=(25,10), dpi=300, ncols=4, nrows=2)

    sc.pl.embedding(adata, basis = "X_emb_2", color = "cluster", ax = ax[0,0], frameon=False, legend_loc='on data')
    ax[0,0].set_xticklabels(ax[0,0].get_xticklabels(), rotation=90)

    sc.pl.embedding(adata, basis = "X_emb_2", color='study', ax = ax[0,1], frameon=False, legend_loc='right margin')
    ax[0,1].set_xticklabels(ax[0,1].get_xticklabels(), rotation=90)

    sc.pl.embedding(adata, basis = "X_emb_2", color='gender', ax = ax[0,2], frameon=False, legend_loc='right margin')
    ax[0,2].set_xticklabels(ax[0,2].get_xticklabels(), rotation=90)

    sc.pl.embedding(adata, basis = "X_emb_2", color='sample', ax = ax[0,3], frameon=False, legend_loc='right margin')
    ax[0,3].set_xticklabels(ax[0,3].get_xticklabels(), rotation=90)

    sc.pl.embedding(adata, basis = "X_emb_2", color = 'specie', ax=ax[1,0], frameon=False, legend_loc='right margin')
    ax[1,0].set_xticklabels(ax[1,0].get_xticklabels(), rotation=90)

    sc.pl.embedding(adata, basis = "X_emb_2", color = 'cell_type', ax=ax[1,1], frameon=False, legend_loc='right margin')
    ax[1,1].set_xticklabels(ax[1,1].get_xticklabels(), rotation=90)

    sc.pl.embedding(adata, basis = "X_emb_2", color = 'cell_annot_comb', ax=ax[1,2], frameon=False, legend_loc='right margin')
    ax[1,2].set_xticklabels(ax[1,2].get_xticklabels(), rotation=90)

    sc.pl.embedding(adata, basis = "X_emb_2", color = 'cell_annot', ax=ax[1,3], frameon=False, legend_loc='right margin')
    ax[1,3].set_xticklabels(ax[1,3].get_xticklabels(), rotation=90)

    #fig.subplots_adjust(wspace=0.2, hspace=0.2)
    plt.tight_layout()
    fig.savefig(dim_red_plot)

else:
    
    # Dimensionality reduction
    fig, ax = plt.subplots(figsize=(20,10), dpi=300, ncols=3, nrows=2)

    sc.pl.embedding(adata, basis = "X_emb_2", color = "cluster", ax = ax[0,0], frameon=False, legend_loc='on data')
    ax[0,0].set_xticklabels(ax[0,0].get_xticklabels(), rotation=90)

    sc.pl.embedding(adata, basis = "X_emb_2", color='study', ax = ax[0,1], frameon=False, legend_loc='right margin')
    ax[0,1].set_xticklabels(ax[0,1].get_xticklabels(), rotation=90)

    sc.pl.embedding(adata, basis = "X_emb_2", color='gender', ax = ax[0,2], frameon=False, legend_loc='right margin')
    ax[0,2].set_xticklabels(ax[0,2].get_xticklabels(), rotation=90)

    sc.pl.embedding(adata, basis = "X_emb_2", color='sample', ax = ax[1,0], frameon=False, legend_loc='right margin')
    ax[1,0].set_xticklabels(ax[1,0].get_xticklabels(), rotation=90)

    sc.pl.embedding(adata, basis = "X_emb_2", color = 'specie', ax=ax[1,1], frameon=False, legend_loc='right margin')
    ax[1,1].set_xticklabels(ax[1,1].get_xticklabels(), rotation=90)

    sc.pl.embedding(adata, basis = "X_emb_2", color = 'cell_type', ax=ax[1,2], frameon=False, legend_loc='right margin')
    ax[1,2].set_xticklabels(ax[1,2].get_xticklabels(), rotation=90)

    #fig.subplots_adjust(wspace=0.2, hspace=0.2)
    plt.tight_layout()
    fig.savefig(dim_red_plot)

# Dendrogram plot
fig, ax = plt.subplots(figsize=(9,5), dpi=300, ncols=1)
sc.pl.dendrogram(adata, groupby = 'cluster', ax=ax)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
fig.savefig(dendrogram_plot)

# Correlation matrix
sc.pl.correlation_matrix(
    adata, 
    groupby = "cluster"
)
fig = plt.gcf()
fig.savefig(corr_mat_plot)

