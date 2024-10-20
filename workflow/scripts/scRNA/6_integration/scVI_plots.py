import scanpy as sc
import matplotlib.pyplot as plt
import louvain
import torch

#########################################################3
# LOAD RAW COUNT MATRIX
#########################################################3

adata = sc.read_h5ad(snakemake.input.adata)

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

# Output plots
dim_red_plot = snakemake.output.dim_red_plot
dendrogram_plot = snakemake.output.dendrogram_plot
corr_mat_plot = snakemake.output.corr_mat_plot
cluster_counts_plot = snakemake.output.cluster_counts_plot

# SCVI
model = torch.load(snakemake.input.model)
elbo_plot = snakemake.output.elbo_plot

train_elbo = model["attr_dict"]["history_"]["elbo_train"][1:]
test_elbo = model["attr_dict"]["history_"]["elbo_validation"][1:]

#########################################################3
# PLOTS
#########################################################3

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

#SCVI plot
fig, ax = plt.subplots()
ax.plot(train_elbo, label="train")
ax.plot(test_elbo, label="test")
ax.spines[['right', 'top']].set_visible(False)
ax.set_xlabel('Epoch')
ax.legend()
fig.savefig(elbo_plot)

# Cluster counts
nrows = 2
ncols = 9

fig, ax = plt.subplots(figsize=(20,10), dpi=300, ncols=ncols, nrows=nrows)
cluster_list = list(adata.obs["cluster"].unique())

i = 0
for row in range(nrows):
    ax[row,0].set_ylabel('Perc cells')
    for col in range(ncols):
        if i < len(cluster_list):
            if row == 0 and col == 0:
                sc.pl.embedding(adata, basis = "X_emb_2", color = "cluster", ax = ax[row,col], frameon=False, legend_loc='on data')
                ax[row,col].set_xticklabels(ax[row,col].get_xticklabels(), rotation=90)
            else:
                cluster = cluster_list[i]
                adata_subset = adata[adata.obs["cluster"] == cluster, :]
                var = adata_subset.obs["sample"]
                var_counts = var.value_counts()
                bars = var_counts.plot(kind='bar', ax=ax[row,col])
                ax[row,col].set_title(cluster)
                ax[row,col].tick_params(axis='x', rotation=45)

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
fig.savefig(cluster_counts_plot)



