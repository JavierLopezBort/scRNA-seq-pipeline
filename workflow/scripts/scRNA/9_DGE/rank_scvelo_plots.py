import scanpy as sc
import scvelo as scv
import pandas as pd
import matplotlib.pyplot as plt

#########################################################3
# LOAD RAW COUNT MATRIX
#########################################################3

adata = sc.read_h5ad(snakemake.input.adata)
rank_genes_table = pd.read_csv(snakemake.input.rank_genes_table)

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

# Params
cat_var_rank = snakemake.params.cat_var_rank

# Output plots
dim_red_rank = snakemake.output.dim_red_rank

#########################################################3
# RANK GENES
#########################################################3

top_gene = 1

top_gene = top_gene
rank_samples = list(adata.obs[cat_var_rank].unique()) 
n_categories = len(rank_samples)

rank_genes = rank_genes_table.iloc[top_gene - 1, 0:n_categories]
rank_scores = rank_genes_table.iloc[top_gene - 1, n_categories:]

logger.info(f"rank genes: {rank_genes}")
logger.info(f"Number of rank genes / samples: {n_categories}")

rank_genes_list = []

for sample, gene, score in zip(rank_samples, rank_genes, rank_scores):
    rank_genes_list.append([gene, sample, score])

#########################################################3
# RANK GENES PLOTS
#########################################################3

ncols = 3
nrows = 2

# Dimensionality reduction
fig, ax = plt.subplots(figsize=(20,10), dpi=300, ncols = ncols, nrows = nrows)

i = 0
for row in range(nrows):
    for col in range(ncols):
        if i < len(rank_genes_list):
            if row == 0 and col == 0:
                scv.pl.velocity_embedding_stream(adata, basis = "X_emb_2", color = cat_var_rank, ax = ax[row,col], frameon=False, legend_loc = "right margin")
                ax[row,col].set_xticklabels(ax[row,col].get_xticklabels(), rotation=90)
            else:
                scv.pl.velocity_embedding_stream(adata, basis = "X_emb_2", color = rank_genes_list[i][0], ax = ax[row,col], frameon=False, title = f"{rank_genes_list[i][0]}, {rank_genes_list[i][1]}, {rank_genes_list[i][2]}")
                ax[row,col].set_xticklabels(ax[row,col].get_xticklabels(), rotation=90)
                i += 1
        else:
            fig.delaxes(ax[row,col])

plt.tight_layout()
fig.savefig(dim_red_rank)
#fig.delaxes(ax[1,2])
#fig.subplots_adjust(wspace=0.2, hspace=0.2)





