import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
from logging_func import handle_exception, setup_logger

#########################################################3
# LOAD RAW COUNT MATRIX
#########################################################3

rank_genes_table = pd.read_csv(snakemake.input.rank_genes_table)
adata = sc.read_h5ad(snakemake.input.adata)

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

logger = setup_logger(snakemake.log[0])

# Params
cat_var_rank = snakemake.params.cat_var_rank
#n_rank_genes = snakemake.params.n_rank_genes

# Output plots
dim_red_rank = snakemake.output.dim_red_rank
#rank_plot = snakemake.output.rank_plot
#rank_dot_plot = snakemake.output.rank_dot_plot

#########################################################3
# RANK GENES
#########################################################3

top_gene = 1

#first_rows = rank_genes_table.groupby('group').first().reset_index()
top_gene_row = rank_genes_table.groupby('group', as_index=False).nth(top_gene - 1)
logger.info(top_gene_row)
rank_genes_list = list(zip(top_gene_row['names'], top_gene_row['group'], top_gene_row['scores']))
rank_genes = top_gene_row['names']

logger.info(rank_genes_list)
logger.info(f"rank genes: {rank_genes}")
logger.info(f"Number of rank genes / samples: {len(rank_genes)}")

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
                sc.pl.embedding(adata, basis = "X_emb_2", color = cat_var_rank, ax = ax[row,col], frameon=False, legend_loc = "right margin")
                ax[row,col].set_xticklabels(ax[row,col].get_xticklabels(), rotation=90)
            else:
                sc.pl.embedding(adata, basis = "X_emb_2", color = rank_genes_list[i][0], ax = ax[row,col], frameon=False, title = f"{rank_genes_list[i][0]}, {rank_genes_list[i][1]}, {rank_genes_list[i][2]}")
                ax[row,col].set_xticklabels(ax[row,col].get_xticklabels(), rotation=90)
                i += 1
        else:
            fig.delaxes(ax[row,col])
        

plt.tight_layout()
fig.savefig(dim_red_rank)
#fig.delaxes(ax[1,2])
#fig.subplots_adjust(wspace=0.2, hspace=0.2)

"""
fig, ax = plt.subplots(figsize=(9,5), dpi=300, ncols=1)
sc.pl.rank_genes_groups(adata, n_genes = n_rank_genes, sharey = False, ax = ax)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
fig.savefig(rank_plot)

fig, ax = plt.subplots(figsize=(9,5), dpi=300, ncols=1)
sc.pl.rank_genes_groups_dotplot(adata, n_genes = n_rank_genes, 
values_to_plot='logfoldchanges', min_logfoldchange=4, vmax=7, vmin=-7, cmap='bwr', ax = ax)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
fig.savefig(rank_dot_plot)
"""
