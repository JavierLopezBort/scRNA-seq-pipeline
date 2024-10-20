import scanpy as sc
#import scvelo as scv
import pandas as pd
import matplotlib.pyplot as plt
#import diffxpy.api as de
import numpy as np

#########################################################3
# LOAD RAW COUNT MATRIX
#########################################################3

adata = sc.read_h5ad(snakemake.input.adata)
rank_genes_table = pd.read_csv(snakemake.input.DGE_table)

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

# Params
cat_var_rank = snakemake.params.cat_var_rank

# Output plots
dim_red_rank = snakemake.output.dim_red_rank
volcano_plot = snakemake.output.volcano_plot
volcano_plot_2 = snakemake.output.volcano_plot_2
DGE_table_filt = snakemake.output.DGE_table_filt

#########################################################3
# RANK GENES
#########################################################3

stat_sig = "qval"
stat_sig_threshold = 0.05
mean_threshold = 1
stat = "log2fc"
stat_upper_threshold = 0.0001
stat_lower_threshold = -0.0001

rank_genes_table = rank_genes_table[rank_genes_table[stat_sig] <= stat_sig_threshold]
rank_genes_table = rank_genes_table[rank_genes_table['mean'] > mean_threshold]
rank_genes_table = rank_genes_table[(rank_genes_table[stat] >= stat_upper_threshold) | (rank_genes_table[stat] <= stat_lower_threshold)]
rank_genes_table = rank_genes_table.sort_values(by=stat, ascending=False)

rank_genes_table.to_csv(DGE_table_filt, index=False)

top_genes = rank_genes_table.head(3)
bottom_genes = rank_genes_table.tail(2)
rank_genes_table_filt = pd.concat([top_genes, bottom_genes], axis=0, ignore_index=True)
rank_genes_table_filt['sample'] = np.where(rank_genes_table_filt[stat] > 0, 'HEK', 'PGCLCs')

rank_genes_list = []

for sample, gene, score in zip(rank_genes_table_filt['sample'].tolist(), rank_genes_table_filt['gene'].tolist(), rank_genes_table_filt[stat].tolist()):
    rank_genes_list.append([gene, sample, score])

logger.info(rank_genes_list)

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
                sc.pl.embedding(adata, basis = "X_emb_2", color = rank_genes_list[i][0], ax = ax[row,col], frameon=False, title = f"{rank_genes_list[i][0]}, {rank_genes_list[i][1]}, {round(rank_genes_list[i][2], 2)}")
                ax[row,col].set_xticklabels(ax[row,col].get_xticklabels(), rotation=90)
                i += 1
        else:
            fig.delaxes(ax[row,col])

plt.tight_layout()
fig.savefig(dim_red_rank)
#fig.delaxes(ax[1,2])
#fig.subplots_adjust(wspace=0.2, hspace=0.2)

# Volcano plot
plt.figure(figsize=(8, 5))
plt.scatter(rank_genes_table[stat], -np.log10(rank_genes_table[stat_sig]), c='blue', label='Genes', s=1)
plt.title(f'Volcano plot of {stat} vs {stat_sig}')
plt.xlabel(stat)
plt.ylabel(f"-log10({stat_sig})")
plt.legend()
plt.grid(False)
plt.show()
plt.savefig(volcano_plot)

# Volcano plot 2
plt.figure(figsize=(8, 5))
plt.scatter(rank_genes_table[stat], rank_genes_table['mean'], c='blue', label='Genes', s=1)
plt.title(f'Scatter plot of {stat} vs mean')
plt.xlabel(stat)
plt.ylabel('mean')
plt.legend()
plt.grid(False)
plt.show()
plt.savefig(volcano_plot_2)


