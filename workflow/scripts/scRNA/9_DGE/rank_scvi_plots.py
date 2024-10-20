import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from logging_func import handle_exception, setup_logger

#########################################################3
# LOAD RAW COUNT MATRIX
#########################################################3

adata = sc.read_h5ad(snakemake.input.adata)
rank_genes_table = pd.read_csv(snakemake.input.rank_genes_table)

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

logger = setup_logger(snakemake.log[0])

# Params
cat_var_rank = snakemake.params.cat_var_rank

# Output plots
dim_red_rank = snakemake.output.dim_red_rank

#########################################################3
# RANK GENES
#########################################################3

lfc_mean = 0
bayes_factor = 3
non_zeros_proportion1 = 0.001 # 0.1
top_gene = 1

cat_var_rank_list = list(adata.obs[cat_var_rank].unique())
rank_genes_list = list()
rank_genes = list()

for i, c in enumerate(cat_var_rank_list):
    
    cid = f"{c} vs Rest"
    cell_type_df = rank_genes_table.loc[rank_genes_table.comparison == cid]
    
    cell_type_df = cell_type_df[cell_type_df.lfc_mean > lfc_mean]
    cell_type_df = cell_type_df[cell_type_df["bayes_factor"] > bayes_factor]
    cell_type_df = cell_type_df[cell_type_df["non_zeros_proportion1"] > non_zeros_proportion1]
    
    top_gene_row = cell_type_df.iloc[top_gene - 1]
    gene = top_gene_row.loc["gene_ids"]
    sample = c
    score = top_gene_row.loc["proba_de"]
    
    rank_genes_list.append([gene, sample, score])
    rank_genes.append(gene)

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