import scanpy as sc
import pandas as pd
import scvelo as scv

#########################################################3
# LOAD RAW COUNT MATRIX
#########################################################3

adata = sc.read_h5ad(snakemake.input.adata)

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

# Rank genes
cat_var_rank = snakemake.params.cat_var_rank
#n_rank_genes = snakemake.params.n_rank_genes

#########################################################3
# RANK GENES
#########################################################3

scv.tl.rank_velocity_genes(
    adata, 
    groupby = cat_var_rank,
    #n_genes = n_rank_genes
)

#########################################################3
# STORE RANK GENES TABLE
#########################################################3

rank_genes_names_df = pd.DataFrame(adata.uns['rank_velocity_genes']['names'])
rank_genes_scores_df = pd.DataFrame(adata.uns['rank_velocity_genes']['scores'])
rank_genes_df = pd.concat([rank_genes_names_df, rank_genes_scores_df], axis=1)

rank_genes_df.to_csv(snakemake.output.rank_genes_table, index=False)

#########################################################3
# WRITE COUNT MATRIX
#########################################################3

adata.write_h5ad(snakemake.output.adata)
