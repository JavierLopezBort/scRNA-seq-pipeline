import scanpy as sc
import pandas as pd

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

stat_method = 't-test'
groups = 'all'        # 'all' ["0"]
reference = 'rest'    # 'rest' "1" 

sc.tl.rank_genes_groups(
    adata, 
    groupby = cat_var_rank, 
    method = stat_method, 
    groups = groups,
    reference = reference,
    key_added = "extract_later"
    #n_rank_genes = n_rank_genes
)

#########################################################3
# STORE RANK GENES TABLE
#########################################################3

rank_genes_df = sc.get.rank_genes_groups_df(adata, None, key = "extract_later")

rank_genes_df.to_csv(snakemake.output.rank_genes_table, index=False)

#########################################################3
# WRITE COUNT MATRIX
#########################################################3

adata.write_h5ad(snakemake.output.adata)
