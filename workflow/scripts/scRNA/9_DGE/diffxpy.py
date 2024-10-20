import anndata
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.stats
import diffxpy.api as de
import pickle
from logging_func import handle_exception, setup_logger

#########################################################3
# LOAD RAW COUNT MATRIX
#########################################################3

adata = sc.read_h5ad(snakemake.input.adata)

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

logger = setup_logger(snakemake.log[0])

# Rank genes
cat_var_rank = snakemake.params.cat_var_rank
#n_rank_genes = snakemake.params.n_rank_genes
cont_covariates_diffxpy = snakemake.params.cont_covariates_diffxpy
dge_method = snakemake.params.dge_method
cell_inducted = snakemake.params.cell_inducted

#########################################################3
# REANNOTATE
#########################################################3

if dge_method == "standard":
    pass

elif dge_method == "meioticinduction":
    
    adata.obs['cell_annot_comb'] = adata.obs['cell_annot_comb'].astype(str)
    adata.obs['cell_annot_comb'] = adata.obs['cell_annot_comb'].replace('OCT', 'SCT')
    adata.obs['cell_annot_comb'] = adata.obs['cell_annot_comb'].replace('PGCLCs', 'PGCs')
    adata.obs['cell_annot_comb'] = pd.Categorical(adata.obs['cell_annot_comb'])

    mask_1 = adata.obs["cell_annot_comb"] == "PGCs"
    mask_2 = adata.obs["cell_annot_comb"] == cell_inducted
    mask_combined = mask_1 | mask_2
    adata = adata[mask_combined, :]

#########################################################3
# WALD TEST
#########################################################3

info = de.utils.preview_coef_names(
    sample_description = adata.obs,
    formula=f"~ 1 + {cat_var_rank}"
)

logger.info(info)

adata.X = adata.layers["counts"]

logger.info(adata.obs["cell_annot"])

test = de.test.wald(
    data = adata,
    #data = adata.layers["counts"]
    #sample_description = adata.obs,
    formula_loc = f"~ 1 + {cat_var_rank} + {cont_covariates_diffxpy[0]}",
    factor_loc_totest = cat_var_rank,
    as_numeric = cont_covariates_diffxpy,
    #coef_to_test = 'PGCLCs_DM+'
)

test_df = test.summary()

#########################################################3
# STORE RANK GENES TABLE
#########################################################3

test_df.to_csv(snakemake.output.DGE_table, index=False)

#########################################################3
# WRITE COUNT MATRIX
#########################################################3

adata.write_h5ad(snakemake.output.adata)



