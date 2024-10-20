import scanpy as sc
import matplotlib.pyplot as plt
from logging_func import handle_exception, setup_logger
from scvi.model import SCVI 

#########################################################3
# LOAD RAW COUNT MATRIX
#########################################################3

adata = sc.read_h5ad(snakemake.input.adata)
model = SCVI.load(snakemake.input.model_dir, adata)
model.to_device('cuda')
logger.info(model)

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

# Params
cat_var_rank = snakemake.params.cat_var_rank

#########################################################3
# CREATE DE MODEL
#########################################################3

import torch
logger.info(f"torch.cuda is available: {torch.cuda.is_available()}")

rank_genes_df = model.differential_expression(
    groupby = cat_var_rank,
)

logger.info(rank_genes_df.head())

rank_genes_df.to_csv(snakemake.output.rank_genes_table)


