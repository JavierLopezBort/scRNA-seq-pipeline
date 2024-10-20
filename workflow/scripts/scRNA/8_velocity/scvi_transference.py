import numpy as np
import pandas as pd

import scvelo as scv
import scanpy as sc

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from logging_func import handle_exception, setup_logger
import louvain

#########################################################3
# LOAD RAW COUNT MATRIX
#########################################################3

adata = sc.read_h5ad(snakemake.input.adata)
fdata = sc.read_h5ad(snakemake.input.fdata)

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

scvi_transference = snakemake.params.scvi_transference

#########################################################3
# SCVI TRANSFERENCE
#########################################################3

if scvi_transference:
    common_idx = np.intersect1d(fdata.obs_names, adata.obs_names)
    logger.info(f'Found {len(common_idx)} common cells between analysis')
    fdata = fdata[common_idx, :]
    adata = adata[common_idx, :]
    # Move over the scVI latent space
    adata.obsm['X_scVI'] = fdata.obsm['X_scVI']
    adata.obsm['X_mde'] = fdata.obsm['X_mde']

#########################################################3
# WRITE COUNT MATRIX
#########################################################3

adata.write_h5ad(snakemake.output.adata)