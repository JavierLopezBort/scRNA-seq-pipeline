import torch
import sys, os
import warnings

import scanpy as sc
from dataclasses import dataclass
from pytorch_lightning import Trainer
import torch
import scvi
from scvi.model import SCVI 
from scvi.model.utils import mde
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numba.core.errors import NumbaDeprecationWarning
from logging_func import handle_exception, setup_logger
import louvain
import numpy as np

#########################################################3
# LOAD RAW COUNT MATRIX
#########################################################3

adata = sc.read_h5ad(snakemake.input.adata)

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

logger = setup_logger(snakemake.log[0])

n_neighbours = snakemake.params.n_neighbours
resolution = snakemake.params.resolution
cluster = snakemake.params.cluster

# scvi
n_latent = snakemake.params.n_latent
batch_size = snakemake.params.batch_size
cat_covariates = snakemake.params.cat_covariates
cont_covariates = snakemake.params.cont_covariates
n_epochs = snakemake.params.n_epochs

#########################################################3
# 0 ADATA UNINTEGRATED
#########################################################3

adata_pre = adata

#########################################################3
# 1 ST DIMENSIONALITY REDUCTION
#########################################################3

# CHECK IF CUDA AVAILABLE
if torch.cuda.is_available():
    logger.info(f"A GPU is available.")
    logger.info(f"GPU Name: {torch.cuda.get_device_name(0)}")
    logger.info(f"Total number of GPUs: {torch.cuda.device_count()}")
else:
    logger.info("No GPU is available.")

# SCVI
warnings.filterwarnings(action="ignore", category=NumbaDeprecationWarning)
warnings.filterwarnings(
    action="ignore", module="scanpy", message="No data for colormapping"
)

SEED = 0
scvi.settings.seed = SEED

logger.info(f'Setting seed to {SEED}')

torch.set_float32_matmul_precision('medium')
Trainer(accelerator='cuda')

# Filter on low batch count samples
batch_counts = adata.obs.value_counts("sample")

# continuous and categorical must be variables highly dependent on non biological reasons, batch effects.
# The study where sample comes from is indeed the batch, but there can be other variables

SCVI.setup_anndata(
    adata, 
    layer='counts', 
    categorical_covariate_keys = cat_covariates,
    continuous_covariate_keys= cont_covariates
)

model = SCVI(adata, n_latent=n_latent)
model.train(
    max_epochs = n_epochs,
    check_val_every_n_epoch=50, 
    #check_val_every_n_epoch=None, 
    early_stopping=True, 
    batch_size=batch_size,
    #train_size=0.8,
    #validation_size=0.2,
    progress_bar_refresh_rate=0,
)
latent = model.get_latent_representation()
# Save only the latent

adata.obsm["X_emb"] = latent
adata.layers["scVI"] = model.get_normalized_expression()
adata.X = adata.layers["scVI"]

# Save model
model.save(snakemake.output.model_dir, overwrite=True, save_ann_data=False)

#########################################################3
# NEIGHBOURS
#########################################################3

sc.pp.neighbors(adata, 
    n_neighbors = n_neighbours, 
    use_rep = "X_emb"
    )

sc.pp.neighbors(adata_pre,
    n_neighbors = n_neighbours)

#########################################################3
# CLUSTERING
#########################################################3

if cluster == "leiden":
    
    sc.tl.leiden(adata, resolution = resolution, key_added = "cluster")
    sc.tl.leiden(adata_pre, resolution = resolution, key_added = "cluster")

elif cluster == "louvain":
    sc.tl.louvain(adata, resolution = resolution, key_added = "cluster")
    sc.tl.louvain(adata_pre, resolution = resolution, key_added = "cluster")
    

sc.tl.dendrogram(adata, 
    groupby = "cluster",
    use_rep = "X_emb")

#########################################################3
# 2 ST DIMENSIONALITY REDUCTION
#########################################################3

#sc.tl.umap(adata, n_components=2)
#sc.tl.diffmap(adata)
#sc.tl.draw_graph(adata)
adata.obsm["X_emb_2"] = mde(adata.obsm["X_emb"])

#########################################################3
# VELOCITY
#########################################################3
"""
if "cell_annot_comb" in adata.obs.columns:
    sc.tl.diffmap(adata_pre)
    sc.tl.diffmap(adata)

    adata_pre.uns["iroot"] = np.flatnonzero(adata.obs["cell_annot_comb"] == "SPG")[0]
    adata.uns["iroot"] = np.flatnonzero(adata.obs["cell_annot_comb"] == "SPG")[0]

    sc.tl.dpt(adata_pre)
    sc.tl.dpt(adata)
"""
#########################################################3
# WRITE COUNT MATRIX
#########################################################3

adata.write_h5ad(snakemake.output.adata_int)
adata_pre.write_h5ad(snakemake.output.adata)
