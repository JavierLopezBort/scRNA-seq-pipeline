#########################################################3
# IMPORT PACKAGES
#########################################################3

import sys
import os
from pathlib import Path
from os import PathLike
from scipy.sparse import csr_matrix

import numpy as np
import pandas as pd
from scipy.stats import mode
import scanpy as sc
from typing import Union, Tuple, Optional, Dict, Any
import sklearn
import warnings
import h5py

sys.path.insert(0, "../")
import scgpt as scg
import faiss
#from build_atlas_index_faiss import load_index, vote, compute_category_proportion

from tqdm import tqdm

#########################################################3
# IMPORT FROM SNAKEMAKE
#########################################################3

adata = sc.read_h5ad(snakemake.input.adata)
model_dir = Path(snakemake.input.model_dir)

n_neighbours = snakemake.params.n_neighbours
resolution = snakemake.params.resolution
cluster = snakemake.params.cluster

#########################################################3
# 0 ADATA UNINTEGRATED
#########################################################3

adata_pre = adata

#########################################################3
# 1 ST DIMENSIONALITY REDUCTION
#########################################################3

gene_col = 'index'

adata_test = adata
adata_test.X = (adata_test.X).todense()

test_embed_adata = scg.tasks.embed_data(
    adata_test,
    model_dir,
    gene_col=gene_col,
    #obs_to_save=cell_type_key,  # optional arg, only for saving metainfo
    batch_size=64,
    return_new_adata=True,
)

test_embed = test_embed_adata.X

adata.obsm['X_emb'] = test_embed

#########################################################3
# NEIGHBOURS
#########################################################3

sc.pp.neighbors(adata_pre,
    n_neighbors = n_neighbours,
    #n_pcs = 3
    )

sc.pp.neighbors(adata, 
    n_neighbors = n_neighbours, 
    use_rep = "X_emb",
    #n_pcs = 3
    )

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
    use_rep = "X_emb",
    n_pcs = 3)

#########################################################3
# 2 ST DIMENSIONALITY REDUCTION
#########################################################3

sc.tl.umap(adata, n_components=2)

#sc.tl.diffmap(adata)
#sc.tl.draw_graph(adata)

adata.obsm['X_emb_2'] = adata.obsm['X_umap']

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