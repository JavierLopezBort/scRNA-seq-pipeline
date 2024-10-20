import scanpy as sc
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

layer_pca = snakemake.params.layer_pca
n_pcs = snakemake.params.n_pcs
cat_covariate = snakemake.params.cat_covariate

#########################################################3
# 0 ADATA UNINTEGRATED
#########################################################3

adata_pre = adata

#########################################################3
# 1 ST DIMENSIONALITY REDUCTION
#########################################################3

# PCA
sc.tl.pca(adata, svd_solver='arpack', layer = layer_pca)

adata.obsm['X_emb'] = adata.obsm['X_pca']

#########################################################3
# NEIGHBOURS
#########################################################3

sc.external.pp.bbknn(adata, 
    neighbors_within_batch = n_neighbours, 
    n_pcs = n_pcs, 
    use_rep = "X_emb",
    batch_key = cat_covariate)

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
    use_rep = "X_emb", 
    n_pcs = n_pcs)

#########################################################3
# 2 ST DIMENSIONALITY REDUCTION
#########################################################3

sc.tl.umap(adata, n_components=2)

adata.obsm['X_emb_2'] = adata.obsm['X_umap']

#sc.tl.draw_graph(adata)
#sc.tl.diffmap(adata)

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
