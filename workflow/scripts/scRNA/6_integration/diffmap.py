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

# common
n_neighbours = snakemake.params.n_neighbours
n_pcs = snakemake.params.n_pcs
resolution = snakemake.params.resolution
cluster = snakemake.params.cluster
layer_pca = snakemake.params.layer_pca

#########################################################3
# 0 ADATA UNINTEGRATED
#########################################################3

adata_pre = adata

#########################################################3
# 1 ST DIMENSIONALITY REDUCTION
#########################################################3

# PCA
sc.tl.pca(adata, svd_solver='arpack', layer = layer_pca, key_added = "X_emb")

#########################################################3
# NEIGHBOURS
#########################################################3

sc.pp.neighbors(adata,
    n_neighbors = n_neighbours,
    n_pcs = n_pcs,
    use_rep = "X_emb")

sc.pp.neighbors(adata_pre,
    n_neighbors = n_neighbours)

#########################################################3
# 1 ST DIMENSIONALITY REDUCTION
#########################################################3

sc.tl.diffmap(adata) # Denoise the graph, it needs previous neighbours.

adata.obsm['X_emb'] = adata.obsm['X_diffmap']

#########################################################3
# NEIGHBOURS
#########################################################3

sc.pp.neighbors(adata, 
    n_neighbors = n_neighbours, 
    #n_pcs = n_pcs,
    use_rep = "X_emb"
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
    n_pcs = n_pcs)

#########################################################3
# 2 ST DIMENSIONALITY REDUCTION
#########################################################3

#sc.tl.umap(adata, n_components=2)
sc.tl.draw_graph(adata)

adata.obsm['X_emb_2'] = adata.obsm['X_draw_graph_fa']

#########################################################3
# WRITE COUNT MATRIX
#########################################################3

adata.write_h5ad(snakemake.output.adata)
