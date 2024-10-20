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

cluster_root = snakemake.params.cluster_root

#########################################################3
# VELOCITY PAGA
#########################################################3

sc.tl.paga(adata, groups = "cluster")
sc.pl.paga(adata, plot=False) 

#########################################################3
# 2 ST DIMENSIONALITY REDUCTION BASED ON PAGA
#########################################################3

sc.tl.draw_graph(adata, init_pos="paga")
#sc.tl.umap(adata, init_pos="paga")

adata.obsm['X_emb_2'] = adata.obsm['X_draw_graph_fa']

#########################################################3
# CLUSTER ROOT
#########################################################3

adata.uns["iroot"] = np.flatnonzero(adata.obs["cell_annot_comb"] == str(cluster_root))[0]

#########################################################3
# PSEUDOTIME
#########################################################3

sc.tl.dpt(adata)

#########################################################3
# WRITE COUNT MATRIX
#########################################################3

adata.write_h5ad(snakemake.output.adata)
