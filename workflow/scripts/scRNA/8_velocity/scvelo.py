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

import scipy

#########################################################3
# LOAD RAW COUNT MATRIX
#########################################################3

adata = sc.read_h5ad(snakemake.input.adata)

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

logger = setup_logger(snakemake.log[0])

logger.info('Starting scvelo.py')
logger.info('Loading velocity data')

#########################################################3
# MATRIX TRANSFORMATION (RETRANSFORMATION)
#########################################################3

print(scipy.__version__)

logger.info(adata)
logger.info(type(adata.X))
logger.info(type(adata.layers["spliced"]))
logger.info(type(adata.layers["unspliced"]))
logger.info(type(adata.obsp["connectivities"]))
logger.info(type(adata.obsp["distances"]))

#logger.info(adata.layers["spliced"].A)

#adata.X = scipy.sparse.csr_matrix(adata.X)
adata.layers["spliced"] = adata.layers["spliced"].toarray()
adata.layers["unspliced"] = adata.layers["unspliced"].toarray()
adata.obsp["connectivities"] = adata.obsp["connectivities"].toarray()
adata.obsp["distances"] = adata.obsp["distances"].toarray()

logger.info(adata)
logger.info(type(adata.X))
logger.info(type(adata.layers["spliced"]))
logger.info(type(adata.layers["unspliced"]))
logger.info(type(adata.obsp["connectivities"]))
logger.info(type(adata.obsp["distances"]))

#########################################################3
# VELOCITY
#########################################################3

adata.obs['clusters'] = adata.obs["cluster"]

scv.tl.recover_dynamics(adata, n_jobs=15)

logger.info(adata)
logger.info(type(adata.X))
logger.info(type(adata.layers["spliced"]))
logger.info(type(adata.layers["unspliced"]))
logger.info(type(adata.obsp["connectivities"]))
logger.info(type(adata.obsp["distances"]))

scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata, n_jobs=20)

#########################################################3
# WRITE COUNT MATRIX
#########################################################3

adata.write_h5ad(snakemake.output.adata)