import pandas as pd
from logging_func import handle_exception, setup_logger
import scanpy as sc
import numpy as np

#########################################################3
# LOAD RAW COUNT MATRIX
#########################################################3

adata = sc.read_h5ad(snakemake.input.adata)

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

logger = setup_logger(snakemake.log[0])
threshold = snakemake.params.threshold

#########################################################3
# CELL ANNOTATION
#########################################################3

boolean = adata.obsm['X_emb_2'][:, 0] < threshold
category = np.where(boolean, 'HEK', 'iPSC')
adata.obs["cell_annot"] = category
adata.obs["cell_annot_comb"] = category

mask_HEK = (adata.obs["cell_annot"] == "HEK").values
mask_iPSC = (adata.obs["cell_annot"] == "iPSC").values

mask_HTO_C0251 = (adata.obs["sample"] == "HTO_C0251").values
mask_HTO_C0252 = (adata.obs["sample"] == "HTO_C0252").values
mask_HTO_C0253 = (adata.obs["sample"] == "HTO_C0253").values
mask_HTO_C0254 = (adata.obs["sample"] == "HTO_C0254").values
mask_HTO_C0255 = (adata.obs["sample"] == "HTO_C0255").values
mask_HTO_C0256 = (adata.obs["sample"] == "HTO_C0256").values

mask_missclass_HEK = (mask_HEK & (mask_HTO_C0254 | mask_HTO_C0255 | mask_HTO_C0256))
mask_missclass_iPSC = (mask_iPSC & (mask_HTO_C0251 | mask_HTO_C0252 | mask_HTO_C0253))

n_HEK = sum(mask_HEK)
n_iPSC = sum(mask_iPSC)
n_HEK_missclass = sum(mask_missclass_HEK)
n_iPSC_missclass = sum(mask_missclass_iPSC)

miss_class_err_HEK = (n_HEK_missclass / n_HEK) * 100
miss_class_err_iPSC = (n_iPSC_missclass / n_iPSC) * 100
miss_class_cell = ((n_iPSC_missclass + n_HEK_missclass) / (n_iPSC + n_HEK)) * 100

logger.info(F"Percentage of iPSC classified as HEK: {miss_class_err_HEK}")
logger.info(F"Percentage of HEK classified as iPSC: {miss_class_err_iPSC}")
logger.info(F"Percentage of HEK and iPSC missclassified: {miss_class_cell}")

#########################################################3
# WRITE COUNT MATRIX
#########################################################3

adata.write_h5ad(snakemake.output.adata)