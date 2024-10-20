import scanpy as sc
import scvelo as scv
from logging_func import handle_exception, setup_logger # the function.py file has to be in the same folder as this python file
#from ivynatal_ref import meiosis_genes
import anndata as an
import numpy as np

#########################################################3
# LOAD RAW COUNT MATRIX
#########################################################3

adata = sc.read_h5ad(snakemake.input.adata)

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

logger = setup_logger(snakemake.log[0])

n_highly_var = snakemake.params.n_highly_var

hvg = snakemake.params.hvg
scaling = snakemake.params.scaling

# Gene list
gene_list = snakemake.params.gene_list

#########################################################3
# INITIAL COUNT
#########################################################3

# GENERAL
logger.info(f'Number of genes: {adata.shape[1]}')

# CHECK GENES
for gene in gene_list:
    logger.info(f'{gene} still in data: {gene in adata.var_names}')

#########################################################3
# NORMALIZATION
#########################################################3

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Save counts data
adata.layers['normalized'] = adata.X.copy()

#########################################################3
# FILTER HIGHLY VARIABLE GENES
#########################################################3

if hvg:

    # in highly variables, you can either specify the number of genes you want to keep, or decide several
    # thresholds (min_mean, max_mean, min_disp) and the number depends on these. However, only one option
    # in scvelo is exactly the same.

    # KEEP GENES BEFORE
    mask = adata.var_names.isin(gene_list)
    keep_genes_data = adata[:, mask]
    adata = adata[:, ~mask]

    # Filter highly variable genes
    #sc.pp.highly_variable_genes(adata, n_top_genes = n_highly_var, min_mean=0.0125, max_mean=3, min_disp=0.5)
    # This function expects log data and no extracts the highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes = n_highly_var, subset = True)

    # KEEP GENES AFTER
    adata_obs = adata.obs
    adata = an.concat([adata, keep_genes_data], axis = 1)
    adata.obs = adata_obs

#########################################################3
# SCALE
#########################################################3

if scaling:

    sc.pp.regress_out(adata, ['n_counts', 'pct_counts_mito'])
    sc.pp.scale(adata, max_value=10)
    # Store
    adata.layers["scaled"] = adata.X.copy()

    # Assign again to normalized, which is the default one for most of the methods (except of PCA, which is scaled)
    adata.X = adata.layers["normalized"]

#########################################################3
# FINAL COUNT
#########################################################3

# GENERAL
logger.info(f'Number of genes: {adata.shape[1]}')

# CHECK GENES
for gene in gene_list:
    logger.info(f'{gene} still in data: {gene in adata.var_names}')

#########################################################3
# WRITE COUNT MATRIX
#########################################################3

adata.write_h5ad(snakemake.output.adata)

