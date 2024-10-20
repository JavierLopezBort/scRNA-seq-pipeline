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

# mt and rb genes
mt_hard_filter = snakemake.params.mt_hard_filter
rb_hard_filter = snakemake.params.rb_hard_filter

# Filter cells
min_genes = snakemake.params.min_genes
min_counts_cells = snakemake.params.min_counts_cells

# Filter genes 
# Based on spliced. Filtered on spliced and unspliced
min_cells = snakemake.params.min_cells
min_counts_genes = snakemake.params.min_counts_genes
# Based on uspliced. Filtered on spliced and unspliced
min_cells_u = snakemake.params.min_cells_u
min_counts_genes_u = snakemake.params.min_counts_genes_u
# Based on spliced and unspliced. Filtered on spliced and unspliced
min_shared_cells = snakemake.params.min_shared_cells
min_shared_counts_genes = snakemake.params.min_shared_counts_genes
# Threshold doublets
threshold_doublets = snakemake.params.threshold_doublets

# Gene list
gene_list = snakemake.params.gene_list

#########################################################3
# INITIAL COUNT
#########################################################3

logger.info(f"\nBEFORE FILTERING\n")

# GENERAL
logger.info(f'Number of cells: {adata.shape[0]}')
logger.info(f'Number of genes: {adata.shape[1]}')
logger.info(f'Total spliced counts: {adata.X.sum()}')
logger.info(f'Total unspliced counts: {adata.layers["unspliced"].sum()}')
logger.info(f'Total counts: {adata.X.sum() + adata.layers["unspliced"].sum()}')
logger.info(f"Sample counts: {adata.obs['sample'].value_counts()}")

# SPLICED
# Counts cells
logger.info(f'Average number of counts per cell (spliced): {adata.obs["total_counts"].mean()}')
# Number of genes per cell
logger.info(f'Average number of genes per cell (spliced): {adata.obs["n_genes_by_counts"].mean()}')
# Counts genes
logger.info(f'Average number of counts per gene (spliced): {adata.var["total_counts"].mean()}')
# Number of cells per gene
logger.info(f'Average number of cells per gene (spliced): {adata.var["n_cells_by_counts"].mean()}')

# UNSPLICED
# Counts cells
logger.info(f'Average number of counts per cell (unspliced): {adata.obs["total_counts_u"].mean()}')
# Number of genes per cell
logger.info(f'Average number of genes per cell (unspliced): {adata.obs["n_genes_by_counts_u"].mean()}')
# Counts genes
logger.info(f'Average number of counts per gene (unspliced): {adata.var["total_counts_u"].mean()}')
# Number of cells per gene
logger.info(f'Average number of cells per gene (unspliced): {adata.var["n_cells_by_counts_u"].mean()}')

# CHECK GENES
for gene in gene_list:
    logger.info(f'{gene} still in data: {gene in adata.var_names}')
    
#########################################################3
# FILTER CELLS BASED ON MITO AND RIBO
#########################################################3

# REMEMBER
# All filterings are applied to all layers (.X and layers) even though the filtered method is based only on .X or a specific layer

adata = adata[adata.obs['pct_counts_mito'] < mt_hard_filter, :]
logger.info(f'Cells left after filtering based on pct_counts_mito: {adata.shape[0]}')

adata = adata[adata.obs['pct_counts_ribo'] < rb_hard_filter, :]
logger.info(f'Cells left after filtering based on pct_counts_ribo: {adata.shape[0]}')

#########################################################3
# FILTER CELLS BASED ON (SPLICED) COUNTS AND GENES
#########################################################3

# Filter cells
# Based on spliced. Filtered on spliced and unspliced
sc.pp.filter_cells(adata, min_genes = min_genes)

logger.info(f'Cells left after filtering based on genes: {adata.shape[0]}')

sc.pp.filter_cells(adata, min_counts = min_counts_cells)

logger.info(f'Cells left after filtering based on (spliced) counts: {adata.shape[0]}')

#########################################################3
# FILTER GENES BASED ON (SPLICED) COUNTS AND CELLS
#########################################################3

# KEEP GENES BEFORE
mask = adata.var_names.isin(gene_list)
keep_genes_data = adata[:, mask]
adata = adata[:, ~mask]

# Filter genes
sc.pp.filter_genes(adata, min_cells = min_cells)

logger.info(f'Genes left after filtering based on cells: {adata.shape[1]}')

sc.pp.filter_genes(adata, min_counts = min_counts_genes)

logger.info(f'Genes left after filtering based on (spliced) counts: {adata.shape[1]}') 

# KEEP GENES AFTER
adata_obs = adata.obs
adata = an.concat([adata, keep_genes_data], axis = 1)
adata.obs = adata_obs

#########################################################3
# FILTER GENES BASED ON (UNSPLICED & SHARED) COUNTS AND CELLS
#########################################################3

# KEEP GENES BEFORE
mask = adata.var_names.isin(gene_list)
keep_genes_data = adata[:, mask]
adata = adata[:, ~mask]

# Filter genes
    # shared: based on spliced and unspliced. Filtered on spliced and unspliced
    # _: based on spliced. Filtered on spliced and unspliced
    # _u: based on unspliced. Filtered on spliced and unspliced
scv.pp.filter_genes(
adata, 
min_shared_counts = min_shared_counts_genes,
min_shared_cells = min_shared_cells,
min_cells_u = min_cells_u,
min_counts_u = min_counts_genes_u
)

# KEEP GENES AFTER
adata_obs = adata.obs
adata = an.concat([adata, keep_genes_data], axis = 1)
adata.obs = adata_obs

logger.info(f'Genes left after filtering based on (unspliced & shared) counts and cells: {adata.shape[1]}') 

#########################################################3
# FILTER CELLS BASED ON DOUBLETS
#########################################################3

# DOUBLETS
# (velo): Based on spliced. Filtered on spliced and unspliced
sc.external.pp.scrublet(adata, 
    #threshold = threshold_doublets
    )

logger.info(f'Predicted Number of doublets: {adata.obs.predicted_doublet.sum()}')

adata = adata[~adata.obs['predicted_doublet'],:]

logger.info(f'Cells left after filtering based on doublets: {adata.shape[0]}')

# Save counts data
adata.layers['counts'] = adata.X.copy()

#########################################################3
# FINAL COUNT
#########################################################3

logger.info(f"\nAFTER FILTERING\n")

# GENERAL
logger.info(f'Number of cells: {adata.shape[0]}')
logger.info(f'Number of genes: {adata.shape[1]}')
logger.info(f'Total spliced counts: {adata.X.sum()}')
logger.info(f'Total unspliced counts: {adata.layers["unspliced"].sum()}')
logger.info(f'Total counts: {adata.X.sum() + adata.layers["unspliced"].sum()}')
logger.info(f"Sample counts: {adata.obs['sample'].value_counts()}")

# SPLICED
# Counts cells
logger.info(f'Average number of counts per cell (spliced): {adata.obs["total_counts"].mean()}')
# Number of genes per cell
logger.info(f'Average number of genes per cell (spliced): {adata.obs["n_genes_by_counts"].mean()}')
# Counts genes
logger.info(f'Average number of counts per gene (spliced): {adata.var["total_counts"].mean()}')
# Number of cells per gene
logger.info(f'Average number of cells per gene (spliced): {adata.var["n_cells_by_counts"].mean()}')

# UNSPLICED
# Counts cells
logger.info(f'Average number of counts per cell (unspliced): {adata.obs["total_counts_u"].mean()}')
# Number of genes per cell
logger.info(f'Average number of genes per cell (unspliced): {adata.obs["n_genes_by_counts_u"].mean()}')
# Counts genes
logger.info(f'Average number of counts per gene (unspliced): {adata.var["total_counts_u"].mean()}')
# Number of cells per gene
logger.info(f'Average number of cells per gene (unspliced): {adata.var["n_cells_by_counts_u"].mean()}')

# CHECK GENES
for gene in gene_list:
    logger.info(f'{gene} still in data: {gene in adata.var_names}')

#########################################################3
# WRITE COUNT MATRIX
#########################################################3

adata.write_h5ad(snakemake.output.adata)

