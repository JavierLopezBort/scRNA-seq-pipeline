import scanpy as sc
from logging_func import handle_exception, setup_logger
import louvain
import anndata
import pandas as pd

#########################################################3
# LOAD ADATA LIST
#########################################################3

adata_path_list = snakemake.input.adata_list

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

logger = setup_logger(snakemake.log[0])

join_method = snakemake.params.join_method

# Check genes list
gene_list = snakemake.params.gene_list
rb_genes_file = pd.read_csv(snakemake.input.rrna_genes, header = None)
mt_genes_file = pd.read_csv(snakemake.input.mt_genes, header = None)

#########################################################3
# JOIN
#########################################################3

#logger.info(f"file_path_list: {adata_path_list}")

adata_list = []
count_list = []
count = 0

for file_path in adata_path_list:
    adata = sc.read_h5ad(file_path)
    adata_list.append(adata)
    logger.info(f"Study: {str(adata.obs['study'])}")
    logger.info(f'Cells before merging: {adata.shape[0]}')
    logger.info(f'Genes before merging: {adata.shape[1]}')
    logger.info(adata.obs["n_genes_by_counts"].mean())
    count_list.append(adata.obs["n_genes_by_counts"].mean())
    count += 1

logger.info(sum(count_list) / count)

#logger.info(f"adata_list: {adata_list}")

adata = anndata.concat(adata_list, join = join_method, index_unique="_", fill_value = 0)

#########################################################3
# INITIAL PARAMETERS
#########################################################3

# Add spliced
# Counts cells
adata.obs["total_counts"] = adata.X.sum(axis=1)
# Number of genes per cell
adata.obs["n_genes_by_counts"] = (adata.X > 0).sum(axis=1)
# Counts genes
adata.var["total_counts"] = adata.X.sum(axis=0).reshape(-1, 1)
# Number of cells per gene
adata.var["n_cells_by_counts"] = (adata.X > 0).sum(axis=0).reshape(-1, 1)

# Add unspliced
# Counts cells
adata.obs["total_counts_u"] = adata.layers["unspliced"].sum(axis=1)
# Number of genes per cell
adata.obs["n_genes_by_counts_u"] = (adata.layers["unspliced"] > 0).sum(axis=1)
# Counts genes
adata.var["total_counts_u"] = adata.layers["unspliced"].sum(axis=0).reshape(-1, 1)
# Number of cells per gene
adata.var["n_cells_by_counts_u"] = (adata.layers["unspliced"] > 0).sum(axis=0).reshape(-1, 1)
    
# Add variables
adata.var['ribo'] = adata.var_names.isin(rb_genes_file[0].values)
adata.var['mito'] = adata.var_names.isin(mt_genes_file[0].values)
sc.pp.calculate_qc_metrics(adata, qc_vars=['ribo', 'mito'], log1p=False, inplace=True)

#########################################################3
# INITIAL COUNT
#########################################################3

sample_counts = adata.obs['sample'].value_counts()
logger.info(f"Sample counts: {sample_counts}")

logger.info("After merging:\n")
logger.info(f'Number of cells: {adata.shape[0]}')
logger.info(f'Number of genes: {adata.shape[1]}')
logger.info(f'Total spliced counts: {adata.X.sum()}')
logger.info(f'Total unspliced counts: {adata.layers["unspliced"].sum()}')
logger.info(f'Total counts: {adata.X.sum() + adata.layers["unspliced"].sum()}')

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