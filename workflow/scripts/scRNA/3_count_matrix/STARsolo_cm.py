#########################################################3
# LOAD LIBRARIES
#########################################################3

import scanpy as sc
import pandas as pd
import scipy

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from logging_func import handle_exception, setup_logger

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

logger = setup_logger(snakemake.log[0])

sample_ids = snakemake.params.sample_ids
study = snakemake.params.study
specie = snakemake.params.specie
gender_dict = snakemake.params.gender_dict
cell_type_dict = snakemake.params.cell_type_dict

rb_genes_file = pd.read_csv(snakemake.input.rrna_genes, header = None)
mt_genes_file = pd.read_csv(snakemake.input.mt_genes, header = None)

# Check genes list
gene_list = snakemake.params.gene_list

#########################################################3
# CREATE ADATA
#########################################################3

adatas = []



for sample_id, spliced, unspliced, ambiguous, features, barcodes in zip(sample_ids, snakemake.input.spliced, snakemake.input.unspliced, snakemake.input.ambiguous, snakemake.input.features, snakemake.input.barcodes):
    
    spliced = scipy.io.mmread(spliced) + scipy.io.mmread(ambiguous)
    unspliced = scipy.io.mmread(unspliced)
    
    spliced = spliced.transpose()
    unspliced = unspliced.transpose()
    
    tdata = sc.AnnData(X = spliced)
    tdata.layers["spliced"] = spliced
    tdata.layers["unspliced"] = unspliced
    
    tdata.obs["study"] = study
    tdata.obs["specie"] = specie
    tdata.obs["sample"] = str(sample_id)
    tdata.obs["gender"] = gender_dict[str(sample_id)]
    tdata.obs["cell_type"] = cell_type_dict[str(sample_id)]
    
    gene_ids = pd.read_csv(features, sep='\t', header = None)[0]
    barcodes = pd.read_csv(barcodes, sep='\t', header = None)[0]
    
    tdata.obs["barcodes"] = barcodes.values
    tdata.obs["barcodes"] = tdata.obs["barcodes"].astype(str)
       
    tdata.obs.index = range(tdata.obs.shape[0])
    tdata.obs.index = tdata.obs.index.astype(str)
    
    tdata.var.index = gene_ids.values
    tdata.var.index = tdata.var.index.astype(str)
    tdata.var.index.name = "gene_ids"
    
    adatas.append(tdata)

adata = sc.concat(adatas, index_unique='_') # All .obs and .var are added here, because I checked tdata and it only contained barcodes and the .X

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

logger.info(adata)

#########################################################3
# INITIAL COUNT
#########################################################3

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
# WRITE RAW COUNT MATRIX
#########################################################3

adata.write_h5ad(snakemake.output.adata)
