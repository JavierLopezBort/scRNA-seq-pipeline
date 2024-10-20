#########################################################3
# LOAD LIBRARIES
#########################################################3

import scanpy as sc
import pandas as pd
import pyroe
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from logging_func import handle_exception, setup_logger

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

logger = setup_logger(snakemake.log[0])

with open(snakemake.input.fbarcodes_dict, 'r') as json_file:
    fbarcodes_dict = json.load(json_file)

frys = snakemake.input.frys
study = snakemake.params.study
specie = snakemake.params.specie
gender_dict = snakemake.params.gender_dict
cell_type_dict = snakemake.params.cell_type_dict

rb_genes_file = pd.read_csv(snakemake.input.rrna_genes, header = None)
mt_genes_file = pd.read_csv(snakemake.input.mt_genes, header = None)

# Check genes list
gene_list = snakemake.params.gene_list

#########################################################3
# CONCATENATE COUNT MATRIXES
#########################################################3

adatas = []
#cell_counts = []
rformat = "velocity"

for i in frys:
    tdata = pyroe.load_fry(i, output_format = rformat, quiet = False)
    tdata.obs.index = range(tdata.obs.shape[0])
    tdata.obs.index = tdata.obs.index.astype(str)
    adatas.append(tdata)
adata = sc.concat(adatas, index_unique='_') # All .obs and .var are added here, because I checked tdata and it only contained barcodes and the .X

#########################################################3
# ADD SAMPLE PARAMETERS
#########################################################3

# Sample
adata.obs['sample'] = adata.obs['barcodes'].apply(lambda x: fbarcodes_dict.get(x, "None"))
adata.obs['sample'] = pd.Categorical(adata.obs['sample'])
sample_counts = adata.obs['sample'].value_counts()
logger.info(f"Sample counts: {sample_counts}")
# Study
adata.obs["study"] = study
# Specie
adata.obs["specie"] = specie
# Gender
adata.obs["gender"] = adata.obs['sample'].apply(lambda x: gender_dict[x])
# Cell type
adata.obs["cell_type"] = adata.obs['sample'].apply(lambda x: cell_type_dict[x])

#########################################################3
# INITIAL PARAMETERS
#########################################################3

# Add unspliced
# Counts cells
adata.obs["total_counts_u"] = adata.layers["unspliced"].sum(axis=1)
# Number of genes per cell
adata.obs["n_genes_by_counts_u"] = (adata.layers["unspliced"] > 0).sum(axis=1)
# Counts genes
adata.var["total_counts_u"] = adata.layers["unspliced"].sum(axis=0)
# Number of cells per gene
adata.var["n_cells_by_counts_u"] = (adata.layers["unspliced"] > 0).sum(axis=0)

# Add variables
adata.var['ribo'] = adata.var_names.isin(rb_genes_file[0].values)
adata.var['mito'] = adata.var_names.isin(mt_genes_file[0].values)
sc.pp.calculate_qc_metrics(adata, qc_vars=['ribo', 'mito'], log1p=False, inplace=True)

#adata.layers["raw"] = adata.X.copy()

#########################################################3
# INITIAL COUNT
#########################################################3

sample_counts = adata.obs['sample'].value_counts()
logger.info(f"Sample counts: {sample_counts}")

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
# WRITE RAW COUNT MATRIX
#########################################################3

adata.write_h5ad(snakemake.output.adata)
