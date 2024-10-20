import anndata
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import gseapy as gp
import json
import pickle
from logging_func import handle_exception, setup_logger

#########################################################3
# LOAD RAW COUNT MATRIX
#########################################################3

adata = sc.read_h5ad(snakemake.input.adata)

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

logger = setup_logger(snakemake.log[0])

gene_set = snakemake.params.gene_set
cat_var_rank = snakemake.params.cat_var_rank

#########################################################3
# GENE SET
#########################################################3

if gene_set == "/mnt/d/ref/gene_set.json":
    gene_set = json.load(open(gene_set, "r"))
    # A bit hacky
    gene_set['oogonia_meiotic'] = ['MEIOSIN' if x=='BHMG1' else x for x in gene_set['oogonia_meiotic']]

#########################################################3
# GSEA
#########################################################3

gsea_res = gp.gsea(data = adata.to_df().T, # row -> genes, column-> samples
        gene_sets = gene_set,
        cls = adata.obs[cat_var_rank],
        permutation_num=1000,
        permutation_type='phenotype',
        outdir=None,
        method='s2n', # signal_to_noise
        threads= 16)

#########################################################3
# WRITE COUNT MATRIX
#########################################################3
"""
# Save the object
with open(snakemake.output.gsea_res, 'wb') as f:
    pickle.dump(gsea_res, f)
"""
gsea_res.res2d.to_csv(snakemake.output.gsea_res_table)
