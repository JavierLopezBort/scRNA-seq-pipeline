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

rank_file = snakemake.input.rank_file

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

logger = setup_logger(snakemake.log[0])

gene_set = snakemake.params.gene_set
rank_stats = snakemake.params.rank_stats

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
        
gene_ranking = pd.read_csv(
        rank_file, 
        #sep="\t", 
        #index_col=0
)

logger.info(gene_ranking)

gsea_res = gp.prerank(
        rnk=gene_ranking.loc[:, rank_stats],
        gene_sets=gene_set,
        #threads=snakemake.threads,
        min_size=5,
        max_size=7000,
        permutation_num=1000,
        outdir=None,
        seed=6, 
        verbose=True
)

#########################################################3
# WRITE COUNT MATRIX
#########################################################3
"""
# Save the object
with open(snakemake.output.gsea_res, 'wb') as f:
    pickle.dump(gsea_res, f)
""" 
gsea_res.res2d.to_csv(snakemake.output.gsea_res_table)
