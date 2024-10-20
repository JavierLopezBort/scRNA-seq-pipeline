# Attempt to determine the cluster IDs 
# Build the reference if possible
import pandas as pd
from logging_func import handle_exception, setup_logger
from cell_cycle_genes import s_genes, g2m_genes
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

#########################################################3
# FILTER GENES LIST
#########################################################3

# Not necessary to filter, the function ignores the genes that are not found
"""
for gene in s_genes:
    if gene in adata.var_names:
        pass
    else:
        s_genes.remove(gene)
        logger.info(f"{gene} not found in adata")

for gene in g2m_genes:
    if gene in adata.var_names:
        pass
    else:
        g2m_genes.remove(gene)
        logger.info(f"{gene} not found in adata")
"""
#########################################################3
# SCORE CELLS
#########################################################3

adata.X = adata.layers["scaled"]
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

#########################################################3
# CALCULATE MEAN EXPRESSION VALUE
#########################################################3

adata_S = adata[adata.obs["phase"] == "S", :]
adata_G2M = adata[adata.obs["phase"] == "G2M", :]

s_genes_expression = list()
g2m_genes_expression = list()

for gene in s_genes:
    expression = np.mean(adata_S[:, gene].X)
    s_genes_expression.append([gene, expression])

for gene in g2m_genes:
    expression = np.mean(adata_G2M[:, gene].X)
    g2m_genes_expression.append([gene, expression])

s_genes_expression_ord = sorted(s_genes_expression, key=lambda x: x[1])
g2m_genes_expression_ord = sorted(g2m_genes_expression, key=lambda x: x[1])

logger.info(f"S genes:\n")
logger.info(s_genes_expression_ord)

logger.info(f"G2M genes:\n")
logger.info(g2m_genes_expression_ord)

#########################################################3
# WRITE COUNT MATRIX
#########################################################3

adata.write_h5ad(snakemake.output.adata)
