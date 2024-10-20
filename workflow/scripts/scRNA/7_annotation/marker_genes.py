# Attempt to determine the cluster IDs 
# Build the reference if possible
import pandas as pd
from logging_func import handle_exception, setup_logger
import scanpy as sc
import numpy as np
from marker_genes_list import Persio, Persio_meiosis, Hermann, Hermann_meiosis, Murat, Murat_meiosis, Taelman, Taelman_meiosis, Irie, Irie_meiosis, Huang, Huang_meiosis, Seita, Seita_meiosis, Sosa, Sosa_meiosis, Merrick, Merrick_meiosis, oocytes, oocytes_meiosis, embryo_female, embryo_female_meiosis, comb_list

#########################################################3
# LOAD RAW COUNT MATRIX
#########################################################3

adata = sc.read_h5ad(snakemake.input.adata)

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

logger = setup_logger(snakemake.log[0])

study_marker_genes = snakemake.params.study_marker_genes

#########################################################3
# FILTER GENE DICT
#########################################################3

if study_marker_genes == "Persio":
    genes_dict = Persio
    genes_dict_meiosis = Persio_meiosis

elif study_marker_genes == "Hermann":
    genes_dict = Hermann
    genes_dict_meiosis = Hermann_meiosis

elif study_marker_genes == "Murat":
    genes_dict = Murat
    genes_dict_meiosis = Murat_meiosis

elif study_marker_genes == "Taelman":
    genes_dict = Taelman
    genes_dict_meiosis = Taelman_meiosis

elif study_marker_genes == "Irie":
    genes_dict = Irie
    genes_dict_meiosis = Irie_meiosis

elif study_marker_genes == "Huang":
    genes_dict = Huang
    genes_dict_meiosis = Huang_meiosis
    
elif study_marker_genes == "Seita":
    genes_dict = Seita
    genes_dict_meiosis = Seita_meiosis
    
elif study_marker_genes == "Sosa":
    genes_dict = Sosa
    genes_dict_meiosis = Sosa_meiosis

elif study_marker_genes == "Merrick":
    genes_dict = Merrick
    genes_dict_meiosis = Merrick_meiosis
    
elif study_marker_genes == "oocytes":
    genes_dict = oocytes
    genes_dict_meiosis = oocytes_meiosis
    
elif study_marker_genes == "embryo_female":
    genes_dict = embryo_female
    genes_dict_meiosis = embryo_female_meiosis

logger.info('Starting cluster_identification')

for cell_type in genes_dict.keys():
    for gene in genes_dict[cell_type]:
        if gene in adata.var_names:
            pass
        else:
            logger.info(f"{gene} not present")
            genes_dict[cell_type].remove(gene)

for cell_type in genes_dict.keys():
    if len(genes_dict[cell_type]) == 0:
        del genes_dict[cell_type]

#########################################################3
# SCORE CELLS
#########################################################3

cell_type_list = list()

for cell_type in genes_dict.keys():
    sc.tl.score_genes(adata, genes_dict[cell_type], score_name = f"{cell_type}")
    cell_type_list.append(f"{cell_type}")

# Calculate the mean score of each cluster using the cell scores (we have already clustered)
member_df = adata.obs.groupby("cluster")[cell_type_list].mean()

# Step 2: Assign each cluster to the cell type with the highest mean score
cluster_to_cell_type = member_df.idxmax(axis=1)

# Step 3: Assign each individual cell to the cell type of the cluster they belong to
adata.obs['cell_annot'] = adata.obs["cluster"].map(cluster_to_cell_type)
adata.obs['cell_annot_comb'] = adata.obs['cell_annot'].map(comb_list)

logger.info(f'Cells before filtering: {adata.shape[0]}')
adata.write_h5ad(snakemake.output.adata_preannot)

#########################################################3
# FILTER MEIOTIC CELLS
#########################################################3

cell_type_list_meiotic = list()

for cell_type in genes_dict_meiosis.keys():
    if genes_dict_meiosis[cell_type] == "meiotic":
        cell_type_list_meiotic.append(cell_type)
        
logger.info(cell_type_list_meiotic)

mask = adata.obs["cell_annot"].isin(cell_type_list_meiotic)
adata = adata[mask, :]  

logger.info(f'Cells after filtering: {adata.shape[0]}')

#########################################################3
# WRITE COUNT MATRIX
#########################################################3

adata.write_h5ad(snakemake.output.adata_annot)


