#########################################################3
# LOAD LIBRARIES
#########################################################3

import scanpy as sc
import pandas as pd
import vpolo.alevin.parser as vp
import anndata as an
import json

from logging_func import handle_exception, setup_logger

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

logger = setup_logger(snakemake.log[0])
fbarcodes_folder = snakemake.input.fbarcodes_folder
threshold = snakemake.params.threshold

#########################################################3
# CREATE FBARCODES DICT
#########################################################3

fbarcodes_folder = str(fbarcodes_folder)

fbarcodes_df = vp.read_quants_bin(fbarcodes_folder)

genes_list = fbarcodes_df.columns.tolist()

genes_dict = {idx: gene for idx, gene in enumerate(genes_list)}

fbarcodes_df_norm = fbarcodes_df / fbarcodes_df.sum()
fbarcodes_df_norm_norm = fbarcodes_df_norm.div(fbarcodes_df_norm.sum(axis=1), axis=0)
fbarcodes_df_filt = fbarcodes_df_norm_norm[fbarcodes_df_norm_norm.gt(threshold).any(axis=1)]

fbarcodes_df_norm_norm.to_csv(snakemake.output.fbarcodes_table)
fbarcodes_df_filt.to_csv(snakemake.output.fbarcodes_table_filt)

logger.info(f"Number of barcode cells before filtering: {fbarcodes_df_norm_norm.shape[0]}")
logger.info(f"Number of barcode cells after filtering: {fbarcodes_df_filt.shape[0]}")

fbarcodes_class = fbarcodes_df_filt.apply(lambda x: genes_dict[x.argmax()], axis=1)
fbarcodes_dict = fbarcodes_class.to_dict()

#########################################################3
# WRITE FBARCODES DICTs
#########################################################3

with open(snakemake.output.fbarcodes_dict, 'w') as json_file:
    json.dump(fbarcodes_dict, json_file, indent=4)
