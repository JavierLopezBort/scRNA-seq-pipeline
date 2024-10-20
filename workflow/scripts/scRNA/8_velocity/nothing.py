import scanpy as sc

#########################################################3
# LOAD RAW COUNT MATRIX
#########################################################3

adata = sc.read_h5ad(snakemake.input.adata)

#########################################################3
# WRITE COUNT MATRIX
#########################################################3

adata.write_h5ad(snakemake.output.adata)
