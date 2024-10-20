import scanpy as sc
import pandas as pd

#########################################################3
# LOAD RAW COUNT MATRIX
#########################################################3

adata = sc.read_h5ad(snakemake.input.adata)
report = open(snakemake.output.report, 'w')

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

# Method
gene_list_corr = snakemake.params.gene_list_corr
gene_list_test = snakemake.params.gene_list_test
n_corr_genes = snakemake.params.n_corr_genes

#########################################################3
# CHECK
#########################################################3

# True/False statements
# correlation
for gene in gene_list_corr:
    logger.info(f'{gene} still in data: {gene in adata.var_names}')

gene_list_corr_statement = all(gene in adata.var_names for gene in gene_list_corr)

if gene_list_corr_statement is False:
    raise ValueError(f"Some correlation genes are not present in the dataset. Please try new ones")
    sys.exit(1)
# test
for gene in gene_list_test:
    logger.info(f'{gene} still in data: {gene in adata.var_names}')

gene_list_test_statement = all(gene in adata.var_names for gene in gene_list_test)

if gene_list_test_statement is False:
    raise ValueError(f"Some tested genes are not present in the dataset. Please try new ones")
    sys.exit(1)

#########################################################3
# EXPRESSION DATA
#########################################################3

expression_data = pd.DataFrame(adata.X.toarray(), columns=adata.var_names, index=adata.obs_names)

#########################################################3
# CORRELATIONS
#########################################################3

method = "pearson"

for gene in gene_list_corr:

    correlations = expression_data.corrwith(expression_data[gene], method = method)
    report.write(f"Gene of interest: {gene}\n")
    report.write(f"Top {n_corr_genes} positive correlated genes: \n{correlations.nlargest(n_corr_genes)}\n")
    report.write(f"Top {n_corr_genes} negative correlated genes: \n{correlations.nsmallest(n_corr_genes)}\n")
    correlations_filt = correlations.loc[gene_list_test]
    report.write(f"Tested genes: \n{correlations_filt}\n")
    report.write("\n")

#########################################################3
# WRITE COUNT MATRIX
#########################################################3

report.close()
