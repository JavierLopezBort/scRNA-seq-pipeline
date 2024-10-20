import scib as sb
import scanpy as sc
import pandas as pd
from logging_func import handle_exception, setup_logger

#########################################################3
# LOAD RAW COUNT MATRIX
#########################################################3

adata = sc.read_h5ad(snakemake.input.adata)
adata_int = sc.read_h5ad(snakemake.input.adata_int)

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

logger = setup_logger(snakemake.log[0])

covariates = snakemake.params.covariates_metrics
scib_method = snakemake.params.scib_method
cell_inducted = snakemake.params.cell_inducted

#########################################################3
# REANNOTATE
#########################################################3

if scib_method == "standard":
    
    pass
    
elif scib_method == "atlas":
    
    mask = adata.obs["study"] == "Ivynatal"
    adata = adata[~mask, :]
    adata_int = adata_int[~mask, :]

    adata.obs['cell_annot_comb'] = adata.obs['cell_annot_comb'].astype(str)
    adata.obs['cell_annot_comb'] = adata.obs['cell_annot_comb'].replace('OCT', 'SCT')
    adata.obs['cell_annot_comb'] = adata.obs['cell_annot_comb'].replace('PGCLCs', 'PGCs')
    adata.obs['cell_annot_comb'] = pd.Categorical(adata.obs['cell_annot_comb'])

    adata_int.obs['cell_annot_comb'] = adata_int.obs['cell_annot_comb'].astype(str)
    adata_int.obs['cell_annot_comb'] = adata_int.obs['cell_annot_comb'].replace('OCT', 'SCT')
    adata_int.obs['cell_annot_comb'] = adata_int.obs['cell_annot_comb'].replace('PGCLCs', 'PGCs')
    adata_int.obs['cell_annot_comb'] = pd.Categorical(adata_int.obs['cell_annot_comb'])

elif scib_method == "meioticinduction":
    
    adata.obs['cell_annot_comb'] = adata.obs['cell_annot_comb'].astype(str)
    adata.obs['cell_annot_comb'] = adata.obs['cell_annot_comb'].replace('OCT', 'SCT')
    adata.obs['cell_annot_comb'] = adata.obs['cell_annot_comb'].replace('PGCLCs', 'PGCs')
    adata.obs['cell_annot_comb'] = pd.Categorical(adata.obs['cell_annot_comb'])

    adata_int.obs['cell_annot_comb'] = adata_int.obs['cell_annot_comb'].astype(str)
    adata_int.obs['cell_annot_comb'] = adata_int.obs['cell_annot_comb'].replace('OCT', 'SCT')
    adata_int.obs['cell_annot_comb'] = adata_int.obs['cell_annot_comb'].replace('PGCLCs', 'PGCs')
    adata_int.obs['cell_annot_comb'] = pd.Categorical(adata_int.obs['cell_annot_comb'])
    
    adata.obs['cell_annot_comb'] = adata.obs['cell_annot_comb'].astype(str)
    adata.obs['cell_annot_comb'] = adata.obs['cell_annot_comb'].replace('PGCs', cell_inducted)
    adata.obs['cell_annot_comb'] = pd.Categorical(adata.obs['cell_annot_comb'])
    
    adata_int.obs['cell_annot_comb'] = adata_int.obs['cell_annot_comb'].astype(str)
    adata_int.obs['cell_annot_comb'] = adata_int.obs['cell_annot_comb'].replace('PGCs', cell_inducted)
    adata_int.obs['cell_annot_comb'] = pd.Categorical(adata_int.obs['cell_annot_comb'])
    
    mask = adata.obs["cell_annot_comb"] == cell_inducted
    adata = adata[mask, :]
    adata_int = adata_int[mask, :]
    
#########################################################3
# CALCULATE METRICS
#########################################################3

df_list = []

for batch_key in covariates:
    #logger.info(batch_key)
    sil_df = sb.metrics.metrics(
        adata = adata,
        adata_int = adata_int,
        batch_key = batch_key,
        label_key = "cell_annot_comb",
        embed = "X_emb",
        cluster_key = "cluster",
        #ari_ = True,
        #nmi_ = True,
        #silhouette_batch_ = True,
        #silhouette_label_ = True,
        silhouette_batch_label_ = True,
        #silhouette_label_batch_ = True,
        #pcr_ = True,
        #cell_cycle_ = True,
        #hvg_score_ = True,
        #isolated_labels_ = True,
        #graph_conn_ = True,
        #trajectory_ = True, 
        #kBET_ = True,
        #lisi_graph_ = True,
        type_ = "embed",
        organism = "human",
        verbose = True,
        return_all = True,
    )
    df_list.append(sil_df)
    
df = pd.concat(df_list, axis=1)
total_mean = df.stack().mean()

df_mean_col = pd.DataFrame(df.mean(axis=1), columns=['Mean'])
df_list_2 = [df, df_mean_col]
df_2 = pd.concat(df_list_2, axis=1)

df_mean_index = pd.DataFrame([(df.mean(axis=0)).values], columns=(df.mean(axis=0)).index, index = ["Mean"])
df_mean_index.index.name = "cell_annot_comb"
df_list_3 = [df_2, df_mean_index]
df_3 = pd.concat(df_list_3, axis=0)
df_3.iloc[-1, -1] = total_mean

logger.info(df_3)

#########################################################3
# WRITE COUNT MATRIX
#########################################################3

df_3.to_csv(snakemake.output.metrics)
