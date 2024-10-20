# Attempt to determine the cluster IDs 
# Build the reference if possible
import pandas as pd
from logging_func import handle_exception, setup_logger
import scanpy as sc
import numpy as np
import json
import numpy as np
from scipy.sparse import csr_matrix

#########################################################3
# LOAD RAW COUNT MATRIX
#########################################################3

adata = sc.read_h5ad(snakemake.input.adata)
#adata_annot = sc.read_h5ad(snakemake.input.adata_annot)

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

logger = setup_logger(snakemake.log[0])
gene_list = snakemake.params.gene_list

sample_corr = snakemake.params.sample_corr
GFP_corr = snakemake.params.GFP_corr
HBB_corr = snakemake.params.HBB_corr

perc_class_error = 9.3 # Percentage of class error
perc_HBB_basal = 0.864
n_020_exp = 104
n_021_exp = 157
n_STRA8_exp = 399
n_GFP_exp = 23
n_020_gfp = 4
n_021_gfp = 4
n_STRA8_gfp = 8
n_GFP_gfp = 277

#########################################################3
# SAMPLE CORRECTION
#########################################################3

if sample_corr:
    
    #########################################################3
    # GET ARRAYS
    #########################################################3

    mask_HEK = (adata.obs["cell_annot"] == "HEK").values
    mask_iPSC = (adata.obs["cell_annot"] == "iPSC").values

    mask_HTO_C0251 = (adata.obs["sample"] == "HTO_C0251").values
    mask_HTO_C0252 = (adata.obs["sample"] == "HTO_C0252").values
    mask_HTO_C0253 = (adata.obs["sample"] == "HTO_C0253").values
    mask_HTO_C0254 = (adata.obs["sample"] == "HTO_C0254").values
    mask_HTO_C0255 = (adata.obs["sample"] == "HTO_C0255").values
    mask_HTO_C0256 = (adata.obs["sample"] == "HTO_C0256").values

    mask_class = (mask_HEK & (mask_HTO_C0251 | mask_HTO_C0252 | mask_HTO_C0253)) | (mask_iPSC & (mask_HTO_C0254 | mask_HTO_C0255 | mask_HTO_C0256))
    mask_missclass = np.logical_not(mask_class)

    X_class = adata.obsm["X_mde"][mask_class,:]
    y_class = (adata.obs['sample'].astype(str).values)[mask_class]
    X_missclass = adata.obsm["X_mde"][mask_missclass,:]
    y_missclass = (adata.obs['sample'].astype(str).values)[mask_missclass]

    num_folds = 5

    stratified_kfold = StratifiedKFold(n_splits=num_folds, shuffle=True, random_state=42)

    #########################################################3
    # TRAIN AND PREDICT
    #########################################################3

    # TRAIN
    acc = 0
    #accuracy = list()

    for idx_train, idx_test in stratified_kfold.split(X = X_class, y = y_class):
        model = svm.SVC()
        model.fit(X_class[idx_train,:], y_class[idx_train])
        y_pred = model.predict(X_class[idx_test,:])
        acc_model = accuracy_score(y_class[idx_test], y_pred)
        logger.info(f"Accuracy: {acc_model}")
        #accuracy.append(acc_model)
        if acc_model > acc:
            best_model = model

    # PREDICT
    y_pred = best_model.predict(X_missclass)

    acc_missclass = accuracy_score(y_missclass, y_pred)
    logger.info(f"Accuracy missclassified: {acc_missclass}")

    sample_uncorrected = adata.obs['sample'].astype(str).values
    #adata.obs['sample_corr'] = adata.obs["sample"]
    adata.obs['sample'][mask_missclass] = y_pred
    sample_corrected = adata.obs['sample'].astype(str).values

    acc_missclass = accuracy_score(sample_uncorrected , sample_corrected)
    logger.info(f"Accuracy missclassified 2: {acc_missclass}")
    
    #adata.obs["sample_pre"] = adata.obs["sample"]
    #adata.obs["sample"] = adata_annot.obs["sample_corr"].combine_first(adata.obs["sample_pre"]).astype('category')

#########################################################3
# GFP CORRECTION
#########################################################3

if GFP_corr:

    boolean = adata.obs["sample"] == "HTO_C0253"
    boolean = boolean.values
    true_indices = np.where(boolean)[0]
    num_to_change = int((perc_class_error/100) * len(true_indices))
    indices_to_change = np.random.choice(true_indices, size=num_to_change, replace=False)
    boolean[indices_to_change] = False

    matrix = adata.X.toarray()
    matrix[boolean, 50008] += matrix[boolean, 50009]
    matrix[boolean, 50009] = 0
    adata.X = csr_matrix(matrix)
    adata.layers["spliced"] = adata.X

    matrix = adata.layers["unspliced"].toarray()
    matrix[boolean, 50008] += matrix[boolean, 50009]
    matrix[boolean, 50009] = 0
    adata.layers["unspliced"] = csr_matrix(matrix)

#########################################################3
# HBB CORRECTION
#########################################################3

if HBB_corr:
    # 253 ####################################
    logger.info("253")
    boolean = adata.obs["sample"] == "HTO_C0253"
    boolean = boolean.values
    true_indices = np.where(boolean)[0]
    n_HTO_C0253 = len(true_indices)
    num_to_change = int((perc_class_error/100) * len(true_indices))
    indices_to_change = np.random.choice(true_indices, size=num_to_change, replace=False)
    boolean[indices_to_change] = False

    boolean_HBB = (adata.X[:,24820].toarray() > 0).flatten()
    boolean_HBB = boolean_HBB & boolean

    logger.info("second checkpoint")
    n_basal = int((perc_HBB_basal / 100) * n_HTO_C0253)
    logger.info(f"n_basal: {n_basal}")
    true_indices = np.where(boolean_HBB)[0]
    logger.info(f"n_HBB total:{len(true_indices)}")
    num_to_change = n_basal
    indices_to_change = np.random.choice(true_indices, size=num_to_change, replace=False)
    boolean_HBB[indices_to_change] = False
    
    true_indices = np.where(boolean_HBB)[0]
    logger.info(f"n_HBB total:{len(true_indices)}")
    n_HBB = len(true_indices)
    n_exp = n_020_exp + n_021_exp + n_STRA8_exp + n_GFP_exp

    n_mRNA_020 = int(round(n_020_exp / n_exp, 3) * n_HBB)
    n_mRNA_021 = int(round(n_021_exp / n_exp, 3) * n_HBB)
    n_mRNA_STRA8 = int(round(n_STRA8_exp / n_exp, 3) * n_HBB)
    n_mRNA_GFP = n_HBB - (n_mRNA_020 + n_mRNA_021 + n_mRNA_STRA8)
    logger.info(true_indices)
    np.random.shuffle(true_indices)

    logger.info(true_indices)
    logger.info(n_mRNA_020)

    idx_mRNA_020 = true_indices[:n_mRNA_020]
    idx_mRNA_021 = true_indices[n_mRNA_020:n_mRNA_020 + n_mRNA_021]
    idx_mRNA_STRA8 = true_indices[n_mRNA_020 + n_mRNA_021:n_mRNA_020 + n_mRNA_021 + n_mRNA_STRA8]
    idx_mRNA_GFP = true_indices[n_mRNA_020 + n_mRNA_021 + n_mRNA_STRA8:]

    matrix = adata.X.toarray()
    matrix[idx_mRNA_020, 50006] += matrix[idx_mRNA_020, 24820]
    matrix[idx_mRNA_020, 24820] = 0
    matrix[idx_mRNA_021, 50007] += matrix[idx_mRNA_021, 24820]
    matrix[idx_mRNA_021, 24820] = 0
    matrix[idx_mRNA_STRA8, 50008] += matrix[idx_mRNA_STRA8, 24820]
    matrix[idx_mRNA_STRA8, 24820] = 0
    matrix[idx_mRNA_GFP, 50009] += matrix[idx_mRNA_GFP, 24820]
    matrix[idx_mRNA_GFP, 24820] = 0
    adata.X = csr_matrix(matrix)
    adata.layers["spliced"] = adata.X

    matrix = adata.layers["unspliced"].toarray()
    matrix[idx_mRNA_020, 50006] += matrix[idx_mRNA_020, 24820]
    matrix[idx_mRNA_020, 24820] = 0
    matrix[idx_mRNA_021, 50007] += matrix[idx_mRNA_021, 24820]
    matrix[idx_mRNA_021, 24820] = 0
    matrix[idx_mRNA_STRA8, 50008] += matrix[idx_mRNA_STRA8, 24820]
    matrix[idx_mRNA_STRA8, 24820] = 0
    matrix[idx_mRNA_GFP, 50009] += matrix[idx_mRNA_GFP, 24820]
    matrix[idx_mRNA_GFP, 24820] = 0
    adata.layers["unspliced"] = csr_matrix(matrix)

    # 252 ####################################
    logger.info("252")
    boolean = adata.obs["sample"] == "HTO_C0252"
    boolean = boolean.values
    true_indices = np.where(boolean)[0]
    n_HTO_C0253 = len(true_indices)
    num_to_change = int((perc_class_error/100) * len(true_indices))
    indices_to_change = np.random.choice(true_indices, size=num_to_change, replace=False)
    boolean[indices_to_change] = False

    boolean_HBB = (adata.X[:,24820].toarray() > 0).flatten()
    boolean_HBB = boolean_HBB & boolean

    n_basal = int((perc_HBB_basal / 100) * n_HTO_C0253)
    true_indices = np.where(boolean_HBB)[0]
    num_to_change = n_basal
    indices_to_change = np.random.choice(true_indices, size=num_to_change, replace=False)
    boolean_HBB[indices_to_change] = False

    true_indices = np.where(boolean_HBB)[0]
    n_HBB = len(true_indices)
    n_gfp = n_020_gfp + n_021_gfp + n_STRA8_gfp + n_GFP_gfp

    n_mRNA_020 = int(round(n_020_gfp / n_gfp, 3) * n_HBB)
    n_mRNA_021 = int(round(n_021_gfp / n_gfp, 3) * n_HBB)
    n_mRNA_STRA8 = int(round(n_STRA8_gfp / n_gfp, 3) * n_HBB)
    n_mRNA_GFP = n_HBB - (n_mRNA_020 + n_mRNA_021 + n_mRNA_STRA8)
    np.random.shuffle(true_indices)

    idx_mRNA_020 = true_indices[:n_mRNA_020]
    idx_mRNA_021 = true_indices[n_mRNA_020:n_mRNA_020 + n_mRNA_021]
    idx_mRNA_STRA8 = true_indices[n_mRNA_020 + n_mRNA_021:n_mRNA_020 + n_mRNA_021 + n_mRNA_STRA8]
    idx_mRNA_GFP = true_indices[n_mRNA_020 + n_mRNA_021 + n_mRNA_STRA8:]

    matrix = adata.X.toarray()
    matrix[idx_mRNA_020, 50006] += matrix[idx_mRNA_020, 24820]
    matrix[idx_mRNA_020, 24820] = 0
    matrix[idx_mRNA_021, 50007] += matrix[idx_mRNA_021, 24820]
    matrix[idx_mRNA_021, 24820] = 0
    matrix[idx_mRNA_STRA8, 50008] += matrix[idx_mRNA_STRA8, 24820]
    matrix[idx_mRNA_STRA8, 24820] = 0
    matrix[idx_mRNA_GFP, 50009] += matrix[idx_mRNA_GFP, 24820]
    matrix[idx_mRNA_GFP, 24820] = 0
    adata.X = csr_matrix(matrix)
    adata.layers["spliced"] = adata.X

    matrix = adata.layers["unspliced"].toarray()
    matrix[idx_mRNA_020, 50006] += matrix[idx_mRNA_020, 24820]
    matrix[idx_mRNA_020, 24820] = 0
    matrix[idx_mRNA_021, 50007] += matrix[idx_mRNA_021, 24820]
    matrix[idx_mRNA_021, 24820] = 0
    matrix[idx_mRNA_STRA8, 50008] += matrix[idx_mRNA_STRA8, 24820]
    matrix[idx_mRNA_STRA8, 24820] = 0
    matrix[idx_mRNA_GFP, 50009] += matrix[idx_mRNA_GFP, 24820]
    matrix[idx_mRNA_GFP, 24820] = 0
    adata.layers["unspliced"] = csr_matrix(matrix)

#########################################################3
# SCORE CELLS PERCENTAGE
#########################################################3

sample_list = list(adata.obs["sample"].unique())
sample_list.sort()
treatments =  ['mRNA_020', 'mRNA_021', 'mRNA_STRA8_Warren', 'mRNA_GFP', 'HBB']

logger.info(f"SCORE CELLS PERCENTAGE\n")
logger.info(f"\n")

data_list_perc = []
for sample in sample_list:
    adata_subset =  adata[adata.obs["sample"] == sample, :]
    n_total = adata_subset.shape[0]
    perc_list = []
    for treatment in treatments:
        mask = (adata_subset[:, treatment].X.toarray() > 0).flatten()
        adata_subset_subset = adata_subset[mask, :]
        n_transfected = adata_subset_subset.shape[0]
        perc_transf = round((n_transfected / n_total)*100, 3)
        logger.info(f"Transfection percentage of treatment {treatment} in sample {sample} is {perc_transf}")
        perc_list.append(perc_transf)
    data_list_perc.append([sample, treatments, perc_list])
    
    logger.info(f"\n")

logger.info(data_list_perc)

#########################################################3
# SCORE CELLS NUMBER
#########################################################3

logger.info(f"SCORE CELLS NUMBER\n")
logger.info(f"\n")

data_list_numb = []
for sample in sample_list:
    adata_subset =  adata[adata.obs["sample"] == sample, :]
    n_total = adata_subset.shape[0]
    perc_list = []
    for treatment in treatments:
        mask = (adata_subset[:, treatment].X.toarray() > 0).flatten()
        adata_subset_subset = adata_subset[mask, :]
        n_transfected = adata_subset_subset.shape[0]
        #perc_transf = round((n_transfected / n_total)*100, 3)
        logger.info(f"Number of cells of treatment {treatment} in sample {sample} is {n_transfected}")
        perc_list.append(n_transfected)
    data_list_numb.append([sample, treatments, perc_list])
    
    logger.info(f"\n")

logger.info(data_list_numb)

#########################################################3
# WRITE OUTPUT
#########################################################3

with open(snakemake.output.data_list_perc, 'w') as json_file:
    json.dump(data_list_perc, json_file, indent=4)

with open(snakemake.output.data_list_numb, 'w') as json_file:
    json.dump(data_list_numb, json_file, indent=4)

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

adata.write_h5ad(snakemake.output.adata)