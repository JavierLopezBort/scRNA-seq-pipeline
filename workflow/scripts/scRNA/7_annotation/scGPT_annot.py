#########################################################3
# IMPORT PACKAGES
#########################################################3

import sys
import os
from pathlib import Path
from os import PathLike

import numpy as np
import pandas as pd
from scipy.stats import mode
import scanpy as sc
from typing import Union, Tuple, Optional, Dict, Any
import sklearn
import warnings
import h5py

sys.path.insert(0, "../")
import scgpt as scg
import faiss
#from build_atlas_index_faiss import load_index, vote, compute_category_proportion

from tqdm import tqdm

from scGPT_annot_list import meiotic_cells, somatic_cells, comb_list

#########################################################3
# FUNCTIONS
#########################################################3

# Load the faiss index
def load_index(
    index_dir: PathLike,
    use_config_file=True,
    use_gpu=False,
    nprobe=None,
) -> Tuple[faiss.Index, np.ndarray]:
    """
    Load index from disk.

    Args:
        index_dir (PathLike): Path to the directory containing the index files.
        use_config_file (bool, optional): Whether to load the index config file. If True, will load the index config file and use the parameters of gpu, nprobe. Defaults to True.
        use_gpu (bool, optional): Whether to use GPU acceleration. Only used when use_config_file is False. Defaults to False.
        nprobe (int, optional): The nprobe to set if index contains :class:`faiss.IndexIVF`. If None, will set based on the number of clusters. Only used when use_config_file is False. Defaults to None.

    Returns:
        faiss.Index: The loaded index and meta labels.
    """

    index_file = os.path.join(index_dir, "index.faiss")
    meta_file = os.path.join(index_dir, "meta.h5ad")
    index_config_file = os.path.join(index_dir, "index_config.json")

    print(f"Loading index and meta from {index_dir} ...")
    index = faiss.read_index(index_file)
    print(f"Index loaded, num_embeddings: {index.ntotal}")
    with h5py.File(meta_file, "r") as f:
        meta_labels = f["meta_labels"].asstr()[:]
    if use_config_file:
        with open(index_config_file, "r") as f:
            config = json.load(f)
        use_gpu = config["gpu"]
        nprobe = config["nprobe"]

    _auto_set_nprobe(index, nprobe=nprobe)

    if use_gpu:
        res = faiss.StandardGpuResources()
        index = faiss.index_cpu_to_gpu(res, 0, index)

    return index, meta_labels

# Used in the load faiss index
def _auto_set_nprobe(index: faiss.Index, nprobe: int = None) -> Optional[int]:
    """
    Set nprobe for IVF index based on the number of clusters.

    Args:
        index (faiss.Index): The index to set nprobe.
        nprobe (int, optional): The nprobe to set. If None, will set based on the number of clusters. Defaults to None.

    Returns:
        int: The nprobe set.
    """

    # set nprobe if IVF index
    index_ivf = faiss.try_extract_index_ivf(index)
    if index_ivf:
        nlist = index_ivf.nlist
        ori_nprobe = index_ivf.nprobe
        index_ivf.nprobe = (
            nprobe
            if nprobe is not None
            else 16
            if nlist <= 1e3
            else 32
            if nlist <= 4e3
            else 64
            if nlist <= 1.6e4
            else 128
        )
        print(
            f"Set nprobe from {ori_nprobe} to {index_ivf.nprobe} for {nlist} clusters"
        )
        return index_ivf.nprobe


def compute_category_proportion(meta_labels) -> Dict[str, float]:
    """
    Compute the proportion of each cell type in the meta_labels, which can be used for weighted voting in the search.

    Args:
        meta_labels (numpy.ndarray): A 1D array of cell type labels.

    Returns:
        dict: A dictionary containing the proportion of each cell type in the input array.
    """
    unique_labels, counts = np.unique(meta_labels, return_counts=True)
    category_proportion = dict(zip(unique_labels, counts / counts.sum()))
    return category_proportion


def weighted_vote(
    predicts_for_query, cell_type_proportion, return_prob=True
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Use the proportion of each cell type as the weight for voting.

    Args:
        predicts_for_query (numpy.ndarray): A 1D array of the predicted cell types for an individual query.
        cell_type_proportion (dict): A dictionary containing the proportion of each cell type in the meta_labels.
        return_prob (bool, optional): Whether to return the probability of each predicted cell type. Defaults to True.

    Returns:
        numpy.ndarray: A 1D array of the predicted cell types for the input query, weighted by the proportion of each cell type and sorted by the proportion.
        numpy.ndarray: A 1D array of the probability of each predicted cell type. Only returned when return_prob is True.
    """
    unique_labels, counts = np.unique(predicts_for_query, return_counts=True)
    weighted_counts = (
        np.clip(counts - 0.01 * counts.sum(), 0, None)
        * 1e-3
        / np.array([cell_type_proportion[l] for l in unique_labels])
    )  # the -1 is to reduce noise
    weighted_counts = weighted_counts / weighted_counts.sum()
    sorted_idx = np.argsort(weighted_counts)[::-1]
    predicts_for_query = unique_labels[sorted_idx]

    if return_prob:
        return predicts_for_query, weighted_counts[sorted_idx]
    return predicts_for_query


def vote(predicts_for_query, return_prob=True) -> Tuple[np.ndarray, np.ndarray]:
    """
    Majority voting for the predicted cell types.

    Args:
        predicts_for_query (numpy.ndarray): A 1D array of the predicted cell types for an individual query.
        return_prob (bool, optional): Whether to return the probability of each predicted cell type. Defaults to True.

    Returns:
        numpy.ndarray: A 1D array of the predicted cell types for the input query, weighted by the proportion of each cell type and sorted by the proportion.
        numpy.ndarray: A 1D array of the probability of each predicted cell type. Only returned when return_prob is True.
    """
    unique_labels, counts = np.unique(predicts_for_query, return_counts=True)
    weighted_counts = counts / counts.sum()
    sorted_idx = np.argsort(weighted_counts)[::-1]
    predicts_for_query = unique_labels[sorted_idx]

    if return_prob:
        return predicts_for_query, weighted_counts[sorted_idx]
    return predicts_for_query

# Those functions are only used when faiss is not installed
def l2_sim(a, b):
    sims = -np.linalg.norm(a - b, axis=1)
    return sims

def get_similar_vectors(vector, ref, top_k=10):
    # sims = cos_sim(vector, ref)
    sims = l2_sim(vector, ref)

    top_k_idx = np.argsort(sims)[::-1][:top_k]
    return top_k_idx, sims[top_k_idx]

#########################################################3
# IMPORT FROM SNAKEMAKE
#########################################################3

adata = sc.read_h5ad(snakemake.input.adata)
index_dir = Path(snakemake.input.index_dir)
model_dir = Path(snakemake.input.model_dir)

faiss_method = snakemake.params.faiss_method
n_neighbors_faiss = snakemake.params.n_neighbors_faiss

#########################################################3
# FAISS INDEX
#########################################################3

if faiss_method == "faiss_build":

    use_gpu = faiss.get_num_gpus() > 0
    logger.info(f"Use GPU: {use_gpu}")

    index, meta_labels = load_index(
        index_dir = index_dir,
        use_config_file=False,
        use_gpu = use_gpu,
    )
    logger.info(f"Loaded index with {index.ntotal} cells")

elif faiss_method == "faiss_add":

    index = faiss.read_index(f"{index_dir}/index.faiss")

elif faiss_method == "no_faiss":

    adata_ref = sc.read_h5ad(f"{index_dir}/adata_ref.h5ad")

    ref_embed_adata = scg.tasks.embed_data(
        adata_ref,
        model_dir,
        gene_col=gene_col,
        obs_to_save=cell_type_key,  # optional arg, only for saving metainfo
        batch_size=64,
        return_new_adata=True,
    )

    ref_cell_embeddings = ref_embed_adata.X


#########################################################3
# GET THE EMBEDDINGS
#########################################################3

gene_col = 'index'

test_embed_adata = scg.tasks.embed_data(
    adata,
    model_dir,
    gene_col=gene_col,
    #obs_to_save=cell_type_key,  # optional arg, only for saving metainfo
    batch_size=64,
    return_new_adata=True,
)

test_embed = test_embed_adata.X

#########################################################3
# GET FINAL DISTANCES
#########################################################3

k = n_neighbors_faiss

if faiss_method == "faiss_build":

    # test with the first 100 cells
    distances, idx = index.search(test_embed, k)
    predict_labels = meta_labels[idx]

    voting = []
    for preds in tqdm(predict_labels):
        voting.append(vote(preds, return_prob=False)[0])
    voting = np.array(voting)
    adata.obs['cell_annot'] = voting

elif faiss_method == "faiss_add":

    # test with the first 100 cells
    distances, idx = index.search(test_embed, k)

    idx_list=[i for i in range(test_emebd.shape[0])]
    preds = []

    for k in idx_list:
        idx = labels[k]
        pred = ref_embed_adata.obs[cell_type_key][idx].value_counts()
        preds.append(pred.index[0])
    
    voting = []
    for preds in tqdm(predict_labels):
        voting.append(vote(preds, return_prob=False)[0])
    voting = np.array(voting)
    adata.obs['cell_annot'] = voting
    

elif faiss_method == "no_faiss":

    idx_list=[i for i in range(test_emebd.shape[0])]
    preds = []

    for k in idx_list:
        idx, sim = get_similar_vectors(test_emebd[k][np.newaxis, ...], ref_cell_embeddings, k)
        pred = ref_embed_adata.obs[cell_type_key][idx].value_counts()
        preds.append(pred.index[0])
    
    voting = []
    for preds in tqdm(predict_labels):
        voting.append(vote(preds, return_prob=False)[0])
    voting = np.array(voting)
    adata.obs['cell_annot'] = voting

#########################################################3
# ASSIGN CELL TYPE TO CLUSTERS
#########################################################3

grouped = adata.obs.groupby(["cluster", 'cell_annot']).size().reset_index(name='count')
majority_cell_type = grouped.sort_values('count', ascending=False).drop_duplicates("cluster").set_index("cluster")['cell_annot']
adata.obs['cell_annot'] = adata.obs["cluster"].map(majority_cell_type)
adata.obs['cell_annot_comb'] = adata.obs['cell_annot'].map(comb_list)

#########################################################3
# WRITE COUNT MATRIX
#########################################################3

logger.info(f'Cells before filtering: {adata.shape[0]}')
adata.write_h5ad(snakemake.output.adata_preannot)

#########################################################3
# WRITE COUNT MATRIX
#########################################################3

mask = adata.obs["cell_annot"].isin(meiotic_cells)
adata = adata[mask, :]  

logger.info(f'Cells after filtering: {adata.shape[0]}')

adata.write_h5ad(snakemake.output.adata_annot)