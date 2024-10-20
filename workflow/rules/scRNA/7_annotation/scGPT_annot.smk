#########################################################3
# PARAMETERS
#########################################################3

if config["index_scGPT_annot"] == "CellXGene":
    index_dir = f"{bf}ref/scGPT/faiss_idx_CellXGene/"
if config["model_scGPT_annot"] == "whole-human":
    model_dir = f"{bf}ref/scGPT/model_whole-human/"

#########################################################3
# RULES
#########################################################3

rule scGPT_annot:
    input:
        adata = f"{bp}/{bm}6_integration{bm_int}/count_matrix_int.h5ad",
        index_dir = index_dir,
        model_dir = model_dir
    output:
        adata_preannot = f"{bp}/{bm}7_annotation{bm_an}/count_matrix_preannot.h5ad",
        adata_annot = f"{bp}/{bm}7_annotation{bm_an}/count_matrix_annot.h5ad"
    threads: 8
    log: f"{bp}/{bm}7_annotation{bm_an}/scGPT.log"
    singularity: "docker://xueerchen/scgpt:0.1.7"
    conda: conda_env_scGPT
    group: "scGPT_annot"
    params:
        faiss_method = config["faiss_method"],
        n_neighbors_faiss = config["n_neighbors_faiss"],
    resources:
        mem_mb=16000
    script:
        "../../../scripts/scRNA/7_annotation/scGPT_annot.py"

rule scGPT_annot_plots:
    input:
        adata_preannot = f"{bp}/{bm}7_annotation{bm_an}/count_matrix_preannot.h5ad",
        adata_annot = f"{bp}/{bm}7_annotation{bm_an}/count_matrix_annot.h5ad",
    output:
        dim_red_preannot_plot = f"{bp}/{bm}7_annotation{bm_an}/dim_red_preannot.png",
        dim_red_annot_plot = f"{bp}/{bm}7_annotation{bm_an}/dim_red_annot.png",
    threads: 8
    log: f"{bp}/{bm}7_annotation{bm_an}/scGPT_annot_plots.log"
    singularity: "docker://xueerchen/scgpt:0.1.7"
    conda: conda_env_scGPT
    group: "scGPT_annot"
    resources:
        mem_mb=16000
    script:
        "../../../scripts/scRNA/7_annotation/annotation_plots.py"