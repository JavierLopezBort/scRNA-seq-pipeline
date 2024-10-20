#########################################################3
# PARAMETERS
#########################################################3

if config["model_scGPT_annot"] == "whole-human":
    model_dir = f"{bf}ref/scGPT/model_whole-human/"

#########################################################3
# RULES
#########################################################3

rule scGPT_int:
    input:
        adata=f"{bp}/{bm}5_normalization{bm_nm}/count_matrix.h5ad",
        model_dir = model_dir
    output:
        adata_int = f"{bp}/{bm}6_integration{bm_int}/count_matrix_int.h5ad",
        adata = f"{bp}/{bm}6_integration{bm_int}/count_matrix.h5ad",
    threads: 8
    log: f"{bp}/{bm}6_integration{bm_int}/scGPT.log"
    singularity: "docker://xueerchen/scgpt:0.1.7"
    conda: conda_env_scGPT
    group: "scGPT_int"
    params:
        n_neighbours = config['n_neighbours'],
        resolution = config['resolution'],
        cluster = config["cluster"]
    resources:
        mem_mb=16000
    script:
        "../../../scripts/scRNA/6_integration/scGPT_int.py"

rule scGPT_int_plots:
    input:
        adata = f"{bp}/{bm}6_integration{bm_int}/count_matrix_int.h5ad",
    output:
        dim_red_plot = f"{bp}/{bm}6_integration{bm_int}/dim_red.png",
    threads: 8
    log: f"{bp}/{bm}6_integration{bm_int}/scGPT_annot_plots.log"
    singularity: "docker://xueerchen/scgpt:0.1.7"
    conda: conda_env_scGPT
    group: "scGPT_int"
    resources:
        mem_mb=16000
    script:
        "../../../scripts/scRNA/6_integration/scGPT_int_plots.py"

rule scib_metrics:
    input:
        adata_int = f"{bp}/{bm}6_integration{bm_int}/count_matrix_int.h5ad",
        adata = f"{bp}/{bm}6_integration{bm_int}/count_matrix.h5ad",
    output:
        metrics = f"{bp}/{bm}6_integration{bm_int}/metrics.csv",
    threads: 8
    log: f"{bp}/{bm}6_integration{bm_int}/scGPT_int_metrics.log"
    conda: conda_env_scib
    #conda: "7c121388fccff9718eafc4203b5d5431_"
    #conda: "f567193f21e4767287e25891e8a1e95f_"
    #conda: "1ea3981b6ebca85c0a9ec45dfdec7333_"
    group: "scGPT_int"
    params:
        covariates_metrics = config["covariates_metrics"],
        scib_method = config["scib_method_int"],
        cell_inducted = config["cell_inducted_int"]
    resources:
        mem_mb=16000
    script:
        "../../../scripts/scRNA/6_integration/scib_metrics.py"
