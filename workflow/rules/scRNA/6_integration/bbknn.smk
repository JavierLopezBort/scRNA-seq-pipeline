#########################################################3
# RULES
#########################################################3

rule bbknn:
    input:
        adata=f"{bp}/{bm}5_normalization{bm_nm}/count_matrix.h5ad"
    output:
        adata_int = f"{bp}/{bm}6_integration{bm_int}/count_matrix_int.h5ad",
        adata = f"{bp}/{bm}6_integration{bm_int}/count_matrix.h5ad",
    threads: 8
    log: f"{bp}/{bm}6_integration{bm_int}/bbknn.log"
    conda: conda_env
    group: "bbknn"
    params:
        n_neighbours = config['n_neighbours'],
        resolution = config['resolution'],
        cluster = config["cluster"],

        n_pcs = config['n_pcs'],
        layer_pca = config["layer_pca_int"],
        cat_covariate = config["cat_covariate"],
    resources:
        mem_mb=16000
    script:
        "../../../scripts/scRNA/6_integration/bbknn.py"

rule bbknn_plots:
    input:
        adata=f"{bp}/{bm}6_integration{bm_int}/count_matrix_int.h5ad"
    output:
        pca_variance_plot = f"{bp}/{bm}6_integration{bm_int}/pca_variance_int.png",
        dim_red_plot = f"{bp}/{bm}6_integration{bm_int}/dim_red_int.png",
        dendrogram_plot = f"{bp}/{bm}6_integration{bm_int}/dendrogram_int.png",
        corr_mat_plot = f"{bp}/{bm}6_integration{bm_int}/corr_mat_int.png"
    threads: 8
    log: f"{bp}/{bm}6_integration{bm_int}/bbknn_plots.log"
    conda: conda_env
    group: "bbknn"
    resources:
        mem_mb=16000
    script:
        "../../../scripts/scRNA/6_integration/bbknn_plots.py"

rule scib_metrics:
    input:
        adata_int = f"{bp}/{bm}6_integration{bm_int}/count_matrix_int.h5ad",
        adata = f"{bp}/{bm}6_integration{bm_int}/count_matrix.h5ad",
    output:
        metrics = f"{bp}/{bm}6_integration{bm_int}/metrics.csv",
    threads: 8
    log: f"{bp}/{bm}6_integration{bm_int}/scvi_metrics.log"
    conda: conda_env_scib
    #conda: "7c121388fccff9718eafc4203b5d5431_"
    #conda: "f567193f21e4767287e25891e8a1e95f_"
    #conda: "1ea3981b6ebca85c0a9ec45dfdec7333_"
    group: "bbknn"
    params:
        covariates_metrics = config["covariates_metrics"],
        scib_method = config["scib_method_int"],
        cell_inducted = config["cell_inducted_int"]
    resources:
        mem_mb=16000
    script:
        "../../../scripts/scRNA/6_integration/scib_metrics.py"