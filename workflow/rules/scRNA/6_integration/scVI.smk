#########################################################3
# RULES
#########################################################3

rule scvi:
    input:
        adata=f"{bp}/{bm}5_normalization{bm_nm}/count_matrix.h5ad"
    output:
        adata_int = f"{bp}/{bm}6_integration{bm_int}/count_matrix_int.h5ad",
        adata = f"{bp}/{bm}6_integration{bm_int}/count_matrix.h5ad",
        model_dir = directory(f"{bp}/{bm}6_integration{bm_int}/model/"),
        model = f"{bp}/{bm}6_integration{bm_int}/model/model.pt"
    threads: 8
    log: f"{bp}/{bm}6_integration{bm_int}/scvi.log"
    conda: conda_env_linux
    group: "scvi"
    params:
        n_neighbours = config['n_neighbours'],
        resolution = config['resolution'],
        cluster = config["cluster"],
        # scvi
        n_latent = config["n_latent"],
        batch_size = config["batch_size"],
        cat_covariates = config["cat_covariates"],
        cont_covariates = config["cont_covariates"],
        n_epochs = config["n_epochs"]
    resources:
        mem_mb=16000
    script:
        "../../../scripts/scRNA/6_integration/scVI.py"

rule scvi_plots:
    input:
        adata=f"{bp}/{bm}6_integration{bm_int}/count_matrix_int.h5ad",
        model = f"{bp}/{bm}6_integration{bm_int}/model/model.pt"
    output:
        dim_red_plot = f"{bp}/{bm}6_integration{bm_int}/dim_red_int.png",
        dendrogram_plot = f"{bp}/{bm}6_integration{bm_int}/dendrogram_int.png",
        corr_mat_plot = f"{bp}/{bm}6_integration{bm_int}/corr_mat_int.png",
        elbo_plot = f"{bp}/{bm}6_integration{bm_int}/elbo_plot.png",
        cluster_counts_plot = f"{bp}/{bm}6_integration{bm_int}/cluster_counts.png"
    threads: 8
    log: f"{bp}/{bm}6_integration{bm_int}/scvi_int_plots.log"
    conda: conda_env_linux
    group: "scvi"
    resources:
        mem_mb=16000
    script:
        "../../../scripts/scRNA/6_integration/scVI_plots.py"

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
    group: "scvi"
    params:
        covariates_metrics = config["covariates_metrics"],
        scib_method = config["scib_method_int"],
        cell_inducted = config["cell_inducted_int"]
    resources:
        mem_mb=16000
    script:
        "../../../scripts/scRNA/6_integration/scib_metrics.py"
