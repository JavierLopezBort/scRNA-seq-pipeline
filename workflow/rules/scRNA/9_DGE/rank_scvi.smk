
rule rank_scvi:
    input:
        adata = f"{bp}/{bm}6_integration{bm_int}/count_matrix_int.h5ad",
        model_dir = f"{bp}/{bm}6_integration{bm_int}/model/"
    output:
        rank_genes_table = f"{bp}/{bm}9_DGE{bm_dge}/rank_genes.csv"
    threads: 8
    log: f"{bp}/{bm}9_DGE{bm_dge}/rank_genes.log"
    conda: conda_env
    params:
        cat_var_rank = config["cat_var_rank"]
    group: "rank_scvi"
    resources:
        mem_mb=48000
    script:
        "../../../scripts/scRNA/9_DGE/rank_scvi.py"

rule rank_scvi_plots:
    input:
        adata = f"{bp}/{bm}6_integration{bm_int}/count_matrix_int.h5ad",
        rank_genes_table = f"{bp}/{bm}9_DGE{bm_dge}/rank_genes.csv"
    output:
        dim_red_rank = f"{bp}/{bm}9_DGE{bm_dge}/dim_red_rank.png",
    threads: 8
    log: f"{bp}/{bm}9_DGE{bm_dge}/rank_genes_plots.log"
    conda: conda_env
    params:
        cat_var_rank = config["cat_var_rank"],
    group: "rank_scvi"
    resources:
        mem_mb=48000
    script:
        "../../../scripts/scRNA/9_DGE/rank_scvi_plots.py"

