
rule rank_scvelo:
    input:
        adata = f"{bp}/{bm}6_integration{bm_int}/count_matrix_int.h5ad"
    output:
        adata = f"{bp}/{bm}9_DGE{bm_dge}/count_matrix_rank.h5ad",
        rank_genes_table = f"{bp}/{bm}9_DGE{bm_dge}/rank_genes.csv"
    threads: 8
    log: f"{bp}/{bm}9_DGE{bm_dge}/rank_genes.log"
    conda: conda_env
    params:
        cat_var_rank = config['cat_var_rank'],
        #n_rank_genes = config['n_rank_genes']
    group: "rank_scvelo"
    resources:
        mem_mb=48000
    script:
        "../../../scripts/scRNA/9_DGE/rank_scvelo.py"

rule rank_scvelo_plots:
    input:
        adata = f"{bp}/{bm}9_DGE{bm_dge}/count_matrix_rank.h5ad",
        rank_genes_table = f"{bp}/{bm}9_DGE{bm_dge}/rank_genes.csv"
    output:
        dim_red_rank = f"{bp}/{bm}9_DGE{bm_dge}/dim_red_rank.png",
    threads: 8
    log: f"{bp}/{bm}9_DGE{bm_dge}/rank_genes_plots.log"
    conda: conda_env
    params:
        cat_var_rank = config["cat_var_rank"],
    group: "rank_scvelo"
    resources:
        mem_mb=48000
    script:
        "../../../scripts/scRNA/9_DGE/rank_scvelo_plots.py"


