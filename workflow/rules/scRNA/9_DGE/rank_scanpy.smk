
rule rank_scanpy:
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
    group: "rank_scanpy"
    resources:
        mem_mb=48000
    script:
        "../../../scripts/scRNA/9_DGE/rank_scanpy.py"

rule rank_scanpy_plots:
    input:
        adata = f"{bp}/{bm}9_DGE{bm_dge}/count_matrix_rank.h5ad",
        rank_genes_table = f"{bp}/{bm}9_DGE{bm_dge}/rank_genes.csv"
    output:
        dim_red_rank = f"{bp}/{bm}9_DGE{bm_dge}/dim_red_rank.png",
        #rank_plot = f"{bp}/{bm}9_DGE{bm_dge}/rank_plot.png",
        #rank_dot_plot = f"{bp}/{bm}9_DGE{bm_dge}/rank_dot_plot.png",
    threads: 8
    log: f"{bp}/{bm}9_DGE{bm_dge}/rank_genes_plots.log"
    conda: conda_env
    resources:
        mem_mb=48000
    group: "rank_scanpy"
    params:
        #n_rank_genes = config['n_rank_genes'],
        #groups_plot_rank_genes = config['groups_plot_rank_genes']
        cat_var_rank = config['cat_var_rank'],
    script:
        "../../../scripts/scRNA/9_DGE/rank_scanpy_plots.py"

