rule diffxpy:
    input:
        adata = f"{bp}/{bm}6_integration{bm_int}/count_matrix_int.h5ad"
    output:
        adata = f"{bp}/{bm}9_DGE{bm_dge}/count_matrix_DGE.h5ad",
        DGE_table = f"{bp}/{bm}9_DGE{bm_dge}/DGE.csv",
        #test = f"{bp}/{bm}9_DGE{bm_dge}/test.pkl",
    threads: 8
    log: f"{bp}/{bm}9_DGE{bm_dge}/DGE.log"
    conda: conda_env_diffxpy_2
    params:
        cat_var_rank = config['cat_var_rank'],
        cont_covariates_diffxpy = config["cont_covariates_diffxpy"],
        dge_method = config["dge_method"],
        cell_inducted = config["cell_inducted_dge"],
        #n_rank_genes = config['n_rank_genes']
    group: "diffxpy"
    resources:
        mem_mb=48000
    script:
        "../../../scripts/scRNA/9_DGE/diffxpy.py"

rule diffxpy_plots:
    input:
        adata = f"{bp}/{bm}9_DGE{bm_dge}/count_matrix_DGE.h5ad",
        DGE_table = f"{bp}/{bm}9_DGE{bm_dge}/DGE.csv",
        #test = f"{bp}/{bm}9_DGE{bm_dge}/test.pkl",
    output:
        dim_red_rank = f"{bp}/{bm}9_DGE{bm_dge}/dim_red_DGE.png",
        volcano_plot = f"{bp}/{bm}9_DGE{bm_dge}/volcano_plot.png",
        volcano_plot_2 = f"{bp}/{bm}9_DGE{bm_dge}/volcano_plot_2.png",
        DGE_table_filt = f"{bp}/{bm}9_DGE{bm_dge}/DGE_filt.csv",
        #ttest_plot = f"{bp}/{bm}/9_DGE{bm_dge}/ttest_plot.png",
    threads: 8
    log: f"{bp}/{bm}9_DGE{bm_dge}/DGE_plots.log"
    conda: conda_env_diffxpy_2
    params:
        cat_var_rank = config["cat_var_rank"]
    group: "diffxpy"
    resources:
        mem_mb=48000
    script:
        "../../../scripts/scRNA/9_DGE/diffxpy_plots.py"


