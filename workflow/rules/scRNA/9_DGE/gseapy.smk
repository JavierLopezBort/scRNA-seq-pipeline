rule gseapy:
    input:
        adata = f"{bp}/{bm}6_integration{bm_int}/count_matrix_int.h5ad"
    output:
        gsea_res_table = f"{bp}/{bm}9_DGE{bm_dge}/gsea_res.csv",
        #gsea_res = f"{bp}/{bm}9_DGE{bm_dge}/gsea_res.pkl",
    threads: 8
    log: f"{bp}/{bm}9_DGE{bm_dge}/rank_genes.log"
    conda: conda_env_gseapy
    params:
        cat_var_rank = config['cat_var_rank'],
        gene_set = config["gene_set"],
    group: "gseapy"
    resources:
        mem_mb=48000
    script:
        "../../../scripts/scRNA/9_DGE/gseapy.py"

rule gseapy_plots:
    input:
        gsea_res_table = f"{bp}/{bm}9_DGE{bm_dge}/gsea_res.csv",
        #gsea_res = f"{bp}/{bm}9_DGE{bm_dge}/gsea_res.pkl",
    output:
        volcano_plot = f"{bp}/{bm}9_DGE{bm_dge}/volcano_plot.png",
        gsea_res_table_filt = f"{bp}/{bm}9_DGE{bm_dge}/gsea_res_filt.csv",
    threads: 8
    log: f"{bp}/{bm}9_DGE{bm_dge}/rank_genes_plots.log"
    conda: conda_env_gseapy  
    group: "gseapy"
    resources:
        mem_mb=48000
    script:
        "../../../scripts/scRNA/9_DGE/gseapy_plots.py"


