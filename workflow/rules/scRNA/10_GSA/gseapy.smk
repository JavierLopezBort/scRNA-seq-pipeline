rule gseapy:
    input:
        rank_file = f"{bp}/{bm}9_DGE{bm_dge}/DGE_filt.csv",
    output:
        gsea_res_table = f"{bp}/{bm}10_GSA{bm_gsa}/gsea_res.csv",
        #gsea_res = f"{bp}/{bm}10_GSA{bm_gsa}/gsea_res.pkl",
    threads: 8
    log: f"{bp}/{bm}10_GSA{bm_gsa}/rank_genes.log"
    conda: conda_env_gseapy
    params:
        gene_set = config["gene_set"],
        rank_stats = config["rank_stats"],
    group: "gseapy"
    resources:
        mem_mb=48000
    script:
        "../../../scripts/scRNA/10_GSA/gseapy.py"

rule gseapy_plots:
    input:
        gsea_res_table = f"{bp}/{bm}10_GSA{bm_gsa}/gsea_res.csv",
        #gsea_res = f"{bp}/{bm}10_GSA{bm_gsa}/gsea_res.pkl",
    output:
        volcano_plot = f"{bp}/{bm}10_GSA{bm_gsa}/volcano_plot.png",
        gsea_res_table_filt = f"{bp}/{bm}10_GSA{bm_gsa}/gsea_res_filt.csv",
    threads: 8
    log: f"{bp}/{bm}10_GSA{bm_gsa}/rank_genes_plots.log"
    conda: conda_env_gseapy
    group: "gseapy"
    resources:
        mem_mb=48000
    script:
        "../../../scripts/scRNA/10_GSA/gseapy_plots.py"


