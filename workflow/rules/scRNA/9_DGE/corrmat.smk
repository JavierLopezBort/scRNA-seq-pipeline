
rule corrmat:
    input:
        adata = f"{bp}/{bm}6_integration{bm_int}/count_matrix_int.h5ad"
    output:
        report = f"{bp}/{bm}9_DGE{bm_dge}/report.txt"
    threads: 8
    log: f"{bp}/{bm}9_DGE{bm_dge}/rank_genes.log"
    conda: conda_env
    params:
        gene_list_corr = config["gene_list_corr"],
        gene_list_test = config["gene_list"],
        n_corr_genes = config["n_corr_genes"]
    group: "corrmat"
    resources:
        mem_mb=48000
    script:
        "../../../scripts/scRNA/9_DGE/corrmat.py"
