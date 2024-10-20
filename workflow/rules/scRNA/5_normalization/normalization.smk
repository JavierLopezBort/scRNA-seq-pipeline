#########################################################3
# RULES
#########################################################3

rule normalization:
    input:
        adata = f"{bp}/{bm}4_filtering{bm_ft}/count_matrix.h5ad"
    output:
        adata = f"{bp}/{bm}5_normalization{bm_nm}/count_matrix.h5ad"
    threads: 8
    log: f"{bp}/{bm}5_normalization{bm_nm}/normalization.log"
    conda: conda_env
    group: "normalization"
    params:
        n_highly_var=config["n_highly_var"],
        gene_list = config["gene_list"],
        hvg = config["hvg"],
        scaling = config["scaling"]
    resources:
        mem_mb=42000
    script:
        "../../../scripts/scRNA/5_normalization/normalization.py"