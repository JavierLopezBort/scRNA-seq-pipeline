#########################################################3
# RULES
#########################################################3

rule ivynatal_annot:
    input:
        adata=f"{bp}/{bm}6_integration{bm_int}/count_matrix_int.h5ad",
    output:
        adata = f"{bp}/{bm}7_annotation{bm_an}/count_matrix_annot.h5ad"
    threads: 8
    log: f"{bp}/{bm}7_annotation{bm_an}/clustering_annot.log"
    conda: conda_env
    group: "ivynatal_annot"
    params:
        threshold = config["threshold_an"]
    resources:
        mem_mb=16000
    script:
        "../../../scripts/scRNA/7_annotation/ivynatal_annot.py"

rule ivynatal_annot_plots:
    input:
        adata = f"{bp}/{bm}7_annotation{bm_an}/count_matrix_annot.h5ad"
    output:
        dim_red_plot = f"{bp}/{bm}7_annotation{bm_an}/dim_red.png"
    threads: 8
    log: f"{bp}/{bm}7_annotation{bm_an}/clustering_annot_plots.log"
    conda: conda_env
    group: "ivynatal_annot"
    params:
        threshold = config["threshold_an"]
    resources:
        mem_mb=16000
    script:
        "../../../scripts/scRNA/7_annotation/ivynatal_annot_plots.py"