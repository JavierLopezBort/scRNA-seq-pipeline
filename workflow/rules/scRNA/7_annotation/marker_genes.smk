#########################################################3
# RULES
#########################################################3

rule marker_genes:
    input:
        adata=f"{bp}/{bm}6_integration{bm_int}/count_matrix_int.h5ad",
    output:
        adata_preannot = f"{bp}/{bm}7_annotation{bm_an}/count_matrix_preannot.h5ad",
        adata_annot = f"{bp}/{bm}7_annotation{bm_an}/count_matrix_annot.h5ad"
    threads: 8
    log: f"{bp}/{bm}7_annotation{bm_an}/cell_state.log"
    conda: conda_env
    group: "marker_genes"
    params:
        study_marker_genes = config["study_marker_genes"]
    resources:
        mem_mb=16000
    script:
        "../../../scripts/scRNA/7_annotation/marker_genes.py"

rule marker_genes_plots:
    input:
        adata_preannot = f"{bp}/{bm}7_annotation{bm_an}/count_matrix_preannot.h5ad",
        adata_annot = f"{bp}/{bm}7_annotation{bm_an}/count_matrix_annot.h5ad",
    output:
        dim_red_preannot_plot = f"{bp}/{bm}7_annotation{bm_an}/dim_red_preannot.png",
        dim_red_annot_plot = f"{bp}/{bm}7_annotation{bm_an}/dim_red_annot.png",
    threads: 8
    log: f"{bp}/{bm}7_annotation{bm_an}/marker_genes_plots.log"
    conda: conda_env
    group: "marker_genes"
    resources:
        mem_mb=16000
    script:
        "../../../scripts/scRNA/7_annotation/annotation_plots.py"