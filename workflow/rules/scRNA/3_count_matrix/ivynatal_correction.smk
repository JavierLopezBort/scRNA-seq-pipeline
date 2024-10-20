#########################################################3
# RULES
#########################################################3

rule ivynatal_correction:
    input:
        adata = f"{bp}/{bm}3_count_matrix{bm_cm}/pre_raw_count_matrix.h5ad",
        #adata_annot = f"{bp}/{bm}3_count_matrix{bm_cm}/count_matrix_annot.h5ad",
    output:
        adata = f"{bp}/{bm}3_count_matrix{bm_cm}/raw_count_matrix.h5ad",
        data_list_perc = f"{bp}/{bm}8_annotation{bm_an}/data_list_perc.json",
        data_list_numb = f"{bp}/{bm}8_annotation{bm_an}/data_list_numb.json",
    threads: 8
    log: f"{bp}/{bm}3_count_matrix{bm_cm}/ivynatal_correction.log"
    conda: conda_env
    group: "correction"
    params:
        sample_corr = config["sample_corr"],
        GFP_corr = config["GFP_corr"],
        HBB_corr = config["HBB_corr"]
    resources:
        mem_mb=48000
    script:
        "../../../scripts/scRNA/3_count_matrix/ivynatal_correction.py"

rule ivynatal_correction_plots:
    input:
        adata = f"{bp}/{bm}3_count_matrix{bm_cm}/count_matrix_annot.h5ad",
        data_list_perc = f"{bp}/{bm}3_count_matrix{bm_cm}/data_list_perc.json",
        data_list_numb = f"{bp}/{bm}3_count_matrix{bm_cm}/data_list_numb.json"
    output:
        dim_red_plot = f"{bp}/{bm}3_count_matrix{bm_cm}/dim_red.png",
        barcodes_counts_perc_plot = f"{bp}/{bm}3_count_matrix{bm_cm}/barcodes_counts_perc.png",
        barcodes_counts_numb_plot = f"{bp}/{bm}3_count_matrix{bm_cm}/barcodes_counts_numb.png",
    threads: 8
    log: f"{bp}/{bm}3_count_matrix{bm_cm}/correction_plots.log"
    conda: conda_env
    group: "correction"
    resources:
        mem_mb=16000
    script:
        "../../../scripts/scRNA/8_annotation/ivynatal_correction_plots.py"