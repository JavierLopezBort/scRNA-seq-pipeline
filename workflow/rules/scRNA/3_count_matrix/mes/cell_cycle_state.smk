#########################################################3
# RULES
#########################################################3

rule cell_cycle_state:
    input:
        adata = f"{bp}/{bm}7_clustering_velo{bm_clvl}/count_matrix_clustering_velo.h5ad"
    output:
        adata = f"{bp}/{bm}8_annotation{bm_an}/count_matrix_annot.h5ad",
    threads: 8
    log: f"{bp}/{bm}8_annotation{bm_an}/cell_cycle_state.log"
    conda: conda_env
    group: "cell_cycle_state" 
    resources:
        mem_mb=16000
    script:
        "../../../scripts/scRNA/8_annotation/cell_cycle_state.py"

rule cell_cycle_state_plots:
    input:
        adata = f"{bp}/{bm}8_annotation{bm_an}/count_matrix_annot.h5ad",
    output:
        dim_red_plot = f"{bp}/{bm}8_annotation{bm_an}/dim_red.png",
        violin_plot = f"{bp}/{bm}8_annotation{bm_an}/violin.png",
    threads: 8
    log: f"{bp}/{bm}8_annotation{bm_an}/cell_cycle_state_plots.log"
    conda: conda_env
    group: "cell_cycle_state"
    resources:
        mem_mb=16000
    params:
        dim_red_plots = config["dim_red_plots_an"],
        cluster_plots = config["cluster_plots_an"],
        cluster_an = config["cluster_an"]
    script:
        "../../../scripts/scRNA/8_annotation/cell_cycle_state_plots.py"