rule paga:
    input:
        adata = f"{bp}/{bm}6_integration{bm_int}/count_matrix_int.h5ad"
    output:
        adata = f"{bp}/{bm}8_velocity{bm_vl}/count_matrix_velo.h5ad"
    threads: 8
    log: f"{bp}/{bm}8_velocity{bm_vl}/velocity.log"
    conda: conda_env
    group: "paga"
    params:
        cluster_root = config["cluster_root"]
    resources:
        mem_mb=16000
    script:
        "../../../scripts/scRNA/8_velocity/paga.py"

rule paga_plots:
    input:
        adata = f"{bp}/{bm}8_velocity{bm_vl}/count_matrix_velo.h5ad"
    output:
        dim_red_plot = f"{bp}/{bm}8_velocity{bm_vl}/dim_red.png",
        paga_plot = f"{bp}/{bm}8_velocity{bm_vl}/paga.png",
    threads: 8
    log: f"{bp}/{bm}8_velocity{bm_vl}/velocity_plots.log"
    conda: conda_env
    group: "paga"
    resources:
        mem_mb=16000
    script:
        "../../../scripts/scRNA/8_velocity/paga_plots.py"