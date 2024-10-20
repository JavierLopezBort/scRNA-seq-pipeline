#########################################################3
# RULES
#########################################################3

rule scvelo:
    input:
        adata = f"{bp}/{bm}6_integration{bm_int}/count_matrix_int.h5ad"
    output:
        adata=f"{bp}/{bm}8_velocity{bm_vl}/count_matrix_velo.h5ad"
    threads: 8
    log: f"{bp}/{bm}8_velocity{bm_vl}/scvelo.log"
    conda: conda_env_scvelo
    group: "scvelo"
    resources:
        mem_mb=32000
    script:
        "../../../scripts/scRNA/8_velocity/scvelo.py"

rule scvelo_plots:
    input:
        adata=f"{bp}/{bm}8_velocity{bm_vl}/count_matrix_velo.h5ad"
    output:
        velocity_plot = f"{bp}/{bm}8_velocity{bm_vl}/velocity.png"
    threads: 8
    log: f"{bp}/{bm}8_velocity{bm_vl}/scvelo_plots.log"
    conda: conda_env_scvelo
    group: "scvelo"
    resources:
        mem_mb=32000
    script:
        "../../../scripts/scRNA/8_velocity/scvelo_plots.py"