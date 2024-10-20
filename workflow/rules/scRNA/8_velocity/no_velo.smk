#########################################################3
# RULES
#########################################################3

rule no_velo:
    input:
        adata = f"{bp}/{bm}6_integration{bm_int}/count_matrix_int.h5ad"
    output:
        adata=f"{bp}/{bm}8_velocity{bm_vl}/count_matrix_velo.h5ad"
    threads: 8
    log: f"{bp}/{bm}8_velocity{bm_vl}/no_velo.log"
    conda: conda_env   
    resources:
        mem_mb=16000
    script:
        "../../../scripts/scRNA/8_velocity/nothing.py"