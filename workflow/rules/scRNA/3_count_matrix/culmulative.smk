rule culmulative:
    input:
        adata_list = expand(f"{bp}/0_culm/{{study}}.h5ad", study = config["culm_studies"]),
        mt_genes = config['mt_genes'],
        rrna_genes = config['rrna_genes'],
    output:
        adata = f"{bp}/{bm}3_count_matrix{bm_cm}/raw_count_matrix.h5ad"
    threads: 8
    log: f"{bp}/{bm}3_count_matrix{bm_cm}/raw_count_matrix.log"
    conda: conda_env
    group: "culmulative"
    params:
        join_method = config["join_method"],
        gene_list = config["gene_list"],
    resources:
        mem_mb=16000
    script:
        f"../../../scripts/scRNA/3_count_matrix/culmulative.py"

rule count_matrix_plots:
    input:
        adata=f"{bp}/{bm}3_count_matrix{bm_cm}/raw_count_matrix.h5ad"
    output:
        violin_qc_plots = f"{bp}/{bm}3_count_matrix{bm_cm}/violin_qc_plots.png",
        highest_expr_genes_plot = f"{bp}/{bm}3_count_matrix{bm_cm}/highest_expr_genes_plot.png",
        n_genes_table = f"{bp}/{bm}3_count_matrix{bm_cm}/n_genes_table.csv",
        counts_plot = f"{bp}/{bm}3_count_matrix{bm_cm}/counts.png",
        violin_genes_plot = f"{bp}/{bm}3_count_matrix{bm_cm}/violin_genes.png",
        proportions_fig = f"{bp}/{bm}3_count_matrix{bm_cm}/proportions_fig.png"
    threads: 8
    log: f"{bp}/{bm}3_count_matrix{bm_cm}/culmulative_plots.log"
    conda: conda_env
    group: "culmulative"
    params:
        gene_list = config["gene_list"],
    resources:
        mem_mb=48000
    script:
        "../../../scripts/scRNA/3_count_matrix/count_matrix_plots.py"