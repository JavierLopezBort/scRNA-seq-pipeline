#########################################################3
# RULES
#########################################################3

rule pyroe:
    input:
        quant_json = expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/quant.json", sample = sns),
        quants_mat = expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/alevin/quants_mat.mtx", sample = sns),
        quants_mat_col = expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/alevin/quants_mat_cols.txt", sample = sns),
        quants_mat_row = expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/alevin/quants_mat_rows.txt", sample = sns),
        mt_genes = config['mt_genes'],
        rrna_genes = config['rrna_genes'],
        frys = expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/", sample = sns)
    output:
        adata=f"{bp}/{bm}3_count_matrix{bm_cm}/raw_count_matrix.h5ad"
    threads: 8
    log: f"{bp}/{bm}3_count_matrix{bm_cm}/raw_count_matrix.log"
    conda: conda_env
    group: "pyroe"
    params:
        sample_ids = sns,
        study = study,
        specie = specie,
        gender_dict = gender_dict,
        cell_type_dict = cell_type_dict,
        gene_list = config["gene_list"],
    resources:
        mem_mb=48000
    script:
        "../../../scripts/scRNA/3_count_matrix/pyroe.py"

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
    log: f"{bp}/{bm}3_count_matrix{bm_cm}/pyroe_plots.log"
    conda: conda_env
    group: "pyroe"
    params:
        gene_list = config["gene_list"],
    resources:
        mem_mb=48000
    script:
        "../../../scripts/scRNA/3_count_matrix/count_matrix_plots.py"