#########################################################3
# RULES
#########################################################3

if config["fbarcodes_study"] == "ivn_scRNA":
    fbarcodes_folder = f"{bf}ivn_scRNA/2_alignment{bm_aln}/poolGEX_fb/alevin_output"

rule fbarcodes:
    input:
        quants_mat = f"{fbarcodes_folder}/alevin/quants_mat.gz",
        quants_mat_col = f"{fbarcodes_folder}/alevin/quants_mat_cols.txt",
        quants_mat_row = f"{fbarcodes_folder}/alevin/quants_mat_rows.txt",
        quants_tier_mat = f"{fbarcodes_folder}/alevin/quants_tier_mat.gz",
        fbarcodes_folder = f"{fbarcodes_folder}/"
    output:
        fbarcodes_dict = f"{bp}/{bm}3_count_matrix{bm_cm}/fbarcodes_dict.json",
        fbarcodes_table = f"{bp}/{bm}3_count_matrix{bm_cm}/fbarcodes_table.csv",
        fbarcodes_table_filt = f"{bp}/{bm}3_count_matrix{bm_cm}/fbarcodes_table_filt.csv"
    threads: 8
    log: f"{bp}/{bm}3_count_matrix{bm_cm}/fbarcodes.log"
    conda: conda_env_fbarcodes
    group: "pyroe_fb"
    params:
        threshold = config["threshold"]
    resources:
        mem_mb=48000
    script:
        "../../../scripts/scRNA/3_count_matrix/fbarcodes.py"

# Sample parameters are defined here and not in sample_collection because in sample collection the sample is single,
# and the new parameters are added to the samples got from barcodes, basically they are got from different sources,
# so to be more correct define here.

rule pyroe_fb:
    input:
        fbarcodes_dict = f"{bp}/{bm}3_count_matrix{bm_cm}/fbarcodes_dict.json",

        quant_json = expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/quant.json", sample = sns),
        #collate_json = expand(f"{bp}alevin_fry/{{sample}}/alevin/collate.json", sample=sns),
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
    group: "pyroe_fb"
    params:
        study = study,
        specie = specie,
        gender_dict = gender_dict,
        cell_type_dict = cell_type_dict,
        gene_list = config["gene_list"],
    resources:
        mem_mb=48000
    script:
        "../../../scripts/scRNA/3_count_matrix/pyroe_fb.py"

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
    log: f"{bp}/{bm}3_count_matrix{bm_cm}/pyroe_fb_plots.log"
    conda: conda_env
    group: "pyroe_fb"
    params:
        gene_list = config["gene_list"],
    resources:
        mem_mb=48000
    script:
        "../../../scripts/scRNA/3_count_matrix/count_matrix_plots.py"