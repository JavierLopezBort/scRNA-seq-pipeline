#########################################################3
# RULES
#########################################################3

if config["fbarcodes_study"] == "ivn_scRNA":
    fbarcodes_folder = f"{bf}ivn_scRNA/2_alignment{bm_aln}/poolGEX_fb"

rule STARsolo_cm_fb:
    input:
        matrix_fb = f"{fbarcodes_folder}/Solo.out/Gene/filtered/matrix.mtx",
        features_fb = f"{fbarcodes_folder}/Solo.out/Gene/filtered/features.tsv",
        barcodes_fb = f"{fbarcodes_folder}/Solo.out/Gene/filtered/barcodes.tsv",

        spliced = expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/Solo.out/Velocyto/filtered/spliced.mtx", sample = sns),
        unspliced = expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/Solo.out/Velocyto/filtered/unspliced.mtx", sample = sns),
        ambiguous = expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/Solo.out/Velocyto/filtered/ambiguous.mtx", sample = sns),
        features = expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/Solo.out/Velocyto/filtered/features.tsv", sample = sns),
        barcodes = expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/Solo.out/Velocyto/filtered/barcodes.tsv", sample = sns),
        mt_genes = config['mt_genes'],
        rrna_genes = config['rrna_genes']
    output:
        #fbarcodes_dict = f"{bp}/{bm}3_count_matrix{bm_cm}/fbarcodes_dict.json",
        adata=f"{bp}/{bm}3_count_matrix{bm_cm}/raw_count_matrix.h5ad"
        fbarcodes_table = f"{bp}/{bm}3_count_matrix{bm_cm}/fbarcodes_table.csv",
        fbarcodes_table_filt = f"{bp}/{bm}3_count_matrix{bm_cm}/fbarcodes_table_filt.csv"
    threads: 8
    log: f"{bp}/{bm}3_count_matrix{bm_cm}/fbarcodes.log"
    conda: conda_env_fbarcodes
    group: "pyroe_fb"
    params:
        sample_ids = sns,
        study = study,
        specie = specie,
        gender_dict = gender_dict,
        cell_type_dict = cell_type_dict,
        gene_list = config["gene_list"],
        threshold = config["threshold"]
    resources:
        mem_mb=48000
    script:
        "../../../scripts/scRNA/3_count_matrix/fbarcodes_STAR.py"

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
    log: f"{bp}/{bm}3_count_matrix{bm_cm}/tsv_plots.log"
    conda: conda_env
    group: "tsv"
    params:
        gene_list = config["gene_list"],
    resources:
        mem_mb=48000
    script:
        "../../../scripts/scRNA/3_count_matrix/count_matrix_plots.py"