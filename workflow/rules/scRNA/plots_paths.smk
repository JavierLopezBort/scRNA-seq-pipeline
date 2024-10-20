#########################################################3
# PATHS
#########################################################3

# INDEX
if index_method == "alevin_fry_idx":
    index_plots = [
        f"{bf}{config['idx_alevin']}/"
    ]
elif index_method == "alevin_fry_tab_idx":
    index_plots = [
        f"{bf}{config['idx_alevin']}/"
    ]
elif index_method == "STARsolo_idx":
    index_plots = [
        f"{bf}{config['idx_STAR']}/"
    ]
elif index_method == "faiss_idx":
    index_plots = [
        f"{index_dir}/index.faiss"
    ]

# RAW
if raw_method == "fastqc":
    if paired_end:
        raw_plots = [
            expand(f"{bp}/FASTQ/fastqc{bm_rw}/{{sample}}{read_preffix}{{read}}{read_suffix}.html", sample = sample_trim_wildcard, read = ["1","2"]),
            expand(f"{bp}/FASTQ/fastqc{bm_rw}/{{sample}}{read_preffix}{{read}}{read_suffix}.zip", sample = sample_trim_wildcard, read = ["1","2"])
        ]
    else:
        raw_plots = [
            expand(f"{bp}/FASTQ/fastqc{bm_rw}/{{sample}}.html", sample = sample_trim_wildcard),
            expand(f"{bp}/FASTQ/fastqc{bm_rw}/{{sample}}.zip", sample = sample_trim_wildcard)
        ] 

# TRIMMING
if trimming_method == "trimmomatic":
    if paired_end:
        trimming_plots = [
            expand(f"{bp}/1_trimming{bm_tr}/fastqc/{{sample}}{read_preffix}{{read}}{read_suffix}.html", sample = sample_trim_wildcard, read = ["1","2"]),
            expand(f"{bp}/1_trimming{bm_tr}/fastqc/{{sample}}{read_preffix}{{read}}{read_suffix}.zip", sample = sample_trim_wildcard, read = ["1","2"])
        ]
    else:
        trimming_plots = [
            expand(f"{bp}/1_trimming{bm_tr}/fastqc/{{sample}}.html", sample = sample_trim_wildcard),
            expand(f"{bp}/1_trimming{bm_tr}/fastqc{bm_tr}/{{sample}}.zip", sample = sample_trim_wildcard)
        ]

if trimming_method == "cutadapt":
    if paired_end:
        trimming_plots = [
            expand(f"{bp}/1_trimming{bm_tr}/fastqc/{{sample}}{read_preffix}{{read}}{read_suffix}.html", sample = sample_trim_wildcard, read = ["1","2"]),
            expand(f"{bp}/1_trimming{bm_tr}/fastqc/{{sample}}{read_preffix}{{read}}{read_suffix}.zip", sample = sample_trim_wildcard, read = ["1","2"])
        ]
    else:
        trimming_plots = [
            expand(f"{bp}/1_trimming{bm_tr}/fastqc/{{sample}}.html", sample = sample_trim_wildcard),
            expand(f"{bp}/1_trimming{bm_tr}/fastqc{bm_tr}/{{sample}}.zip", sample = sample_trim_wildcard)
        ] 

if trimming_method == "fastp":
    if paired_end == True and trimming == True:
        trimming_plots = [
            expand(f"{bp}/1_trimming{bm_tr}/{{sample}}.html", sample = sample_trim_wildcard),
            expand(f"{bp}/1_trimming{bm_tr}/{{sample}}.json", sample = sample_trim_wildcard),
            expand(f"{bp}/1_trimming{bm_tr}/fastqc/{{sample}}{read_preffix}{{read}}{read_suffix}.html", sample = sample_trim_wildcard, read = ["1","2"]),
            expand(f"{bp}/1_trimming{bm_tr}/fastqc/{{sample}}{read_preffix}{{read}}{read_suffix}.zip", sample = sample_trim_wildcard, read = ["1","2"])
        ]
    elif paired_end == False:
        trimming_plots = [
            expand(f"{bp}/1_trimming{bm_tr}/{{sample}}.html", sample = sample_trim_wildcard),
            expand(f"{bp}/1_trimming{bm_tr}/{{sample}}.json", sample = sample_trim_wildcard),
            expand(f"{bp}/1_trimming{bm_tr}/fastqc/{{sample}}.html", sample = sample_trim_wildcard),
            expand(f"{bp}/1_trimming{bm_tr}/fastqc/{{sample}}.zip", sample = sample_trim_wildcard)
        ]
    elif trimming == False:
        trimming_plots = [
            expand(f"{bp}/1_trimming{bm_tr}/{{sample}}.html", sample = sample_trim_wildcard),
            expand(f"{bp}/1_trimming{bm_tr}/{{sample}}.json", sample = sample_trim_wildcard)
        ]
            

if trimming_method == "no_trim":
    if paired_end:
        trimming_plots = [
            expand(f"{bp}/1_trimming{bm_tr}/fastqc/{{sample}}{read_preffix}{{read}}{read_suffix}.html", sample = sample_trim_wildcard, read = ["1","2"]),
            expand(f"{bp}/1_trimming{bm_tr}/fastqc/{{sample}}{read_preffix}{{read}}{read_suffix}.zip", sample = sample_trim_wildcard, read = ["1","2"])
        ]
    else:
        trimming_plots = [
            expand(f"{bp}/1_trimming{bm_tr}/fastqc/{{sample}}.html", sample = sample_trim_wildcard),
            expand(f"{bp}/1_trimming{bm_tr}/fastqc{bm_tr}/{{sample}}.zip", sample = sample_trim_wildcard)
        ]       
     
# ALIGNMENT
if alignment_method == "alevin_fry" or alignment_method == "alevin_fry_tab":
    alignment_plots = [
        expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/report.txt", sample = sns)
    ]
elif alignment_method == "STARsolo" or alignment_method == "STARsolo_tab":
    alignment_plots = [
        expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/Log.final.out", sample = sns),
        expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/Log.out", sample = sns),
        expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/Log.progress.out", sample = sns),
        expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/SJ.out.tab", sample = sns),
        expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/Solo.out/Barcodes.stats", sample = sns),
        expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/Solo.out/Velocyto/Features.stats", sample = sns),
        expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/Solo.out/Velocyto/Summary.csv", sample = sns)
    ]

# COUNT_MATRIX
if count_matrix_method == "pyroe" or count_matrix_method == "pyroe_fb" or count_matrix_method == "culmulative" or count_matrix_method == "STARsolo_cm" or count_matrix_method == "STARsolo_cm_nosplici":
    count_matrix_plots = [
        f"{bp}/{bm}3_count_matrix{bm_cm}/violin_qc_plots.png",
        f"{bp}/{bm}3_count_matrix{bm_cm}/highest_expr_genes_plot.png",
        f"{bp}/{bm}3_count_matrix{bm_cm}/n_genes_table.csv",
        f"{bp}/{bm}3_count_matrix{bm_cm}/counts.png",
        f"{bp}/{bm}3_count_matrix{bm_cm}/violin_genes.png",
        f"{bp}/{bm}3_count_matrix{bm_cm}/proportions_fig.png",
    ]

elif count_matrix_method == "ivynatal_correction":
    count_matrix_plots = [
        f"{bp}/{bm}3_count_matrix{bm_cm}/dim_red.png",
        f"{bp}/{bm}3_count_matrix{bm_cm}/barcodes_counts_perc.png",
        f"{bp}/{bm}3_count_matrix{bm_cm}/barcodes_counts_numb.png",
    ]

# FILTERING
if filtering_method == "filtering":
    filtering_plots = [
        f"{bp}/{bm}4_filtering{bm_ft}/violin_qc_plots.png",
        f"{bp}/{bm}4_filtering{bm_ft}/highest_expr_genes_plot.png",
        f"{bp}/{bm}4_filtering{bm_ft}/n_genes_table.csv",
        f"{bp}/{bm}4_filtering{bm_ft}/counts.png",
        f"{bp}/{bm}4_filtering{bm_ft}/violin_genes.png",
        f"{bp}/{bm}4_filtering{bm_ft}/proportions_fig.png",
        f"{bp}/{bm}4_filtering{bm_ft}/doublets.png",
    ]
elif filtering_method == "filtering_nosplici":
    filtering_plots = [
        f"{bp}/{bm}4_filtering{bm_ft}/violin_qc_plots.png",
        f"{bp}/{bm}4_filtering{bm_ft}/highest_expr_genes_plot.png",
        f"{bp}/{bm}4_filtering{bm_ft}/n_genes_table.csv",
        f"{bp}/{bm}4_filtering{bm_ft}/counts.png",
        f"{bp}/{bm}4_filtering{bm_ft}/violin_genes.png",
        #f"{bp}/{bm}4_filtering{bm_ft}/proportions_fig.png",
        f"{bp}/{bm}4_filtering{bm_ft}/doublets.png",
    ]

# NORMALIZATION
if normalization_method == "normalization":
    normalization_plots = [
        f"{bp}/{bm}5_normalization{bm_nm}/count_matrix.h5ad"
    ]

# INTEGRATION
if integration_method == "bbknn" or integration_method == "diffmap":
    integration_plots = [
        f"{bp}/{bm}6_integration{bm_int}/pca_variance_int.png",
        f"{bp}/{bm}6_integration{bm_int}/dim_red_int.png",
        f"{bp}/{bm}6_integration{bm_int}/dendrogram_int.png",
        f"{bp}/{bm}6_integration{bm_int}/corr_mat_int.png",
        f"{bp}/{bm}6_integration{bm_int}/metrics.csv",
        ]
elif integration_method == "scVI":
    integration_plots = [
        f"{bp}/{bm}6_integration{bm_int}/dim_red_int.png",
        f"{bp}/{bm}6_integration{bm_int}/dendrogram_int.png",
        f"{bp}/{bm}6_integration{bm_int}/corr_mat_int.png",
        f"{bp}/{bm}6_integration{bm_int}/elbo_plot.png",
        f"{bp}/{bm}6_integration{bm_int}/cluster_counts.png",
        f"{bp}/{bm}6_integration{bm_int}/metrics.csv",
        ]

elif integration_method == "scGPT_int":
    integration_plots = [
        f"{bp}/{bm}6_integration{bm_int}/dim_red.png",
        f"{bp}/{bm}6_integration{bm_int}/metrics.csv"
        ]

# ANNOTATION
if annotation_method == "ivynatal_annot":
    annotation_plots = [
         f"{bp}/{bm}7_annotation{bm_an}/dim_red.png",
    ]

elif annotation_method == "scGPT_annot":
    annotation_plots = [
         f"{bp}/{bm}7_annotation{bm_an}/dim_red_preannot.png",
         f"{bp}/{bm}7_annotation{bm_an}/dim_red_annot.png",
    ]
elif annotation_method == "marker_genes":
    annotation_plots = [
         f"{bp}/{bm}7_annotation{bm_an}/dim_red_preannot.png",
         f"{bp}/{bm}7_annotation{bm_an}/dim_red_annot.png",
    ]

# VELOCITY
if velocity_method == "scvelo":
    velocity_plots = [
        f"{bp}/{bm}8_velocity{bm_vl}/velocity.png"
    ]
elif velocity_method == "paga":
    velocity_plots = [
        f"{bp}/{bm}8_velocity{bm_vl}/dim_red.png",
        f"{bp}/{bm}8_velocity{bm_vl}/paga.png",
    ]
elif velocity_method == "no_velo":
    velocity_plots = [
        f"{bp}/{bm}8_velocity{bm_vl}/count_matrix_velo.h5ad",
    ]

# DGE
if DGE_method == "diffxpy":
    DGE_plots = [
        f"{bp}/{bm}9_DGE{bm_dge}/dim_red_DGE.png",
        f"{bp}/{bm}9_DGE{bm_dge}/volcano_plot.png",
        f"{bp}/{bm}9_DGE{bm_dge}/volcano_plot_2.png",
        f"{bp}/{bm}9_DGE{bm_dge}/DGE_filt.csv",
        #f"{bp}/{bm}9_DGE{bm_dge}/volcano_plot.png",
        #f"{bp}/{bm}9_DGE{bm_dge}/ttest_plot.png",
    ]
elif DGE_method == "rank_scvi":
    DGE_plots = [
        f"{bp}/{bm}9_DGE{bm_dge}/dim_red_rank.png"
    ]
elif DGE_method == "rank_scvelo":
    DGE_plots = [
        f"{bp}/{bm}9_DGE{bm_dge}/dim_red_rank.png"
    ]
elif DGE_method == "rank_scanpy":
    DGE_plots = [
        f"{bp}/{bm}9_DGE{bm_dge}/dim_red_rank.png"
    ]
elif DGE_method == "corrmat":
    DGE_plots = [
        f"{bp}/{bm}9_DGE{bm_dge}/report.txt"
    ]
elif DGE_method == "gseapy":
    DGE_plots = [
        f"{bp}/{bm}9_DGE{bm_dge}/volcano_plot.png",
        f"{bp}/{bm}9_DGE{bm_dge}/gsea_res_filt.csv",
    ]

# GSA
if GSA_method == "gseapy":
    GSA_plots = [
        f"{bp}/{bm}10_GSA{bm_gsa}/volcano_plot.png",
        f"{bp}/{bm}10_GSA{bm_gsa}/gsea_res_filt.csv",
    ]

