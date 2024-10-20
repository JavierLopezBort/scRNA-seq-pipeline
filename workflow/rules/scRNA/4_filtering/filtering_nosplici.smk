#########################################################3
# RULES
#########################################################3

rule filtering_nosplici:
    input:
        adata = f"{bp}/{bm}3_count_matrix{bm_cm}/raw_count_matrix.h5ad"
    output:
        adata = f"{bp}/{bm}4_filtering{bm_ft}/count_matrix.h5ad"
    threads: 8
    log: f"{bp}/{bm}4_filtering{bm_ft}/filtering.log"
    conda: conda_env
    group: "filtering_nosplici"
    params:
        gene_list = config["gene_list"],
        # mt and rb genes
        mt_hard_filter = config['mt_hard_filter'],
        rb_hard_filter = config['rb_hard_filter'],
        # Filter cells
        min_genes = config["min_genes"],
        min_counts_cells = config["min_counts_cells"],
        # Filter genes 
        # (velo) Based on spliced. Filtered on spliced and unspliced
        min_cells = config["min_cells"],
        min_counts_genes = config["min_counts_genes"],
        # (velo) Based on uspliced. Filtered on spliced and unspliced
        min_cells_u = config["min_cells_u"],
        min_counts_genes_u = config["min_counts_genes_u"],
        # (velo) Based on spliced and unspliced. Filtered on spliced and unspliced
        min_shared_cells = config["min_shared_cells"],
        min_shared_counts_genes = config["min_shared_counts_genes"],
        # Threshold doublets
        threshold_doublets = config["threshold_doublets"],
    resources:
        mem_mb=42000
    script:
        f"../../../scripts/scRNA/4_filtering/filtering_nosplici.py"

rule filtering_nosplici_plots:
    input:
        adata=f"{bp}/{bm}4_filtering{bm_ft}/count_matrix.h5ad"
    output:
        violin_qc_plots = f"{bp}/{bm}4_filtering{bm_ft}/violin_qc_plots.png",
        highest_expr_genes_plot = f"{bp}/{bm}4_filtering{bm_ft}/highest_expr_genes_plot.png",
        n_genes_table = f"{bp}/{bm}4_filtering{bm_ft}/n_genes_table.csv",
        counts_plot = f"{bp}/{bm}4_filtering{bm_ft}/counts.png",
        violin_genes_plot = f"{bp}/{bm}4_filtering{bm_ft}/violin_genes.png",
        #proportions_fig = f"{bp}/{bm}4_filtering{bm_ft}/proportions_fig.png",
        doublets_plot = f"{bp}/{bm}4_filtering{bm_ft}/doublets.png",
    threads: 8
    log: f"{bp}/{bm}4_filtering{bm_ft}/pyroe_plots.log"
    conda: conda_env
    group: "filtering_nosplici"
    params:
        gene_list = config["gene_list"],
    resources:
        mem_mb=48000
    script:
        "../../../scripts/scRNA/4_filtering/filtering_nosplici_plots.py"