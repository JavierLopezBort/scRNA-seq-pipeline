import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt
#from ivysnake.logging import handle_exception, setup_logger
from logging_func import handle_exception, setup_logger
import seaborn as sns

#########################################################3
# LOAD RAW COUNT MATRIX
#########################################################3

adata = sc.read_h5ad(snakemake.input.adata)

#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

logger = setup_logger(snakemake.log[0])

# Violin genes plots
gene_list = snakemake.params.gene_list

highest_expr_genes_plot = snakemake.output.highest_expr_genes_plot
violin_qc_plots = snakemake.output.violin_qc_plots
n_genes_table = snakemake.output.n_genes_table
counts_plot = snakemake.output.counts_plot
violin_genes_plot = snakemake.output.violin_genes_plot
proportions_fig = snakemake.output.proportions_fig
doublets_plot = snakemake.output.doublets_plot

#########################################################3
# QUALITY CONTROL PLOTS
#########################################################3

# VIOLIN PLOTS
fig, ax = plt.subplots(figsize=(9,5), dpi=300, nrows = 2, ncols=5)

sc.pl.violin(
    adata, 
    "total_counts", 
    #groupby=snakemake.params.batch_key, 
    jitter=0.4,
    size = 0, 
    log=True,
    ax=ax[0,0]
)
ax[0,0].set_xticklabels(ax[0,0].get_xticklabels(), rotation=90)

sc.pl.violin(
    adata, 
    "n_genes_by_counts", 
    #groupby=snakemake.params.batch_key, 
    jitter=0.4,
    size = 0,
    ax=ax[0,1]
)
ax[0,1].set_xticklabels(ax[0,1].get_xticklabels(), rotation=90)

sns.violinplot(data=[adata.var["total_counts"]], ax=ax[0,2])

sns.violinplot(data=[adata.var["n_cells_by_counts"]], ax=ax[0,3])

sc.pl.violin(
    adata, 
    "pct_counts_mito",
    size = 0,
    #groupby=snakemake.params.batch_key, 
    jitter=0.4, 
    ax=ax[0,4]
)
ax[0,4].set_xticklabels(ax[0,4].get_xticklabels(), rotation=90)

sc.pl.violin(
    adata, 
    "total_counts_u",
    size = 0,
    #groupby=snakemake.params.batch_key, 
    jitter=0.4, 
    ax=ax[1,0]
)
ax[1,0].set_xticklabels(ax[1,0].get_xticklabels(), rotation=90)

sc.pl.violin(
    adata, 
    "n_genes_by_counts_u",
    size = 0,
    #groupby=snakemake.params.batch_key, 
    jitter=0.4, 
    ax=ax[1,1]
)
ax[1,1].set_xticklabels(ax[1,1].get_xticklabels(), rotation=90)

sns.violinplot(data=[adata.var["total_counts_u"]], ax=ax[1,2])

sns.violinplot(data=[adata.var["n_cells_by_counts_u"]], ax=ax[1,3])

sc.pl.violin(
    adata, 
    "pct_counts_ribo",
    size = 0,
    #groupby=snakemake.params.batch_key, 
    jitter=0.4, 
    ax=ax[1,4]
)
ax[1,4].set_xticklabels(ax[1,4].get_xticklabels(), rotation=90)

plt.tight_layout()
fig.savefig(violin_qc_plots)

# HIGHEST EXPRESSED GENES PLOT
n_high_expr_genes = 30

fig, ax = plt.subplots(figsize=(9,5), dpi=300, ncols=1)
sc.pl.highest_expr_genes(
    adata, 
    n_top = n_high_expr_genes,
    ax=ax
)
fig.savefig(highest_expr_genes_plot)

# Store the n genes
adata.obs['n_genes_by_counts'].describe().to_csv(n_genes_table)

# COUNTS PLOT
nrows_counts = 2
ncols_counts = 3
cat_var = ["sample","study","gender","cell_type","specie"]
fig, ax = plt.subplots(figsize=(20,10), dpi=300, ncols=ncols_counts, nrows=nrows_counts)

i = 0
for row in range(nrows_counts):
    ax[row,0].set_ylabel('N cells')
    for col in range(ncols_counts):
        if i < len(cat_var):
            var_name = cat_var[i]
            var = adata.obs[var_name]
            var_counts = var.value_counts()
            bars = var_counts.plot(kind='bar', ax=ax[row,col])
            ax[row,col].set_xlabel(var_name)
            ax[row,col].tick_params(axis='x', rotation=45)

            # Adding text above each bar
            for bar in bars.containers[0]:
                height = bar.get_height()
                ax[row, col].text(bar.get_x() + bar.get_width() / 2., height, 
                                  f'{height}', ha='center', va='bottom')
        else:
            fig.delaxes(ax[row,col])
        i += 1 

plt.tight_layout()
#fig.subplots_adjust(wspace=0.25, hspace=0)
fig.savefig(counts_plot)

# VIOLIN GENE PLOTS
nrows_violin = 4
ncols_violin = 10

fig, ax = plt.subplots(figsize=(16,10), dpi=300, ncols = ncols_violin, nrows = nrows_violin)

i = 0
for row in range(nrows_violin):
    col = 0
    while col < ncols_violin:
        if i < len(gene_list):
            if gene_list[i] in adata.var_names:
                sc.pl.violin(adata, gene_list[i], size = 0, jitter = 0.4, ax = ax[row,col])
                col += 1
            i += 1
        else:
            fig.delaxes(ax[row,col])
            col += 1
        
plt.tight_layout()
fig.savefig(violin_genes_plot)

# Proportions
scv.pl.proportions(
    adata, 
    groupby = "sample"
)
fig = plt.gcf()
fig.savefig(proportions_fig)

# Doublet score plot
sc.pl.scrublet_score_distribution(adata)
fig = plt.gcf()
fig.savefig(doublets_plot)
