import scanpy as sc
import scvelo as scv
import pandas as pd
import matplotlib.pyplot as plt
#import diffxpy.api as de
import numpy as np
import gseapy as gp

#########################################################3
# LOAD RAW COUNT MATRIX
#########################################################3

gsea_res_table = pd.read_csv(snakemake.input.gsea_res_table)
"""
with open(snakemake.input.gsea_res, 'rb') as f:
    gsea_res = pickle.load(f)
"""
#########################################################3
# LOAD DATA FROM SNAKEMAKE
#########################################################3

# Output plots
volcano_plot = snakemake.output.volcano_plot

#########################################################3
# RANK GENES TABLE
#########################################################3

stat_sig = "NOM p-val" # "FDR q-val" "NOM p-val"
stat_sig_threshold = 0.05
stat = "NES"
stat_upper_threshold = 0.0001
stat_lower_threshold = -0.0001

gsea_res_table = gsea_res_table[gsea_res_table[stat_sig] <= stat_sig_threshold]
gsea_res_table = gsea_res_table[(gsea_res_table[stat] >= stat_upper_threshold) | (gsea_res_table[stat] <= stat_lower_threshold)]
gsea_res_table = gsea_res_table.sort_values(by=stat, ascending=False)

gsea_res_table.to_csv(snakemake.output.gsea_res_table_filt, index=False)

#########################################################3
# RANK GENES PLOTS
#########################################################3

logger.info(len(gsea_res_table[stat]))
logger.info(len(-np.log10(gsea_res_table[stat_sig])))
logger.info(gsea_res_table[stat])
logger.info(-np.log10(gsea_res_table[stat_sig]))
# Plotting
plt.figure(figsize=(8, 5))
plt.scatter(gsea_res_table[stat], -np.log10(gsea_res_table[stat_sig]), c='blue', label='Pathways', s = 1)
plt.title(f'Volcano plot of {stat} vs {stat_sig}')
plt.xlabel(stat)
plt.ylabel(f"-log10({stat_sig})")
plt.legend()
#plt.xlim([-3, 3])
#plt.ylim([0, 3])
plt.grid(False)
plt.show()
plt.savefig(volcano_plot)

"""
plt.rcParams['figure.dpi'] = 300

i = 0
terms = gsea_res.res2d.Term
gp.gseaplot(
    rank_metric=gsea_res.ranking,
    term=terms[i],
    ofname=snakemake.output[1],
    **gsea_res.results[terms[i]]
)
"""

