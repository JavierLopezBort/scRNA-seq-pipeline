#########################################################3
# PARAMETERS
#########################################################3

run = "local" # local, gcp, hpc, computerome, lambda

# SAMPLE
samplex = "culm_annot" 
# "Hermann", "Persio", "Murat", "Murat_chimp", "Huang", "Murat_chicken", "Murat_macaque", "Murat_marmoset", 
# "Murat_mouse", "Reich", "Ge", "Ivynatal", "Ivynatal_fb", "Irie", "Seita", "Smela", "Sosa", "Taelman"
read_preffix = "_R"          # example: "_R"
read_suffix = "_001"         # example: ""
fastq_suffix = ".fastq.gz"   # example : ".fastq.gz"

# BENCHMARK
bm = "/" # Remember the "/"
bm_rw = ""
bm_tr = ""
bm_aln = ""
bm_cm = ""
bm_ft = ""
bm_nm = ""
bm_int = ""
bm_an = ""
bm_vl = ""
bm_dge = ""
bm_gsa = ""

#########################################################3
# PREREQUISITES
#########################################################3

# PATHS
include: "rules/scRNA/paths.smk"
# CONDA ENVS
include: "rules/scRNA/conda_envs.smk"
# METHODS
include: "rules/scRNA/methods.smk"
# SAMPLE BUILDER
include: "rules/scRNA/sample_collection.smk"
# PLOTS PATHS
include: "rules/scRNA/plots_paths.smk"

logger.info(f"bf: {bf}")
logger.info(f"bp {bp}")

#########################################################3
# RULES
#########################################################3

# IDX
#include: f"rules/scRNA/0_index/{index_method}.smk"

# RAW
#include: f"rules/scRNA/0_raw/{raw_method}.smk"

# TRIMMING
#include: f"rules/scRNA/1_trimming/{trimming_method}.smk"

# ALIGNMENT
#include: f"rules/scRNA/2_alignment/{alignment_method}.smk"

# COUNT MATRIX
#include: f"rules/scRNA/3_count_matrix/{count_matrix_method}.smk"

# FILTERING
#include: f"rules/scRNA/4_filtering/{filtering_method}.smk"

# NORMALIZATION
#include: f"rules/scRNA/5_normalization/{normalization_method}.smk"

# INTEGRATION
#include: f"rules/scRNA/6_integration/{integration_method}.smk"

# ANNOTATION
#include: f"rules/scRNA/7_annotation/{annotation_method}.smk"

# VELOCITY
#include: f"rules/scRNA/8_velocity/{velocity_method}.smk"

# DGE
#include: f"rules/scRNA/9_DGE/{DGE_method}.smk"

# GSA
#include: f"rules/scRNA/10_GSA/{GSA_method}.smk"

#########################################################3
# RULE ALL
#########################################################3

rule all:
    input:
        #index_plots,
        #raw_plots,
        #trimming_plots,
        #alignment_plots,
        #count_matrix_plots,
        #filtering_plots,
        #normalization_plots,
        #integration_plots,
        #annotation_plots,
        #velocity_plots,
        #DGE_plots,
        #GSA_plots,