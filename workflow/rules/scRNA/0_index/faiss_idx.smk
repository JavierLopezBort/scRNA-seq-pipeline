#########################################################3
# PARAMETERS
#########################################################3

index_dir = {bf}{config["index_dir_0"]}
model_dir = {bf}{config["model_dir_0"]}

#########################################################3
# RULES
#########################################################3

rule faiss_index:
    input:
        index_dir = index_dir,
        model_dir = model_dir
    output:
        index_faiss = f"{index_dir}/index.faiss"
    threads: 8
    log: f"{index_dir}/faiss_index.log"
    conda: conda_env_scGPT
    group: "faiss_index"
    params:
        faiss_index_method = config["faiss_index_method"],
    resources:
        mem_mb=16000
    script:
        "../../../scripts/scRNA/0_index/faiss_index.py"