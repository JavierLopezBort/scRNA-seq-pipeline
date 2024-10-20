#########################################################3
# PARAMETERS
#########################################################3

transcriptome = f"{bf}{config["transcriptome_tsv"]}"
idx = f"{bf}{config['idx_alevin']}"

#########################################################3
# SALMON INDEX TAB
#########################################################3

rule alevin_fry_tab_idx:
    input:
        tab_file = transcriptome,
    output:
        idx = f"{idx}/",
        idx_files = [
            f"{idx}/complete_ref_lens.bin",
            f"{idx}/ctable.bin",
            f"{idx}/ctg_offsets.bin",
            f"{idx}/duplicate_clusters.tsv",
            f"{idx}/info.json",
            f"{idx}/mphf.bin",
            f"{idx}/pos.bin",
            f"{idx}/pre_indexing.log",
            f"{idx}/rank.bin",
            f"{idx}/refAccumLengths.bin",
            f"{idx}/ref_indexing.log",
            f"{idx}/reflengths.bin", # lacking
            f"{idx}/refseq.bin", 
            f"{idx}/seq.bin",
            f"{idx}/versionInfo.json",
        ]
    log: f"{idx}/salmon_index_tab.log"
    conda: conda_env
    shell:
        """
        salmon index -t {input.tab_file} -i {output.idx} --features -k7
        """
