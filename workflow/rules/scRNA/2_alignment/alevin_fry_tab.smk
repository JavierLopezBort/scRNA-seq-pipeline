# scRNA-seq reads from SRA seem to use the UMI as read 1 and the actual sequence as read 2

#########################################################3
# PARAMETERS
#########################################################3

idx = f"{bf}{config['idx_alevin']}"

#########################################################3
# RULES
#########################################################3

rule alevin_fry_tab:
    input:
        mate1 = lambda w : [i.format(read='1') for i in sample_dict[w.sample]],
        mate2 = lambda w : [i.format(read='2') for i in sample_dict[w.sample]],
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
            f"{idx}/reflengths.bin",
            f"{idx}/refseq.bin", # only one diff from the previous one
            f"{idx}/seq.bin",
            f"{idx}/versionInfo.json",
        ]
    output:
        outdir = directory(f"{bp}/2_alignment{bm_aln}/{{sample}}/"),
        quants_mat = f"{bp}/2_alignment{bm_aln}/{{sample}}/alevin/quants_mat.gz",
        quants_mat_col = f"{bp}/2_alignment{bm_aln}/{{sample}}/alevin/quants_mat_cols.txt",
        quants_mat_row = f"{bp}/2_alignment{bm_aln}/{{sample}}/alevin/quants_mat_rows.txt",
        quants_tier_mat = f"{bp}/2_alignment{bm_aln}/{{sample}}/alevin/quants_tier_mat.gz",
        alevin_log = f"{bp}/2_alignment{bm_aln}/{{sample}}/alevin/alevin.log",
        log = f"{bp}/2_alignment{bm_aln}/{{sample}}/alevin_tab.log",
        featureDump = f"{bp}/2_alignment{bm_aln}/{{sample}}/alevin/featureDump.txt",
        alevin_meta_info = f"{bp}/2_alignment{bm_aln}/{{sample}}/aux_info/alevin_meta_info.json",
        ambig_info = f"{bp}/2_alignment{bm_aln}/{{sample}}/aux_info/ambig_info.tsv",
        expected_bias = f"{bp}/2_alignment{bm_aln}/{{sample}}/aux_info/expected_bias.gz",
        fld = f"{bp}/2_alignment{bm_aln}/{{sample}}/aux_info/fld.gz",
        meta_info = f"{bp}/2_alignment{bm_aln}/{{sample}}/aux_info/meta_info.json",
        observed_bias = f"{bp}/2_alignment{bm_aln}/{{sample}}/aux_info/observed_bias.gz",
        observed_bias_3p = f"{bp}/2_alignment{bm_aln}/{{sample}}/aux_info/observed_bias_3p.gz",
        flenDist = f"{bp}/2_alignment{bm_aln}/{{sample}}/libParams/flenDist.txt",
        salmon_quant = f"{bp}/2_alignment{bm_aln}/{{sample}}/logs/salmon_quant.log",
        cmd_info = f"{bp}/2_alignment{bm_aln}/{{sample}}/cmd_info.json",
        lib_format_counts = f"{bp}/2_alignment{bm_aln}/{{sample}}/lib_format_counts.json",
    threads: 16
    conda: conda_env
    resources:
        mem_mb=32000
    params:
        libtype = config["libtype_tab"],
        tech = config["tech_tab"],
        featureStart = config["featureStart"],
        featureLength = config["featureLength"],
        whitelist = config.get('barcode_whitelist', ''),
        expected_cells=config.get('expected_cells', ''),
    log: f"{bp}/2_alignment{bm_aln}/{{sample}}/alevin_tab_fake.log"
    group: "alevin_fry_tab"
    shell:
        """
        salmon alevin -l {params.libtype} -i {input.idx}\
        -1 {input.mate1} -2 {input.mate2}\
        --output {output.outdir} \
        -p {threads}\
        --{params.tech} --featureStart {params.featureStart} --featureLength {params.featureLength} --naiveEqclass\
        > {output.log} 2>&1
        """

rule report:
    input:
        quants_mat_folder = expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/", sample = sns),
        log = expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/alevin_tab.log", sample = sns),
    output:
        report = expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/report.txt", sample = sns)
    threads: 8
    log: expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/report_log.log", sample = sns)
    conda: conda_env_fbarcodes
    resources:
        mem_mb=48000
    group: "alevin_tab"
    script:
        "../../../scripts/scRNA/2_alignment/report_fb.py"


