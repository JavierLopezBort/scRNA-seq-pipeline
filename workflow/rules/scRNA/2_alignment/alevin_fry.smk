# scRNA-seq reads from SRA seem to use the UMI as read 1 and the actual sequence as read 2

#########################################################3
# PARAMETERS
#########################################################3

idx = f"{bf}{config['idx_alevin']}"
tx2gene = f"{bf}{config['tx2gene']}"
mt_genes = f"{bf}{config['mt_genes']}"
rrna_genes = f"{bf}{config['rrna_genes']}"

#########################################################3
# RULES
#########################################################3

if config["paired_end"]:
    rule rad_pe:
        input:
            mate1 = lambda w : [i.format(read='1') for i in sample_dict[w.sample]],
            mate2 = lambda w : [i.format(read='2') for i in sample_dict[w.sample]],
            mt_genes = mt_genes,
            rrna_genes = rrna_genes,
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
                f"{idx}/refseq.bin",
                #f"{idx}/reflengths.bin",
                f"{idx}/seq.bin",
                f"{idx}/versionInfo.json"
            ],
            tx2gene = tx2gene,
        output:
            outdir = directory(f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/"),
            map_file = f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/map.rad",
            unmapped = f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/unmapped_bc_count.bin",
            log = f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/rad.log",
            cmd_info = f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/cmd_info.json",
            alevin_log = f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/alevin/alevin.log",
            meta_info = f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/aux_info/meta_info.json",
            libParams = directory(f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/libParams/"),
            salmon_quant = f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/logs/salmon_quant.log",
        #threads: 16
        conda: conda_env_salmon
        #resources:
            #mem_mb=40000
        params:
            libtype = config["libtype_fry"],
            tech = config["tech_fry"],
            whitelist = config.get('barcode_whitelist', ''),
            expected_cells=config.get('expected_cells', ''),
        log: f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/rad_fake.log"
        group: "alevin_fry"
        shell:
            """
            salmon alevin -l {params.libtype} -i {input.idx}\
                -1 {input.mate1} -2 {input.mate2}\
                --{params.tech} --output {output.outdir} --tgMap {input.tx2gene}\
                -p {threads}\
                --sketch\
                --mrna {input.mt_genes}\
                --rrna {input.rrna_genes}\
                --dumpFeatures\
                > {output.log} 2>&1
            """
else:
    rule rad_se:
        input:
            read = lambda w : [i for i in sample_dict[w.sample]],
            mt_genes = mt_genes,
            rrna_genes = rrna_genes,
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
                f"{idx}/refseq.bin",
                #f"{idx}/reflengths.bin",
                f"{idx}/seq.bin",
                f"{idx}/versionInfo.json"
            ],
            tx2gene = tx2gene,
        output:
            outdir = directory(f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/"),
            map_file = f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/map.rad",
            unmapped = f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/unmapped_bc_count.bin",
            log = f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/rad.log",
            cmd_info = f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/cmd_info.json",
            alevin_log = f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/alevin/alevin.log",
            meta_info = f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/aux_info/meta_info.json",
            libParams = directory(f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/libParams/"),
            salmon_quant = f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/logs/salmon_quant.log",
        #threads: 16
        conda: conda_env_salmon
        #resources:
            #mem_mb=40000
        params:
            libtype = config["libtype_fry"],
            tech = config["tech_fry"],
            whitelist = config.get('barcode_whitelist', ''),
            expected_cells=config.get('expected_cells', ''),
        log: f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/rad_fake.log"
        group: "alevin_fry"
        shell:
            """
            salmon alevin -l {params.libtype} -i {input.idx}\
                -r {input.read}\
                --{params.tech} --output {output.outdir} --tgMap {input.tx2gene}\
                -p {threads}\
                --sketch\
                --mrna {input.mt_genes}\
                --rrna {input.rrna_genes}\
                --dumpFeatures\
                > {output.log} 2>&1
            """

barcode_list_type = config["barcode_list_type"]

print(sns)

if barcode_list_type == "external":

    rule permit_collate_quant_ext_bc:
        input:
            map_file = f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/map.rad",
            unmapped = f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/unmapped_bc_count.bin",
            tx2gene = tx2gene,
            indir = f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/",
            barcode_list = config["barcode_list"]
        output:
            outdir = directory(f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/"),
            quant_json = f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/quant.json",
            quants_mat = f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/alevin/quants_mat.mtx",
            quants_mat_col = f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/alevin/quants_mat_cols.txt",
            quants_mat_row = f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/alevin/quants_mat_rows.txt",
            log = f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/permit_collate_quant.log",
            #all_freq = f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/all_freq.bin",
            collate = f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/collate.json",
            feature_dump = f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/featureDump.txt",
            generate_permit_list = f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/generate_permit_list.json",
            map_collated = f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/map.collated.rad",
            permit_freq = f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/permit_freq.bin",
            permit_map = f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/permit_map.bin",
            unmapped_collated = f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/unmapped_bc_count_collated.bin",
        #threads: 16
        log: f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant_fake.log"
        conda: conda_env_salmon
        group: "alevin_fry"
        params:
            max_records = config["max_records"],
            expected_ori = config["expected_ori"]
        #resources:
            #mem_mb=40000
        shell:
            """
            alevin-fry generate-permit-list -d {params.expected_ori} -i {input.indir} -o {output.outdir} --unfiltered-pl {input.barcode_list} --min-reads 1 > {output.log} 2>&1;
            alevin-fry collate -t {threads} -i {output.outdir} --rad-dir {input.indir} -m {params.max_records} >> {output.log} 2>&1;
            alevin-fry quant -t {threads} -i {output.outdir} -o {output.outdir} --tg-map {input.tx2gene} --resolution cr-like --use-mtx >> {output.log} 2>&1;
            """
        # In the quant step, the duplicates are removed (same UMI).

    rule report_ext_bc:
        input:
            permit_freq = expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/permit_freq.bin", sample = sns),
            feature_dump = expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/featureDump.txt", sample = sns),
            rad_log = expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/rad.log", sample = sns),
            quants_mat = expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/alevin/quants_mat.mtx", sample = sns)
        output:
            report = expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/report.txt", sample = sns)
        #threads: 8
        log: expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/report/report_log.log", sample = sns)
        conda: conda_env_linux
        group: "alevin_fry"
        #resources:
            #mem_mb=48000
        script:
            "../../../scripts/scRNA/2_alignment/report_ext_bc.py"


elif barcode_list_type == "internal":

    rule permit_collate_quant_int_bc:
        input:
            map_file = f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/map.rad",
            unmapped = f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/unmapped_bc_count.bin",
            tx2gene = tx2gene,
            indir = f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/",
            #barcode_list = config["barcode_list"]
        output:
            outdir = directory(f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/"),
            quant_json = f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/quant.json",
            quants_mat = f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/alevin/quants_mat.mtx",
            quants_mat_col = f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/alevin/quants_mat_cols.txt",
            quants_mat_row = f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/alevin/quants_mat_rows.txt",
            log = f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/permit_collate_quant.log",
            all_freq = f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/all_freq.bin",
            collate = f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/collate.json",
            feature_dump = f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/featureDump.txt",
            generate_permit_list = f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/generate_permit_list.json",
            map_collated = f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/map.collated.rad",
            permit_freq = f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/permit_freq.bin",
            permit_map = f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/permit_map.bin",
            unmapped_collated = f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/unmapped_bc_count_collated.bin",
        #threads: 16
        log: f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant_fake.log"
        conda: conda_env_salmon
        group: "alevin_fry"
        params:
            max_records = config["max_records"],
            expected_ori = config["expected_ori"]
        #resources:
            #mem_mb=40000
        shell:
            """
            alevin-fry generate-permit-list -d {params.expected_ori} -k -i {input.indir} -o {output.outdir} > {output.log} 2>&1;
            alevin-fry collate -t {threads} -i {output.outdir} --rad-dir {input.indir} -m {params.max_records} >> {output.log} 2>&1;
            alevin-fry quant -t {threads} -i {output.outdir} -o {output.outdir} --tg-map {input.tx2gene} --resolution cr-like --use-mtx >> {output.log} 2>&1;
            """

    rule report_int_bc:
        input:
            all_freq = expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/all_freq.bin", sample = sns),
            permit_freq = expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/permit_freq.bin", sample = sns),
            feature_dump = expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/featureDump.txt", sample = sns),
            rad_log = expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/rad/rad.log", sample = sns),
            quants_mat = expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/permit_collate_quant/alevin/quants_mat.mtx", sample = sns)
        output:
            report = expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/report.txt", sample = sns)
        #threads: 8
        log: expand(f"{bp}/2_alignment{bm_aln}/{{sample}}/report/report_log.log", sample = sns)
        conda: conda_env_linux
        group: "alevin_fry"
        #resources:
            #mem_mb=16000
        script:
            "../../../scripts/scRNA/2_alignment/report_int_bc.py"
    
        

