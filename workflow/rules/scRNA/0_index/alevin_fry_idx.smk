#########################################################3
# PARAMETERS
#########################################################3

# CREATE INDEX
transcriptome_gz = f"{bf}{config['transcriptome']}"
gtf_gz = f"{bf}{config['gtf']}"
# Dynamically generate a new index if one doesn't exist for that read_length
read_length = config["read_length"]
splici_flank_trim_len = config["splici_flank_trim_len"]
transcriptome_extra = config.get("transcriptome_extra")

transcriptome = transcriptome_gz.rstrip(".gz")
gtf_unfixed = gtf_gz.rstrip(".gz")

# USE INDEX
alevin_folder = f"{bf}{config["alevin_folder"]}"
transcriptome_splici = f"{bf}{config["transcriptome_splici"]}"
tx2gene_alevin = f"{bf}{config["tx2gene"]}"
gene_id_to_name_alevin = f"{bf}{config["gene_id_to_name"]}"
idx = f"{bf}{config["idx_alevin"]}"

#########################################################3
# UNZIP
#########################################################3

if transcriptome_gz != transcriptome:
    rule unzip_trasncriptome:
        input:
            transcriptome_gz = transcriptome_gz,
        output:
            transcriptome = transcriptome,
        shell:
            "gunzip -c {input.transcriptome_gz} > {output.transcriptome}"

if gtf_gz != gtf_unfixed:
    rule unzip_gtf:
        input:
            gtf_gz = gtf_gz,
        output:
            gtf_unfixed = gtf_unfixed,
        shell:
            "gunzip -c {input.gtf_gz} > {output.gtf_unfixed}"

#########################################################3
# GTF FIXED
#########################################################3

gtf = gtf_unfixed.replace(".gtf", ".fixed.gtf")

# ncbi gtf file is technically correct having ### at the end
# however, pyranges cannot handle this
# :TODO check out using the new RUST interval tree package
# Fixed semicolon inside quotes error, skips ### at end - Yejun

if gtf_unfixed != gtf:
    rule fix_pyranges_ncbi_gtf:
        input:
            gtf_unfixed = gtf_unfixed,
        output:
            gtf = gtf
        shell:
            r"""
            awk '
            BEGIN {{ FS="\t"; OFS="\t" }}
            {{
                # Skip lines starting with ###
                if ($0 ~ /^###/) next

                # Print all fields except the last one (attributes) as-is
                for (i=1; i<NF; i++)
                    printf "%s\t", $i
                
                # Process the attributes field
                attr = $NF
                output = ""
                in_quotes = 0
                for (i=1; i<=length(attr); i++) {{
                    c = substr(attr, i, 1)
                    if (c == "\"") {{
                        in_quotes = !in_quotes
                        output = output c
                    }} else if (c == ";" && in_quotes) {{
                        # Skip semicolons inside quotes
                        continue
                    }} else {{
                        output = output c
                    }}
                }}
                
                print output
            }}
            ' {input.gtf_unfixed} > {output.gtf}
            """

#########################################################3
# SPLICI TRANSCRIPTOME
#########################################################3

if transcriptome_extra is None:

    rule pyroe_makesplicedintronic:
        input:
            fasta = transcriptome,
            gtf = gtf
            # spliced="extra_spliced.fa", # Optional
            # unspliced="extra_unspliced.fa", # Optional
        output:
            fasta = transcriptome_splici,
            t2g = tx2gene_alevin,
            gene_id_to_name = gene_id_to_name_alevin,
        params:
            read_length=read_length,
            flank_trim_length=splici_flank_trim_len,
            #extra="--dedup-seqs"
            extra=""
        threads: 4
        log:
            f"{config["alevin_folder"]}/log"
        wrapper:
            "v3.9.0/bio/pyroe/makesplicedintronic"

else:

    rule pyroe_makesplicedintronic_synthetic:
        input:
            fasta = transcriptome,
            gtf = gtf,
            spliced = f"{bf}{transcriptome_extra}"
            # unspliced="extra_unspliced.fa", # Optional
        output:
            fasta = transcriptome_splici,
            t2g = tx2gene_alevin,
            gene_id_to_name = gene_id_to_name_alevin,
        params:
            read_length=read_length,
            flank_trim_length=splici_flank_trim_len,
            #extra="--dedup-seqs"
            extra=""
        threads: 4
        log:
            f"{config["alevin_folder"]}/log"
        wrapper:
            "v3.9.0/bio/pyroe/makesplicedintronic"

#########################################################3
# SALMON INDEX
#########################################################3

rule salmon_index:
    input:
        sequences = transcriptome_splici
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
            f"{idx}/reflengths.bin",
            f"{idx}/refseq.bin", # only one diff from the previous one
            f"{idx}/seq.bin",
            f"{idx}/versionInfo.json",
        ]
        tmp = f"{idx}/tmp/",
    log:
        f"{idx}/alevin_fry.log",
    threads: 2
    params:
        # optional parameters
        extra="",
    shell:
        "salmon index "
        "--transcripts {input.sequences} "
        "--index {output.idx} "
        "--threads {threads} "
        "--tmpdir {output.tmp} "
        "{extra} "
        " > {log} 2>&1"  # Logging