#########################################################3
# PARAMETERS
#########################################################3

# CREATE INDEX
transcriptome_gz = f"{bf}{config['transcriptome']}"
gtf_gz = f"{bf}{config['gtf']}"

transcriptome = transcriptome_gz.rstrip(".gz")
gtf_unfixed = gtf_gz.rstrip(".gz")

# USE INDEX
idx = f"{bf}{config['idx_STAR']}"

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
# STARSOLO INDEX
#########################################################3

rule STARsolo_index:
    input:
        fasta = transcriptome,
        gtf = gtf,
    output:
        idx = f"{idx}/",
        idx_files = [
            f"{idx}/Genome",
            f"{idx}/Log.out",
            f"{idx}/SA",
            f"{idx}/SAindex",
            f"{idx}/chrLength.txt",
            f"{idx}/chrName.txt",
            f"{idx}/chrNameLength.txt",
            f"{idx}/chrStart.txt",
            f"{idx}/exonGeTrInfo.tab",
            f"{idx}/exonInfo.tab",
            f"{idx}/geneInfo.tab",
            f"{idx}/genomeParameters.txt",
            f"{idx}/sjdbInfo.txt",
            f"{idx}/sjdbList.fromGTF.out.tab",
            f"{idx}/sjdbList.out.tab",
            f"{idx}/transcriptInfo.tab"
        ],
        tmp = f"{idx}/tmp/",
    message:
        "Indexing genome with STAR"
    resources:
        mem_mb=32000,
    conda: conda_env_star
    threads: 8
    log: f"{idx}/STAR.log"
    params:
        extra='--genomeSAindexNbases 2',
    shell:
        "STAR"
        " --runThreadN {threads}"  # Number of threads
        " --runMode genomeGenerate"  # Indexation mode
        " --genomeFastaFiles {input.fasta}"  # Path to fasta files
        " --sjdbGTFfile {input.gtf}"  # Highly recommended GTF
        " {params.extra}"  # Optional parameters
        " --outTmpDir {output.tmp}"  # Temp dir
        " --genomeDir {output.idx}"  # Path to output
        " > {log} 2>&1"  # Logging