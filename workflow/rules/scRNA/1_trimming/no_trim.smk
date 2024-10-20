#########################################################3
# PARAMETERS
#########################################################3

paired_end = bool(config['paired_end'])

#########################################################3
# RULES
#########################################################3

if paired_end:
    #trimmer_string = ""
    rule no_trim_pe:
        input:
            r1=f"{bp}/FASTQ/{{sample}}{read_preffix}1{read_suffix}{fastq_suffix}",
            r2=f"{bp}/FASTQ/{{sample}}{read_preffix}2{read_suffix}{fastq_suffix}"
        output:
            r1 = f"{bp}/1_trimming{bm_tr}/{{sample}}{read_preffix}1{read_suffix}{fastq_suffix}",
            r2 = f"{bp}/1_trimming{bm_tr}/{{sample}}{read_preffix}2{read_suffix}{fastq_suffix}",
        log:
            f"{bp}/1_trimming{bm_tr}/{{sample}}.log"
        benchmark:
            f"{bp}/1_trimming{bm_tr}/{{sample}}_benchmark.txt"
        group:
            "no_trim"
        threads:
            8
        resources:
            mem_mb=32000
        shell:
            """
            cp {input.r1} {output.r1}\

            cp {input.r2} {output.r2}\
            """

else:
    rule no_trim_se:
        input:
            f"{bp}/FASTQ/{{sample}}{fastq_suffix}"
        output:
            f"{bp}/1_trimming{bm_tr}/{{sample}}{fastq_suffix}"
        log:
            f"{bp}/1_trimming{bm_tr}/{{sample}}.log"
        benchmark:
            f"{bp}/1_trimming{bm_tr}/{{sample}}_benchmark.txt"
        group:
            "no_trim"
        threads:
            8
        resources:
            mem_mb=8192
        shell:
            """
            cp {input} {output}
            """

if paired_end:
    rule fastqc_pe:
        input:
            f"{bp}/1_trimming{bm_tr}/{{sample}}{read_preffix}{{read}}{read_suffix}{fastq_suffix}"
        output:
            html=f"{bp}/1_trimming{bm_tr}/fastqc/{{sample}}{read_preffix}{{read}}{read_suffix}.html",
            zip=f"{bp}/1_trimming{bm_tr}/fastqc/{{sample}}{read_preffix}{{read}}{read_suffix}.zip",
        params: "--verbose"
        log:
            f"{bp}/1_trimming{bm_tr}/fastqc/{{sample}}{read_preffix}{{read}}{read_suffix}.log"
        benchmark:
            f"{bp}/1_trimming{bm_tr}/fastqc/{{sample}}{read_preffix}{{read}}{read_suffix}_benchmark.txt"
        threads: 1
        group: "no_trim"
        resources:
            mem_mb=16000
        wrapper:
            "v1.21.2/bio/fastqc"

else:
    rule fastqc_se:
        input:
            f"{bp}/1_trimming{bm_tr}/{{sample}}{fastq_suffix}"
        output:
            html=f"{bp}/1_trimming{bm_tr}/fastqc/{{sample}}.html",
            zip=f"{bp}/1_trimming{bm_tr}/fastqc/{{sample}}.zip",
        params: "--verbose"
        log:
            f"{bp}/1_trimming{bm_tr}/fastqc/{{sample}}_fastqc.log"
        benchmark:
            f"{bp}/1_trimming{bm_tr}/fastqc/{{sample}}_fastqc_benchmark.txt"
        threads: 1
        group: "no_trim"
        resources:
            mem_mb=16000
        wrapper:
            "v1.21.2/bio/fastqc"