#########################################################3
# RULES
#########################################################3

if paired_end:
    rule fastqc_raw_pe:
        input:
            f"{bp}/FASTQ/{{sample}}{read_preffix}{{read}}{read_suffix}{fastq_suffix}"
        output:
            html=f"{bp}/FASTQ/fastqc{bm_rw}/{{sample}}{read_preffix}{{read}}{read_suffix}.html",
            zip=f"{bp}/FASTQ/fastqc{bm_rw}/{{sample}}{read_preffix}{{read}}{read_suffix}.zip",
        params: "--verbose"
        log:
            f"{bp}/FASTQ/fastqc{bm_rw}/{{sample}}{read_preffix}{{read}}{read_suffix}.log"
        threads: 1
        resources:
            mem_mb=8000
        wrapper:
            "v3.5.3/bio/fastqc"

else:
    rule fastqc_raw_se:
        input:
            f"{bp}/FASTQ/{{sample}}{fastq_suffix}"
        output:
            html=f"{bp}/FASTQ/fastqc{bm_rw}/{{sample}}.html",
            zip=f"{bp}/FASTQ/fastqc{bm_rw}/{{sample}}.zip",
        params: "--verbose"
        log:
            f"{bp}/FASTQ/fastqc{bm_rw}/{{sample}}_fastqc.log"
        threads: 1
        resources:
            mem_mb=16000
        wrapper:
            "v3.5.3/bio/fastqc"