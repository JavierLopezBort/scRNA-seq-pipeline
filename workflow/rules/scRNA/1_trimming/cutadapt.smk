#########################################################3
# PARAMETERS
#########################################################3

paired_end = bool(config['paired_end'])

#########################################################3
# RULES
#########################################################3

if paired_end:
    rule cutadapt_pe:
        input:
            [f"{bp}/FASTQ/{{sample}}{read_preffix}1{read_suffix}{fastq_suffix}",
            f"{bp}/FASTQ/{{sample}}{read_preffix}2{read_suffix}{fastq_suffix}"],
        output:
            fastq1=f"{bp}/1_trimming{bm_tr}/{{sample}}{read_preffix}1{read_suffix}{fastq_suffix}",
            fastq2=f"{bp}/1_trimming{bm_tr}/{{sample}}{read_preffix}2{read_suffix}{fastq_suffix}",
            qc=f"{bp}/1_trimming{bm_tr}/{{sample}}.qc.txt",
        params:
            # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
            adapters = config["adapters_ca"],
            # https://cutadapt.readthedocs.io/en/stable/guide.html#
            extra = config["extra_ca"],
        log:
            f"{bp}/1_trimming{bm_tr}/{{sample}}.log",
        benchmark:
            f"{bp}/1_trimming{bm_tr}/{{sample}}_benchmark.txt",
        threads: 8  # set desired number of threads here
        group: "cutadapt"
        wrapper:
            "v3.3.6/bio/cutadapt/pe"

else:
    rule cutadapt_se:
        input:
            f"{bp}/FASTQ/{{sample}}{fastq_suffix}",
        output:
            fastq=f"{bp}/1_trimming{bm_tr}/{{sample}}{fastq_suffix}",
            qc=f"{bp}/1_trimming{bm_tr}/{{sample}}.qc.txt",
        params:
            adapters = config["adapters_ca"],
            extra = config["extra_ca"],
        log:
            f"{bp}/1_trimming{bm_tr}/{{sample}}.log",
        benchmark:
            f"{bp}/1_trimming{bm_tr}/{{sample}}_benchmark.txt",
        threads: 8  # set desired number of threads here
        group: "cutadapt"
        wrapper:
            "v3.3.6/bio/cutadapt/se"

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
        group: "cutadapt"
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
        group: "cutadapt"
        resources:
            mem_mb=16000
        wrapper:
            "v1.21.2/bio/fastqc"