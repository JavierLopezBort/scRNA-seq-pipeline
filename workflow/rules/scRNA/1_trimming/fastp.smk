#########################################################3
# PARAMETERS
#########################################################3

paired_end = bool(config['paired_end'])
trimming = config["trimming"]

#########################################################3
# RULES
#########################################################3

if paired_end == True and trimming == True:
    rule fastp_pe:
        input:
            sample = [
                f"{bp}/FASTQ/{{sample}}{read_preffix}1{read_suffix}{fastq_suffix}", 
                f"{bp}/FASTQ/{{sample}}{read_preffix}2{read_suffix}{fastq_suffix}"
                ]
        output:
            trimmed = [
                f"{bp}/1_trimming{bm_tr}/{{sample}}{read_preffix}1{read_suffix}{fastq_suffix}",
                f"{bp}/1_trimming{bm_tr}/{{sample}}{read_preffix}2{read_suffix}{fastq_suffix}"
                ],
            # Unpaired reads separately
            unpaired1=f"{bp}/1_trimming{bm_tr}/{{sample}}{read_preffix}1{read_suffix}_unpaired{fastq_suffix}",
            unpaired2=f"{bp}/1_trimming{bm_tr}/{{sample}}{read_preffix}2{read_suffix}_unpaired{fastq_suffix}",
            # or in a single file
            # unpaired="trimmed/fastp/{sample}.singletons.fastq",
            #merged="trimmed/pe/{sample}.merged.fastq",
            failed=f"{bp}/1_trimming{bm_tr}/{{sample}}_failed{fastq_suffix}",
            html=f"{bp}/1_trimming{bm_tr}/{{sample}}.html",
            json=f"{bp}/1_trimming{bm_tr}/{{sample}}.json"
        log: f"{bp}/1_trimming{bm_tr}/{{sample}}.log"
        benchmark: f"{bp}/1_trimming{bm_tr}/{{sample}}_benchmark.txt"
        params:
            adapters = config["adapters_fp"],
            extra = config["extra"]
        group: "fastp"
        threads: 10
        wrapper:
            "v1.21.2-1-g4c64b964/bio/fastp"

elif paired_end != True:
    rule fastp_se:
        input:
            sample=[f"{bp}/FASTQ/{{sample}}{fastq_suffix}"]
        output:
            trimmed=f"{bp}/1_trimming{bm_tr}/{{sample}}{fastq_suffix}",
            failed=f"{bp}/1_trimming{bm_tr}/{{sample}}_failed{fastq_suffix}",
            html=f"{bp}/1_trimming{bm_tr}/{{sample}}.html",
            json=f"{bp}/1_trimming{bm_tr}/{{sample}}.json"
        log:
            f"{bp}/1_trimming{bm_tr}/{{sample}}.log"
        benchmark:
            f"{bp}/1_trimming{bm_tr}/{{sample}}_benchmark.txt"
        params:
            adapters = config["adapters_fp"],
            extra = config["extra_fp"]
        threads: 8
        group: "fastp"
        wrapper:
            "v3.3.6/bio/fastp"

elif trimming != True:
    rule fastp_pe_wo_trimming:
        input:
            sample=[f"{bp}/FASTQ/{{sample}}{read_preffix}1{read_suffix}{fastq_suffix}",
                    f"{bp}/FASTQ/{{sample}}{read_preffix}2{read_suffix}{fastq_suffix}"]
        output:
            html=f"{bp}/1_trimming{bm_tr}/{{sample}}.html",
            json=f"{bp}/1_trimming{bm_tr}/{{sample}}.json"
        log:
            f"{bp}/1_trimming{bm_tr}/{{sample}}.log"
        benchmark:
            f"{bp}/1_trimming{bm_tr}/{{sample}}_benchmark.txt"
        params:
            extra = config["extra_fp"]
        threads: 8
        group: "fastp"
        wrapper:
            "v3.3.6/bio/fastp"

if paired_end:
    rule fastqc_pe:
        input:
            f"{bp}/1_trimming{bm_tr}/{{sample}}{read_preffix}{{read}}{read_suffix}{fastq_suffix}"
        output:
            html=f"{bp}/1_trimming{bm_tr}/fastqc/{{sample}}{read_preffix}{{read}}{read_suffix}.html",
            zip=f"{bp}/1_trimming{bm_tr}/fastqc/{{sample}}{read_preffix}{{read}}{read_suffix}.zip",
        #params: "--verbose"
        log:
            f"{bp}/1_trimming{bm_tr}/fastqc/{{sample}}{read_preffix}{{read}}{read_suffix}.log"
        benchmark:
            f"{bp}/1_trimming{bm_tr}/fastqc/{{sample}}{read_preffix}{{read}}{read_suffix}_benchmark.txt"
        threads: 1
        group: "fastp"
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
        #params: "--verbose"
        log:
            f"{bp}/1_trimming{bm_tr}/fastqc/{{sample}}_fastqc.log"
        benchmark:
            f"{bp}/1_trimming{bm_tr}/fastqc/{{sample}}_fastqc_benchmark.txt"
        threads: 1
        group: "fastp"
        resources:
            mem_mb=16000
        wrapper:
            "v1.21.2/bio/fastqc"