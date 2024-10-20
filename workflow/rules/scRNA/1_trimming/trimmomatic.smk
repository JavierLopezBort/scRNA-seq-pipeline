#########################################################3
# PARAMETERS
#########################################################3

paired_end = bool(config['paired_end'])

#########################################################3
# TRIMMER PARAMETERS
#########################################################3

adapters = config["adapters_tr"]
adapters_params = config["adapters_params"]
trailing = config["trailing"]
trailing_params = config["trailing_params"]
leading = config["leading"]
leading_params = config["leading_params"]
sliding_window = config["sliding_window"]
sliding_window_params = config["sliding_window_params"]
min_len = config["min_len"]
min_len_params = config["min_len_params"]

def create_trimmer(adapters, adapters_params, trailing, trailing_params, leading, leading_params,
sliding_window, sliding_window_params, min_len, min_len_params):
    trimmer = []
    trimmer.append(f"ILLUMINACLIP:{adapters}:{adapters_params}")
    if trailing:
        trimmer.append(f"TRAILING:{trailing_params}")
    if leading:
        trimmer.append(f"LEADING:{leading_params}")
    if sliding_window:
        trimmer.append(f"SLIDINGWINDOW:{sliding_window_params}")
    if min_len:
        trimmer.append(f"MINLEN:{min_len_params}")
    return trimmer
    
trimmer = create_trimmer(adapters, adapters_params, trailing, trailing_params, leading, leading_params,
sliding_window, sliding_window_params, min_len, min_len_params)

#########################################################3
# RULES
#########################################################3

if paired_end:
    rule trimmomatic_pe:
        input:
            r1=f"{bp}/FASTQ/{{sample}}{read_preffix}1{read_suffix}{fastq_suffix}",
            r2=f"{bp}/FASTQ/{{sample}}{read_preffix}2{read_suffix}{fastq_suffix}",
            adapters=adapters,
        output:
            r1 = f"{bp}/1_trimming{bm_tr}/{{sample}}{read_preffix}1{read_suffix}{fastq_suffix}",
            r2 = f"{bp}/1_trimming{bm_tr}/{{sample}}{read_preffix}2{read_suffix}{fastq_suffix}",
            # reads where trimming entirely removed the mate
            r1_unpaired = f"{bp}/1_trimming{bm_tr}/{{sample}}{read_preffix}1{read_suffix}_unpaired{fastq_suffix}",
            r2_unpaired = f"{bp}/1_trimming{bm_tr}/{{sample}}{read_preffix}2{read_suffix}_unpaired{fastq_suffix}"
        log:
            f"{bp}/1_trimming{bm_tr}/{{sample}}.log"
        benchmark:
            f"{bp}/1_trimming{bm_tr}/{{sample}}_benchmark.txt"
        group:
            "trimmomatic"
        params:
            trimmer = trimmer,
            extra="",
            compression_level = config["compression_level"],
        threads:
            8
        resources:
            mem_mb=32000
        wrapper:
            "v3.3.6/bio/trimmomatic/pe"

else:
    rule trimmomatic_se:
        input:
            reads = f"{bp}/FASTQ/{{sample}}{fastq_suffix}"
        output:
            f"{bp}/1_trimming{bm_tr}/{{sample}}{fastq_suffix}"
        log:
            f"{bp}/1_trimming{bm_tr}/{{sample}}.log"
        benchmark:
            f"{bp}/1_trimming{bm_tr}/{{sample}}_benchmark.txt"
        group:
            "trimmomatic"
        params:
            trimmer,
            extra="",
            compression_level = config["compression_level"],
        threads:
            8
        resources:
            mem_mb=8192
        wrapper:
            "v3.3.6/bio/trimmomatic/se"

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
        group: "trimmomatic"
        resources:
            mem_mb=9000
        wrapper:
            "v3.3.6/bio/fastqc"

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
        group: "trimmomatic"
        resources:
            mem_mb=16000
        wrapper:
            "v3.3.6/bio/fastqc"