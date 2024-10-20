#########################################################3
# PARAMETERS
#########################################################3

idx = f"{bf}{config['idx_STAR']}"
transcriptome = f"{bf}{config['transcriptome_fa']}"

#########################################################3
# RULES
#########################################################3

rule STARsolo:
    input:
        fq2 = lambda w : [i.format(read='1') for i in sample_dict[w.sample]],
        fq1 = lambda w : [i.format(read='2') for i in sample_dict[w.sample]],
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
            ]
    output:
        tmp = f"{bp}/2_alignment{bm_aln}/{{sample}}/tmp/",
        FileNamePrefix = directory(f"{bp}/2_alignment{bm_aln}/{{sample}}/"),

        log_final = f"{bp}/2_alignment{bm_aln}/{{sample}}/Log.final.out",
        log = f"{bp}/2_alignment{bm_aln}/{{sample}}/Log.out",
        log_progress = f"{bp}/2_alignment{bm_aln}/{{sample}}/Log.progress.out",
        sj = f"{bp}/2_alignment{bm_aln}/{{sample}}/SJ.out.tab",

        barcode_stats = f"{bp}/2_alignment{bm_aln}/{{sample}}/Solo.out/Barcodes.stats",
        feature_stats = f"{bp}/2_alignment{bm_aln}/{{sample}}/Solo.out/Velocyto/Features.stats",
        summary = f"{bp}/2_alignment{bm_aln}/{{sample}}/Solo.out/Velocyto/Summary.csv",

        features = f"{bp}/2_alignment{bm_aln}/{{sample}}/Solo.out/Velocyto/filtered/features.tsv",
        barcodes = f"{bp}/2_alignment{bm_aln}/{{sample}}/Solo.out/Velocyto/filtered/barcodes.tsv",
        splided = f"{bp}/2_alignment{bm_aln}/{{sample}}/Solo.out/Velocyto/filtered/spliced.mtx",
        unsplided = f"{bp}/2_alignment{bm_aln}/{{sample}}/Solo.out/Velocyto/filtered/unspliced.mtx",
        ambiguous = f"{bp}/2_alignment{bm_aln}/{{sample}}/Solo.out/Velocyto/filtered/ambiguous.mtx",
        
        aln = f"{bp}/2_alignment{bm_aln}/{{sample}}/Aligned.sortedByCoord.out.bam",
        # unmapped = [
        #     f"{bp}/2_alignment{bm_aln}/{{sample}}/Unmapped.out.mate1", 
        #     f"{bp}/2_alignment{bm_aln}/{{sample}}/Unmapped.out.mate2"
        #     ],
    log:
        f"{bp}/2_alignment{bm_aln}/{{sample}}/{{sample}}.log",
    params:
        extra = f"--soloType CB_UMI_Simple --soloCBwhitelist None --genomeFastaFiles {transcriptome} --soloUMIlen 12 --soloCBlen 16 --soloUMIstart 17 --soloCBstart 1 --soloBarcodeReadLength 0 --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM --outSAMtype BAM SortedByCoordinate --soloUMIfiltering MultiGeneUMI_CR --soloCellFilter  CellRanger2.2 --soloUMIdedup 1MM_CR --outFilterScoreMin 30 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloFeatures Gene Velocyto"
        # extra="--soloType CB_UMI_Simple --soloCBwhitelist None --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 --outSAMtype BAM SortedByCoordinate --outSAMattributes CR UR CY UY CB UB --soloFeatures Gene Velocyto --soloBarcodeReadLength 0" # tenX_v3 (3' v3) # buffalo is tenX_v3
        # extra="--soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 --outSAMtype BAM SortedByCoordinate --outSAMattributes CR UR CY UY CB UB --soloFeatures Gene Velocyto"                                                  # tenX_v2 (3' v2)
        # extra="--soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloStrand Reverse --soloUMIdedup 1MM_CR --outFilterScoreMin 30 --outSAMtype BAM SortedByCoordinate --outSAMattributes CR UR CY UY CB UB --soloFeatures Gene Velocyto"                                                           # tenX_5p (5')
    threads: 8
    group: "STARsolo"
    conda: conda_env_star
    shell:
        """
        STAR --runThreadN 8\
        --genomeDir {input.idx}\
        --readFilesIn {input.fq1} {input.fq2}\
        --readFilesCommand gunzip -c\
        {params.extra}\
        --outReadsUnmapped Fastx\
        --outTmpDir {output.tmp}\
        --outFileNamePrefix {output.FileNamePrefix}
        --outStd BAM_SortedByCoordinate\
        > {output.aln}\
        {log}
        """