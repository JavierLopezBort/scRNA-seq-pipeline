import struct
import os
import pandas as pd
import re
import scipy.io

for permit_freq, rad_log, feature_dump, quants_mat, report in zip(snakemake.input.permit_freq, snakemake.input.rad_log, snakemake.input.feature_dump, snakemake.input.quants_mat, snakemake.output.report): 
    
    # Read files

    report = open(report, 'w')
    rad_log = open(rad_log, 'r')
    feature_dump_df = pd.read_csv(feature_dump, sep='\t')
    permit_freq = open(permit_freq, 'rb')
    quants_mat = scipy.io.mmread(quants_mat)

    # Rad log
    report.write(f"RAD LOG\n")
    flag = False
    for line in rad_log:
        if re.search(r"Number uniquely mapped", line):
	        flag = True
        if flag:
            report.write(line)

    # Permit freq
    report.write(f"\n")
    report.write(f"PERMIT FREQ BIN\n")

    header = struct.unpack("<Q", permit_freq.read(8))[0]
    report.write(f"Header: {header}\n")

    version = struct.unpack("<Q", permit_freq.read(8))[0]
    report.write(f"Version: {version}\n")

    n_barcodes = struct.unpack("<Q", permit_freq.read(8))[0]
    report.write(f"N cells / barcodes: {n_barcodes}\n")

    freq_dict = dict()

    for i in range(n_barcodes):
        key = struct.unpack("<Q", permit_freq.read(8))[0]
        value = struct.unpack("<Q", permit_freq.read(8))[0]
        freq_dict[key] = value

    n_reads = sum(freq_dict.values())
    report.write(f"N reads: {n_reads}\n\n")

    # Featuredump
    report.write(f"FEATURE DUMP\n")

    n_cells = feature_dump_df.shape[0]
    CorrectedReads = feature_dump_df['CorrectedReads'].sum()
    MappedReads = feature_dump_df['MappedReads'].sum()
    DeduplicatedReads = feature_dump_df['DeduplicatedReads'].sum()

    report.write(f"N cells: {n_cells}\n")
    report.write(f"Corrected Reads: {CorrectedReads}\n")
    report.write(f"Mapped Reads: {MappedReads}\n")
    report.write(f"Deduplicated Reads: {DeduplicatedReads}\n")

    # Quants mat
    report.write(f"\n")
    report.write(f"QUANTS MAT\n")

    quants_mat = quants_mat.todense()
    quants_mat = pd.DataFrame(quants_mat)
    read_count = quants_mat.sum().sum()
    n_cells_quant = quants_mat.shape[0]

    report.write(f"Read count: {read_count}\n")
    report.write(f"N cells: {n_cells_quant}\n")

    # Close files
    permit_freq.close()
    rad_log.close()
    report.close()