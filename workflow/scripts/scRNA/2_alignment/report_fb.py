import os
import pandas as pd
import re
import scipy.io
import scanpy as sc
import pandas as pd
import vpolo.alevin.parser as vp
import anndata as an

for log, quants_mat_folder, report in zip(snakemake.input.log, snakemake.input.quants_mat_folder, snakemake.output.report): 

    # Read files

    log = open(log, 'r')
    quants_mat_folder = str(quants_mat_folder)
    quants_mat = vp.read_quants_bin(quants_mat_folder)
    report = open(report, 'w')

    # Log
    report.write(f"LOG\n")
    for line in log:
        if re.search(r"Found \d+ transcripts", line):
            line = re.search(r"\[alevinLog\] \[info\] (.+)",line).group(1)
            report.write(f"{line}\n")
        if re.search(r"Filled with \d+ txp", line):
            line = re.search(r"\[alevinLog\] \[info\] (.+)",line).group(1)
            report.write(f"{line}\n")
        if re.search(r"Found .+ transcripts to gene mappings", line):
            line = re.search(r"\[alevinLog\] \[info\] (.+)",line).group(1)
            report.write(f"{line}\n")
        if re.search(r"# Barcodes Used:", line):
            line = re.search(r"\[alevinLog\] \[info\] (.+)",line).group(1)
            report.write(f"{line}\n")
        if re.search(r"Forcing to use", line):
            line = re.search(r"\[alevinLog\] \[info\] (.+)",line).group(1)
            report.write(f"{line}\n")
        if re.search(r"Throwing \d+ barcodes", line):
            line = re.search(r"\[alevinLog\] \[info\] (.+)",line).group(1)
            report.write(f"{line}\n")
        if re.search(r"low confidence\) barcodes", line):
            line = re.search(r"\[alevinLog\] \[info\] (.+)",line).group(1)
            report.write(f"{line}\n")
        if re.search(r"Total .+ reads will be thrown", line):
            line = re.search(r"\[alevinLog\] \[info\] (.+)",line).group(1)
            report.write(f"{line}\n")
        if re.search(r"Total \d+ CB got", line):
            line = re.search(r"\[alevinLog\] \[info\] (.+)",line).group(1)
            report.write(f"{line}\n")
        if re.search(r"Total Unique barcodes found", line):
            line = re.search(r"\[alevinLog\] \[info\] (.+)",line).group(1)
            report.write(f"{line}\n")
        if re.search(r"Index contained \d+ targets", line):
            line = re.search(r"\[jointLog\] \[info\] (.+)",line).group(1)
            report.write(f"{line}\n")
        if re.search(r"Number of decoys :", line):
            line = re.search(r"\[jointLog\] \[info\] (.+)",line).group(1)
            report.write(f"{line}\n")
        if re.search(r"Computed \d+ rich equivalence classes", line):
            line = re.search(r"\[jointLog\] \[info\] (.+)",line).group(1)
            report.write(f"{line}\n")
        if re.search(r"Counted \d+ total reads", line):
            line = re.search(r"\[jointLog\] \[info\] (.+)",line).group(1)
            report.write(f"{line}\n")
        if re.search(r"Number of fragments discarded because", line):
            line = re.search(r"\[jointLog\] \[info\] (.+)",line).group(1)
            report.write(f"{line}\n")
        if re.search(r"Found \d+ reads with", line):
            line = re.search(r"\[jointLog\] \[warning\] (.+)",line).group(1)
            report.write(f"{line}\n")
        if re.search(r"Mapping rate =", line):
            line = re.search(r"\[jointLog\] \[info\] (.+)",line).group(1)
            report.write(f"{line}\n")
        if re.search(r"Total .+ UMI after deduplicating.", line):
            line = re.search(r"\[alevinLog\] \[info\] (.+)",line).group(1)
            report.write(f"{line}\n")
        if re.search(r"Total \d+ BiDirected Edges.", line):
            line = re.search(r"\[alevinLog\] \[info\] (.+)",line).group(1)
            report.write(f"{line}\n")
        if re.search(r"Total \d+ UniDirected Edges.", line):
            line = re.search(r"\[alevinLog\] \[info\] (.+)",line).group(1)
            report.write(f"{line}\n")
        if re.search(r"Skipped \d+ barcodes due", line):
            line = re.search(r"\[alevinLog\] \[warning\] (.+)",line).group(1)
            report.write(f"{line}\n")
            
    # Quants mat
    report.write(f"\n")
    report.write(f"QUANTS MAT\n")

    read_count = quants_mat.sum().sum()
    n_cells_quant = quants_mat.shape[0]

    report.write(f"Read count: {read_count}\n")
    report.write(f"N cells: {n_cells_quant}\n")

    # Close files
    log.close()
    report.close()