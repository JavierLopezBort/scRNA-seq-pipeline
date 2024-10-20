# INSTRUCTIONS
"""
1- scRNA_Jav.smk: uncomment the step you want to execute in rule all. Plots are always the final outputs of each step
2- sample_collection.smk: define fastq sample files. There are more instructions inside the file
3- scRNA_Jav.yaml: decide bp (base path for the sample files) and bm (try used for benchmarking)
4- command line: $ ivipipe.py -p scRNA_Jav -c scRNA_Jav.yaml -g GRCh38_cloud.yaml --cloud --keep_files
"""

# BENCHMARK
# INSTRUCTIONS
# bm is a variable added at the end of each step folder. If "None", an empty string is returned and pipeline is runned
# normally. If you want to benchmark or try another parameters in a specific step without deleting the previous results,
# you can add a new try changing "None" by "_1" for example, and this string will be added to the step folder. Please,
# every time you add a new try, use the same bm variable for all downstream steps to keep the benchmarking try and not
# mixing results

# INSTRUCTIONS
"""
This sample construction only works with the following structure, which usually matches most of the samples:
- sample_name
- read preffix (if the sample is paired end)
- read number
- read suffix (sometimes)
- fastq suffix (extension of the fastq file)

# REMEMBER
This sample construction DOES NOT WORK if there is some sample strings AFTER the read number that are different
between different samples (wildcards will not match, new combinations will be created). Therefore, we assume the read
number is the last part of the sample name. However, we can also add a read suffix, IN CASE is fixed for every sample


Fill the read_preffix, read_suffix and fastq_suffix variables. Read preffix contains the preffix
that separates the read number (1 or 2) from the rest of sample name, including _. If the sample
is not paired_end, leave it empty. Read_suffix contains every FIXED string between the read number and the fastq
extension. Fastq suffix refers to the fastq file extension (including compression). It MUST contain the dot

For the sample names, you need to create sample_dict_raw by hand. Basically, each key corresponds to a sample and each
value from the same key are different lanes, sequencing runs, etc from the sample. Decide the name of each sample for each
key. However, the values should match with the name of the files. Each value should contain all the sample string since 
read_preffix, which we assume is the same for every sample. If each sample is stored in ONLY one file, each key will only
contain one value, but the dictionary structure is kept.
"""

# INSTRUCTIONS
"""
Uncomment step you want to execute. Plots are always the final output of each step
"""

