# scRNA-seq-pipeline
scRNA-seq pipeline used for my Master thesis called: "Optimisation of a scRNA-seq bioinformatic pipeline to assess efficiency of meiotic cell induction"

## Prerequisites

### Installation

The conda environment used for running ivypipe.py:

$ conda config --add channels bioconda

$ conda config --add channels defaults

$ conda config --add channels pytorch

$ conda config --add channels nvidia

$ conda config --add channels r

$ mamba create -n ivysnake

$ mamba activate ivysnake

$ mamba install pip

$ sudo apt update

$ sudo apt install gcc                                               # To install packages through pip

$ export PATH=/usr/bin/gcc:$PATH

$ pip install snakemake

$ pip install snakemake-storage-plugin-gcs                           # Optional, only for GKE

$ pip install snakemake-executor-plugin-kubernetes                   # Optional, only for GKE

$ curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh     # To install packages through github links

$ export PATH=$HOME/.cargo/bin:$PATH


The rest of environments are contained in workflow/envs folder, and are installed directly through conda argument in snakemake rules.

### Extra installation

To run scib and scvelo environments, the files contained in code_files need to subsitute the original ones for each package. The location of each file is stored in the folder_location.txt files.

### Resources

Highly recommended to run Snakemake on Linux environment because some programs only work on this environment.


Some programs need GPU, such as scVI or scGPT. scGPT specifically needs a GPU with a special architecture.


## Usage

python ivypipe.py -p {pipeline} -c {config file} -g {genome file} -- {additional commands}

## Details

Note that all the params in the config and genome file get spit out.  This is mainly due to the fact that a singular config file cannot be passed. Here there is a brief explanation about the commands used:

Pipeline: massive sequencing pipeline you want to run. For example: scRNA-seq, ATAC-seq, etc. It will run a Snakefile containing the main rules for a specific pipeline. This file is located at ivysnake/workflow

Genome config file: genome reference used for alignment. It contains all the files necessary to build the index from the genome reference and the other genome-related files used for alignment. Each genome config file contains one genome reference. For example: genome_human_splici, genome_mouse_splici, etc. This config file is a .yaml file, it is located at ivysnake/configs and the name is genome_specie_extrainfo.yaml.

Config file: it contains all parameters related to the pipeline. For example: method used for each step (scvi for integration, alevin-fry for alignment, etc), parameters (number of neighbours in clustering, dimension reduced method used, etc). This config file is a .yaml file, it is located at ivysnake/configs and the name is config_pipeline.yaml.

ADDITIONAL COMMANDS (Most important ones)

Cloud: it adds the necessary commands to run the pipeline using Google Kubernetes Engine (GKE)

Docker: it adds the necessary commands to run the pipeline using Docker files instead of Conda environments (you must define the container / singularity directive in a Snakemake rule beforehand) 

MORE

You can check how ivypipe package works looking at ivypipe.py file and gen_cmd.py file.

Pipeline

The Snakefile “scRNA-seq” contains several parts, delimited by comments. 

First, you need to define some parameters:

Run: device you want to use to run the pipeline. To add a new one, go to paths.smk file and define the folder where your samples are located.

Samplex: sample or dataset you want to run with the pipeline. You must also define the read prefix, read suffix and fastq suffix of the FASTQ sample files. To add a new one, go to sample_collection.smk, where you need to configure five parameters manually:

bp: sample location folder

sample_dict_raw: define all the Fastq files that belong to the same sample. The key name is up to you, but it must be consistent because you will use it after. The name of the fastq files must be trimmed from the end until L00* (read prefix, read suffix and fastq suffix are defined before, as we saw).

gender_dict: define the gender of each sample (KEY NAMES MUST BE THE SAME). Can be male, female, mixed or None.

cell_type_dict: define the cell type of each sample (KEY NAMES MUST BE THE SAME). Can be (at this moment, but it can be updated) adult, embryo or induced.

study: define the name of the study or dataset. By default, is the same name as “samplex”.

specie: define the specie of the samples from the study or dataset.

The last five parameters are important because they are also defined in the adata columns (adata.obs) when you construct the count matrix at step 3_count_matrix

Benchmark: benchmark section are some tags that you can add to each step. They are added to the end of the step name. This is done in order to run different methods in the same step, so you can benchmark them. In snakemake, if the name of the step or folder in this case, is the same, the Snakemake cannot run, so you need to delete the method or change the name. There is one for each step (bm must be empty or some name with a / at the end). You can check how it works inside individual rules for each step.

Prerequisites: it contains some smk files with important information but they do not need to be modified or defined as in the previous steps. All of them are located at ivysnake/workflow/rules/scRNA:

paths: main folder location. It defines the variable run

conda_envs: it defines the path of the conda environment yaml files. They are located at ivysnake/workflow/envs. The conda env yaml files define a conda environment that is created automatically by Snakemake when the rules run and it will be located at miniforge3/envs. In conda_envs.smk file the environment variables are defined and they are used in each rule depending on the resources need for that rule in the directive conda.

methods: for each pipeline step, you choose the method you want to use in the config.yaml file and it will run the smk file associated with this method in this method.smk file. As said, each method has an associated smk file where some rules are defined (the first one is usually the main one and the second one the plots). All these smk methods files are located at ivysnake/workflow/rules/scRNA/step (methods from the same step are grouped inside the same folder). Each smk method uses different python scripts, and they are located at ivysnake/workflow/scripts/scRNA/step. Again, python scripts from the same step are grouped in the same folder.

Sample_collection: define the sample and some important features related.

Plots_paths: define the paths variables for plots of every step. If you add a new method, with a corresponding new smk file, you need to define the paths of the plots generated in this smk file. These path variable are they used in rule all of pipeline.smk, to initiate the run if some of the plots are missing.

Rules: 

Rules (methods): each step is defined by a smk file method, that has been defined in methods.smk. You must uncomment the step you want to run.

Rules (all): each step is defined by the final plots you want to generate, and the location folder has been defined in plots_paths.smk. You must uncomment the plots you want to generate (must be the same as rule (methods))

Config file (config_scRNA.yaml)

The config file contains a lot of parameters defined by: “parameter: value”. These parameters are passed to the snakemake command line after “--config” and they are stored as config[“parameter”], so they can be used in the smk files. 

First, you have some command line arguments you must define. They are some resources parameters, and the conda_prefix and apptainer_prefix location folders. These are the folders where the conda environments / docker files or images are stored. By default, they are miniforge3/envs, and docker respectively

Secondly, for each step you must choose the method you want to use (they are limited, so if you want to add a new one you need to create a new smk rule as explained before) and define the value of some parameters associated with each method. Some parameters are exclusive for one method and others are shared (the parameters related to plots are usually shared by several methods).

Genome file

Each genome file is defined by only one reference genome (different species, version, modification, etc). It usually contains the index you want to construct, other parameters related to this such as gcf or tx2gene file, etc. It also contains some reference files, such as mitochondrial and ribosomal genes list. For each reference genome / genome file, you can have several indexes, one for each alignment method (for example alevin-fry or STARsolo), however they are in the same file because they have been created using the same reference genome. They are mostly used in the index and alignment step.

 
