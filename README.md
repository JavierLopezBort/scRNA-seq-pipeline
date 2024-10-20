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

python ivypipe.py -p {pipeline} -c {config file} -g {genome file}

## Details
