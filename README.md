# scRNA-seq-pipeline
scRNA-seq pipeline used for my Master thesis called: "Optimisation of a scRNA-seq bioinformatic pipeline to assess efficiency of meiotic cell induction"

## Prerequisites

The conda environment used for running ivypipe.py:

conda config --add channels bioconda
conda config --add channels defaults
conda config --add channels pytorch
conda config --add channels nvidia
conda config --add channels r

mamba create -n ivysnake
mamba activate ivysnake

mamba install pip
sudo apt update
sudo apt install gcc
pip install snakemake
pip install snakemake-storage-plugin-gcs             # Optional, only for GKE
pip install snakemake-executor-plugin-kubernetes     # Optional, only for GKE

The rest of environments are contained in workflow/envs folder, and are installed directly through conda argument in snakemake rules. Highly recommended to run Snakemake on Linux environment

## Usage

python ivypipe.py -p {pipeline} -c {config file} -g {genome file}

## Details
