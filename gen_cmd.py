""" Generates snakemake commands
"""
import os
import sys
from pathlib import Path
import yaml

from typing import List, Dict, Union

def sanity_dict_clean(snk_dict: Dict[str, Union[str, List[str]]]) -> Dict[str, Union[str, List[str]]]:
    """ Cleans dict since it gets loaded into the global namespace
    """
    unwanted_keys = ['maindir', 'workflow']
    for k in unwanted_keys:
        if snk_dict and k in snk_dict:
            del snk_dict[k]
    return snk_dict


def parse_common_YAML(args, workflow) -> str:
    """Parses the YAML config and the specifics annotations for the pipeline"""
    configfile = f"configs/{args.configfile}"
    with open(configfile, 'r') as f:
        cfg = yaml.load(f, Loader=yaml.FullLoader)
    if args.pipeline == 'imaging':
        pass
    else:
        try:
            genome = f"configs/{args.genome}"
            with open(genome, 'r') as f:
                genome = yaml.load(f, Loader=yaml.FullLoader)
            cfg.update(genome)
        except AttributeError as e:
            print(e)
            sys.exit(1)
    
    # BACKGROUND
    if args.nohup != "hola":
        snakemake_cmd = (
            f"nohup "
        )
    else:
        snakemake_cmd = (
            f""
        )
    
    # GENERAL
    snakemake_cmd += (
        f"snakemake "
        f"--snakefile {workflow['config']} "
    )
    
    # RESOURCES
    snakemake_cmd += (
        "--keep-going "
        "--latency-wait 20 "
        f"--jobs {str(cfg['maxJobs'])} "
        f"--keep-incomplete "
        f"--rerun-incomplete "
    )
    
    # CONDA ENVS / DOCKERS
    if args.docker:
        snakemake_cmd += (
        f"--software-deployment-method conda apptainer "
        f"--conda-prefix {cfg['conda_prefix']} "
        f"--apptainer-prefix {cfg['apptainer_prefix']} " 
        )
    else:
        snakemake_cmd += (
        f"--software-deployment-method conda "
        f"--conda-prefix {cfg['conda_prefix']} " 
        )
    
    # LOCAL / CLOUD
    if args.cloud:
        snakemake_cmd += (
        "--storage-gcs-project "
        "--default-storage-provider gcs "
        "--default-storage-prefix "
        "--executor kubernetes "
        )
    
    # RULE
    if args.rule is not None:
        snakemake_cmd += (
            f"{args.rule} "
        )
    
    # OTHERS
    if args.miss_output:
        snakemake_cmd += (
        "--rerun-triggers code " 
        )
    if args.unlock:
        snakemake_cmd += (
        "--unlock " 
        )
    if args.dryrun:
        snakemake_cmd += (
        "--dryrun "
        )
    if args.unlock:
        # For restarting the pipline if there are errors and snakemake creates a lock
        snakemake_cmd += (
        "--unlock "
        )
    
    # CONFIG + GENOME
    snakemake_cmd += (
        "--config "
    )
    for i in cfg.keys():
        if isinstance(cfg[i], int): 
            snakemake_cmd += f"{i}={cfg[i]} "
        else:
            snakemake_cmd += f"{i}='{cfg[i]}' "
    
    # BACKGROUND
    if args.nohup != "hola":
        snakemake_cmd += (
            f"> {args.nohup}.log 2>&1 & "
        )
         
    return snakemake_cmd