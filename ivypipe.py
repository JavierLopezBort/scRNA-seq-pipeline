__description__="""
Generates snakemake command for various bioinformatics pipelines

usage example:
    ivypipe.py -p pipeline -d working-dir -g genome -c config.yaml
"""

import argparse
import sys
import textwrap
import gen_cmd as igc
import os


def make_argparser() -> argparse.ArgumentParser:

    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(__description__),
        add_help=False,
    )
    parser.add_argument(
        "-c",
        "--configfile",
        dest="configfile",
    )
    parser.add_argument(
        "-g",
        "--genome",
        dest="genome",
    )
    parser.add_argument(
        "-p",
        "--pipeline",
        dest="pipeline",
        help='index, rnaseq or scRNA',
    )
    parser.add_argument(
        "-r",
        "--rule",
        dest="rule",
        help='run a specific rule',
        default=None
    )

    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true")
    optional = parser.add_argument_group("Options")
    optional.add_argument(
        "--dryrun", 
        dest="dryrun", 
        action="store_true", 
        help="Dry run of the pipeline"
    )
    optional.add_argument(
        "--unlock",
        dest="unlock",
        action="store_true",
        help="Unlock file lock",
    )
    optional.add_argument(
        "--cloud",
        dest="cloud",
        action="store_true",
        help="Add arguments for the cloud, note make sure that genome config is the cloud path",
    )
    optional.add_argument(
        "--keep_files",
        dest="keep_files",
        action="store_true",
        help="Keep intermediate and already created files even though the job was incompleted",
    )
    optional.add_argument(
        "--miss_output",
        dest="miss_output",
        action="store_true",
        help="Keep intermediate and already created files even though the job was incompleted",
    )
    optional.add_argument(
        "--nohup",
        dest="nohup",
        default="hola",
        help="Keep intermediate and already created files even though the job was incompleted",
    )
    optional.add_argument(
        "--docker",
        dest="docker",
        action="store_true",
        help="Keep intermediate and already created files even though the job was incompleted",
    )
    return parser


def main():
    pipeline_dict = {
        'scRNA-seq': 'workflow/scRNA-seq.smk',
    }
    parser = make_argparser()
    args = parser.parse_args()
    # :TODO add sensible cluster configs
    #cluster_config = "configs/cluster.yaml"
    try:
        smk_file = pipeline_dict[args.pipeline]
    except KeyError:
        print("Pipeline not found")
        sys.exit(1)
    workflow = {"config": smk_file}

    snakemake_cmd = igc.parse_common_YAML(args, workflow)
    
    if args.dryrun:
        snakemake_cmd += "--dryrun "
    if args.unlock:
        # For restarting the pipline if there are errors and snakemake creates a lock
        snakemake_cmd += "--unlock "
    print(snakemake_cmd)
    os.system(snakemake_cmd)


if __name__ == "__main__":
    main()