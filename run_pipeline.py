import sys
import argparse
from snakemake import snakemake

"""Pipeline calling script"""

def args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', "--container", action="store_true")
    parser.add_argument('-w', "--workdir", required=True, type=str)
    parser.add_argument('-j', "--cores", required=True, type=int)
    return parser.parse_args()


if __name__ == "__main__":
    params = args()
    if params.container:
        print("Running containerized version")

        snakemake("Limix_QTL/rules/limix_qtl_singularity.smk",
                  workdir = params.workdir,
                  cores = params.cores)
    else:
        print("Running local version")

        snakemake("Limix_QTL/rules/limix_qtl.smk",
                  workdir = params.workdir,
                  cores = params.cores)
