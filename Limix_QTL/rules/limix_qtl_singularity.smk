import glob
import os
from pathlib import Path
from subprocess import run
import pandas as pd
import re
from os.path import join

shell.prefix("set -euo pipefail;")

def _multi_arg_start(flag, files):
    flag += " "
    return " ".join(flag + f for f in files)

def _multi_arg_end(flag, files):
    flag = " "+flag
    return " ".join(f + flag for f in files)

def _multi_arg_both_ends(flag1, flag2, files):
    flag1 += " "
    flag2 = " "+flag2
    return " ".join(flag1 + f + flag2 for f in files)

def flatenChunk(chunk):
    return chunk.replace(":", "_").replace("-", "_")

def extendChunk(chunk):
    relChunk = chunk.pop()
    chunkSplitted = relChunk.split("_")
    return chunkSplitted[0]+":"+chunkSplitted[1]+"-"+chunkSplitted[2]

## Configuration file handling
configfile: "../config_files/config.yaml"
OUTPUT_FOLDER = Path(config["outputFolder"])

finalQTLRun = OUTPUT_FOLDER/ 'qtl_results_all.txt'
topQTL = OUTPUT_FOLDER/ 'top_qtl_results_all.txt'

with open(config["chunkFile"],'r') as f:
    chunks = [x.strip() for x in f.readlines()]

qtlOutput = []
for chunk in chunks:
    processedChunk = flatenChunk(chunk)
    processedChunk=expand(OUTPUT_FOLDER+'{chunk}.finished',chunk=processedChunk )
    qtlOutput.append(processedChunk)

## flatten these lists
qtlOutput = [filename for elem in qtlOutput for filename in elem]

rule all:
    input:
        qtlOutput,finalQTLRun,topQTL

rule run_qtl_mapping:
    input:
        af = config["annotationFile"],
        pf = config["phenotypeFile"],
        cf = config["covariateFile"],
        kf = config["kinshipFile"],
        smf = config["sampleMappingFile"]
    output:
        OUTPUT_FOLDER / '{chunk}.finished'
    params:
        gen = config["genotypeFile"],
        od = OUTPUT_FOLDER,
        np = config["parameters"]["numberOfPermutations"],
        maf = config["parameters"]["minorAlleleFrequency"],
        hwe = config["parameters"]["hwequilibrium"],
        w = config["parameters"]["windowSize"],
    run:
        chunkFull = extendChunk({wildcards.chunk})
        shell(
            " singularity exec --bind ~ ~/limix.simg python /limix_qtl/Limix_QTL/run_QTL_analysis.py  "
            " --bgen {params.gen} "
            " -af {input.af} "
            " -pf {input.pf} "
            " -cf {input.cf} "
            " -od {params.od} "
            " -rf {input.kf} "
            " --sample_mapping_file {input.smf} "
            " -gr {chunkFull} "
            " -np {params.np} "
            " -maf {params.maf} "
            " -hwe {params.hwe} "
            " -w {params.w} "
            " -c -gm gaussnorm -bs 500 -rs 0.95 ")
        shell("touch {output}")

rule aggregate_qtl_results:
    input:
        IF = OUTPUT_FOLDER,
        OF = OUTPUT_FOLDER,
        finalFiles = qtlOutput
    output:
        finalQTLRun
    run:
        shell(
            " singularity exec --bind ~ ~/limix.simg python /limix_qtl/Limix_QTL/post-processing_QTL/minimal_postprocess.py "
            " -id {input.IF} "
            " -od {input.OF} "
            " -sfo ")

rule top_feature:
    input:
        IF = OUTPUT_FOLDER,
        OF = OUTPUT_FOLDER,
        finalFile = finalQTLRun
    output:
        topQTL
    run:
        shell(
            " singularity exec --bind ~ ~/limix.simg python /limix_qtl/Limix_QTL/post-processing_QTL/minimal_postprocess.py "
            "-id {input.IF} "
            "-od {input.OF} "
            "-tfb "
            "-sfo ")