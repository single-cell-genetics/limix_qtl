import glob
import os
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

## Variables ##
##Files, these are populated now with the test examples and we use here plink genotypes.
chunkFile = './Chunks.txt' ##Not provided in the repo, instructions can be found on the repo wiki.
genotypeFile = '/limix_qtl/Limix_QTL/test_data/Genotypes/Geuvadis'  ##Genotype without file extension. Please update flag in the runner to reflect --plink or --bgen 
kinshipFile = '/limix_qtl/Limix_QTL/test_data/Genotypes/Geuvadis_chr1_kinship'
annotationFile = '/limix_qtl/Limix_QTL/test_data/Expression/Geuvadis_CEU_Annot.txt'
phenotypeFile = '/limix_qtl/Limix_QTL/test_data/Expression/Geuvadis_CEU_YRI_Expr.txt'
covariateFile = '/limix_qtl/Limix_QTL/test_data/Expression/Geuvadis_CEU_YRI_covariates.txt'
randomEffFile = '' #no second random effect in this example
sampleMappingFile = '/limix_qtl/Limix_QTL/test_data/Expression/Geuvadis_CEU_Annot.txt' 
outputFolder = './OutGeuvadis/'

##Settings
numberOfPermutations = '1000'
minorAlleleFrequency = '0.1'
windowSize = '1000000'
hwequilibrium = '0.000001'
FDR = '0.05'
## End Variables ##

finalQTLRun = outputFolder+'qtl_results_all.txt'
topQTL = outputFolder+'top_qtl_results_all.txt'


with open(chunkFile,'r') as f:
    chunks = [x.strip() for x in f.readlines()]

qtlOutput = []
for chunk in chunks:
    #print(chunk)
    processedChunk = flatenChunk(chunk)
    #print(processedChunk)
    processedChunk=expand(outputFolder+'{chunk}.finished',chunk=processedChunk )
    #print(processedChunk)
    qtlOutput.append(processedChunk)

## flatten these lists
qtlOutput = [filename for elem in qtlOutput for filename in elem]

rule all:
    input:
        qtlOutput,finalQTLRun,topQTL

rule run_qtl_mapping:
    input:
        af = annotationFile,
        pf = phenotypeFile,
        cf = covariateFile,
        kf = kinshipFile,
        smf = sampleMappingFile
    output:
        outputFolder + '{chunk}.finished'
    params:
        gen=genotypeFile,
        od = outputFolder,
        np = numberOfPermutations,
        maf = minorAlleleFrequency,
        hwe = hwequilibrium,
        w = windowSize,
    run:
        chunkFull = extendChunk({wildcards.chunk})
        shell(
            " singularity exec --bind ~ ~/limix.simg python /limix_qtl/Limix_QTL/post-processing_QTL/minimal_postprocess.py  "
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
        IF = outputFolder,
        OF = outputFolder,
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
        IF = outputFolder,
        OF = outputFolder,
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

