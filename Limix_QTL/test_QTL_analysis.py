#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#from depricated_run_QTL_analysis_limix_1 import run_QTL_analysis
from run_QTL_analysis import run_QTL_analysis
from qtl_utilities import merge_QTL_results
import subprocess
import numpy as np
import pandas as pd
import pytest

def hdf5_results_checking(filename,fun=lambda df: np.mean(df['beta']) ):
    '''For a given hdf5 results file, returns a value derived from the first dataframe
    in the file.'''
    with pd.HDFStore(filename,'r') as h5file:
        feature_id = sorted(h5file.keys())[0]
        df = h5file.select(feature_id)
    return fun(df)

def results_checking(results_checking_dict,error_tolerance=1e-6):
    for f in results_checking_dict.keys():
        check_value = hdf5_results_checking(f)
        print(f,check_value)
        assert(abs(results_checking_dict[f]-check_value)<error_tolerance)


def test_QTL_analysis():
    '''Run a set of test cases'''
    data_path = '../geuvadis_CEU_test_data/'
    covariates_filename = data_path+'Expression/Geuvadis_CEU_YRI_covariates.txt'
    geno_prefix = data_path+'Genotypes/Geuvadis'
    pheno_filename = data_path+'Expression/Geuvadis_CEU_YRI_Expr.txt.gz'
    anno_filename = data_path+'Expression/Geuvadis_CEU_Annot.txt'
    kinship_filename= data_path+'Genotypes/Geuvadis_chr1_kinship.normalized.txt'
    individual2sample_filename = data_path + 'Geuvadis_CEU_gte.txt'
    min_maf = 0.05
    min_hwe_P=0.001
    min_call_rate =0.95
    blocksize = 1000
    output_dir = data_path+'limix_QTL_results_kinship_covs/'
    randomSeed = 73
    genetic_range = '1'
    
    ws = 2500000
    
    run_QTL_analysis(pheno_filename, anno_filename, geno_prefix, True, output_dir, ws,
                     min_maf, min_hwe_P, min_call_rate, blocksize, cis_mode=True,
                     seed=randomSeed, n_perm=100, snps_filename=None, feature_filename = None,
                     genetic_range=genetic_range, covariates_filename=covariates_filename,
                     kinship_filename=kinship_filename, sample_mapping_filename=individual2sample_filename)

    #results_checking_dict = {output_dir+'qtl_results_1.h5':-0.015720008359251764}
    #results_checking(results_checking_dict)
    
    output_dir = data_path+'limix_QTL_results_kinship_covs_cmd_line/'
    subprocess.call('python run_QTL_analysis.py '
                    '--plink {geno_prefix} '
                    '-af {anno_file} '
                    '-pf {pheno_file} '
                    '-od {output_dir} '
                    '--window {ws} '
                    '-gr {genetic_range} '
                    '--covariates_file {covariates_file} '
                    '--kinship_file {kinship_file} '
                    '-smf {samplemap_file} '
                    '-c'
                    .format(geno_prefix=geno_prefix,
                            anno_file=anno_filename,
                            pheno_file=pheno_filename,
                            ws=ws,
                            plinkGenotype = True,
                            output_dir=output_dir,
                            genetic_range=genetic_range,
                            covariates_file=covariates_filename,
                            kinship_file=kinship_filename,
                            samplemap_file=individual2sample_filename
                            ),
                    shell=True)

    #results_checking_dict = {output_dir+'qtl_results_1.h5':-0.015720008359251764}
    #results_checking(results_checking_dict)

    #run again, without specifying chromosome
    subprocess.call('python run_QTL_analysis.py '
                    '--plink {geno_prefix} '
                    '-af {anno_file} '
                    '-pf {pheno_file} '
                    '-od {output_dir} '
                    '--window {ws} '
                    '-cf {covariates_file} '
                    '-kf {kinship_file} '
                    '-smf {samplemap_file} '
                    '-c'
                    .format(geno_prefix=geno_prefix,
                            anno_file=anno_filename,
                            pheno_file=pheno_filename,
                            output_dir=output_dir,
                            ws=ws,
                            covariates_file=covariates_filename,
                            kinship_file=kinship_filename,
                            samplemap_file=individual2sample_filename
                            ),
                    shell=True)

    results_checking_dict = {output_dir+'qtl_results_all.h5':-0.015720008359251764}
    #results_checking(results_checking_dict)


    #run another test case
    
    geno_prefix = data_path+'Genotypes/Geuvadis'
    pheno_filename = data_path+'Expression/Geuvadis_CEU_Expr.txt'
    anno_filename = data_path+'Expression/Geuvadis_CEU_Annot.txt'
    
    output_dir = data_path+'TestOutput/limix_QTL_results/'
        
    ws = 250000
    
    for chromosome in ['1','2']:
        run_QTL_analysis(pheno_filename,anno_filename,geno_prefix,True,output_dir,ws,min_maf, min_hwe_P,min_call_rate,blocksize,cis_mode=True, seed=randomSeed, n_perm=50,chromosome=chromosome)

    results_checking_dict = {output_dir+'qtl_results_1.h5':0.034497,
                        output_dir+'qtl_results_2.h5':0.002150}
    results_checking(results_checking_dict)


if __name__=='__main__':
    test_QTL_analysis()
