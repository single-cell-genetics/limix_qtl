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
    anno_filename = data_path+'Expression/Geuvadis_CEU_Annot_small.txt'
    kinship_filename= data_path+'Genotypes/Geuvadis_chr1_kinship.normalized.txt'
    individual2sample_filename = data_path + 'Geuvadis_CEU_gte.txt'
    min_maf = 0.05
    min_hwe_P=0.001
    min_call_rate =0.95
    blocksize = 50
    output_dir = data_path+'limix_QTL_results_kinship_covs/'
    randomSeed = 73
    chromosome = '1'
    
    ws = 2500000
    
    run_QTL_analysis(pheno_filename,anno_filename,geno_prefix,True,output_dir,ws,
                     min_maf, min_hwe_P, min_call_rate,
                     blocksize,cis_mode=True, seed=randomSeed, n_perm=100, snps_filename=None,feature_filename = None,
                     genetic_range=chromosome,
                     covariates_filename=covariates_filename,
                     kinship_filename=kinship_filename, write_permutations = True,
                     sample_mapping_filename=individual2sample_filename)

if __name__=='__main__':
    test_QTL_analysis()
