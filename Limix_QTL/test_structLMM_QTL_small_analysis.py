#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from run_structLMM_QTL_analysis import run_structLMM_QTL_analysis
from qtl_utilities import merge_QTL_results
import subprocess
import numpy as np
import pandas as pd
import pytest
import h5py

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
    env_filename = data_path+'Expression/Geuvadis_CEU_YRI_covariates.txt'
    geno_prefix = data_path+'Genotypes/Geuvadis'
    pheno_filename = data_path+'Expression/Geuvadis_CEU_YRI_Expr.txt.gz'
    anno_filename = data_path+'Expression/Geuvadis_CEU_Annot_small.txt'
    kinship_filename= data_path+'Genotypes/Geuvadis_chr1_kinship.normalized.txt'
    individual2sample_filename = data_path + 'Geuvadis_CEU_gte.txt'
    min_maf = 0.25
    min_hwe_P=0.01
    min_call_rate =0.95
    blocksize = 50
    #output_dir = data_path+'limix_QTL_results_kinship_covs_struct/'
    output_dir = '/nfs/leia/research/stegle/acuomo/structlmm_test_results/'
    randomSeed = 73
    chromosome = '1'
    ws = 25000 
    run_structLMM_QTL_analysis(pheno_filename = pheno_filename, anno_filename = anno_filename, 
                    env_filename = env_filename, geno_prefix = geno_prefix, plinkGenotype = True, 
                    output_dir = output_dir, window_size = ws, min_maf = min_maf, 
                    min_hwe_P = min_hwe_P, min_call_rate = min_call_rate, association_mode = True,
                    cis_mode = True, blocksize = blocksize, 
                    seed = randomSeed,  n_perm = 10, write_permutations = True, 
                    genetic_range = chromosome, covariates_filename = covariates_filename, kinship_filename = kinship_filename, 
                    sample_mapping_filename = individual2sample_filename)


if __name__=='__main__':
    test_QTL_analysis()
