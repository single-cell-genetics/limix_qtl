import os
import glob
import sys
import time
from optparse import OptionParser

import dask.dataframe as dd
import h5py
import pandas as pd
import numpy as np
import scipy as sp
import scipy.stats as scst

from struct_lmm import StructLMM
import limix
import qtl_output
import qtl_loader_utils
import qtl_parse_args
import qtl_utilities as utils
from qtl_snp_qc import do_snp_qc
from sklearn.preprocessing import Imputer

def struct_lmm(geno_df,
                pheno,
                env,
                covs = None,
                rhos = None,
                no_association_test = False,
                no_interaction_test = False):
    """
    Utility function to run StructLMM

    Parameters
    ----------
    geno_df: (`N`, `S`) pandas.DataFrame
        genotype matrix for `N` samples, `S` SNPis
    pheno : (`N`, 1) ndarray
        phenotype vector
    env : (`N`, `K`)
          Environmental matrix (indviduals by number of environments)
    covs : (`N`, L) ndarray
        fixed effect design for covariates `N` samples and `L` covariates.
    rhos : list
        list of ``rho`` values.
        ``rho=0`` correspond to no persistent effect (only GxE);
        ``rho=1`` corresponds to only persitent effect (no GxE);
        By default, ``rho=[0, 0.2, 0.4, 0.6, 0.8, 1.]``
    no_association_test : bool
        if True the association test is not considered.
        The default value is False.
    no_interaction_test : bool
        if True the interaction test is not considered.
        The default value is False.
    
    Returns
    -------
    res : *:class:`pandas.DataFrame`*
        contains pv of joint test (if no_association_test is False), 
        pv of interaction test (if no_interaction_test is False).
    """
    #import pdb; pdb.set_trace()
    if covs is None:
        covs = sp.ones((env.shape[0], 1))
    if rhos is None:
        rhos = [0, .2, .4, .6, .8, 1.]
    if not no_association_test:
        # slmm fit null
        slmm = StructLMM(pheno, env, W=env, rho_list=rhos)
        null = slmm.fit_null(F=covs, verbose=False)
    if not no_interaction_test:
        # slmm int
        slmm_int = StructLMM(pheno, env, W=env, rho_list=[0])
    t0 = time.time()
    #import pdb; pdb.set_trace()
    res = pd.DataFrame(data = [], index = geno_df.columns.values)    
    pv = sp.zeros(geno_df.shape[1])
    pv_int = sp.zeros(geno_df.shape[1])
    
    for snp in xrange(geno_df.shape[1]):
        x = geno_df.values[:, [snp]]
        if not no_association_test:
            # association test
            _p, _rho = slmm.score_2_dof(x) #second arg: optimal rho
            pv[snp] = _p
        if not no_interaction_test:
            # interaction test
            covs1 = sp.hstack((covs, x))
            null = slmm_int.fit_null(F=covs1, verbose=False)
            _p, _rho = slmm_int.score_2_dof(x)
            pv_int[snp] = _p
    # add pvalues to res
    if not no_association_test:
        res['pv'] = pv 
    if not no_interaction_test:
        res['pv_int'] = pv_int
    t = time.time() - t0
    #print '%.2f s elapsed' % t
    return res

def run_structLMM_QTL_analysis(pheno_filename, anno_filename, env_filename, geno_prefix, plinkGenotype, output_dir, window_size = 250000, 
                        min_maf = 0.05, min_hwe_P = 0.001, min_call_rate = 0.95, association_mode = True, cis_mode = True,
                        blocksize = 1000, skipAutosomeFiltering = False, gaussianize_method = None, minimum_test_samples = 10, 
                        seed = np.random.randint(40000), n_perm = 0, write_permutations = False, relatedness_score = 0.95, 
                        feature_variant_covariate_filename = None, snps_filename = None, feature_filename = None, snp_feature_filename = None, 
                        genetic_range = 'all', covariates_filename = None, kinship_filename = None, sample_mapping_filename = None, 
                        extended_anno_filename = None):
    fill_NaN = Imputer(missing_values = np.nan, strategy = 'mean', axis = 0)
    print('Running structLMM QTL analysis.')
    '''Core function to take input and run QTL tests on a given chromosome.'''
        
    [phenotype_df, kinship_df, covariate_df, environment_df, sample2individual_df, complete_annotation_df, annotation_df, snp_filter_df, snp_feature_filter_df, 
    genetically_unique_individuals, minimum_test_samples, feature_list,bim,fam,bed, chromosome, selectionStart, selectionEnd, feature_variant_covariate_df]=\
    utils.run_structLMM_QTL_analysis_load_intersect_phenotype_environments_covariates_kinship_sample_mapping(pheno_filename = pheno_filename, anno_filename = anno_filename, 
                env_filename = env_filename, geno_prefix = geno_prefix, plinkGenotype = plinkGenotype, cis_mode = cis_mode, association_mode = association_mode, skipAutosomeFiltering = skipAutosomeFiltering, 
                minimum_test_samples = minimum_test_samples, relatedness_score = relatedness_score, snps_filename = snps_filename, 
                feature_filename = feature_filename, snp_feature_filename = snp_feature_filename, selection = genetic_range, covariates_filename = covariates_filename, 
                kinship_filename = kinship_filename, sample_mapping_filename = sample_mapping_filename, extended_anno_filename = extended_anno_filename, 
                feature_variant_covariate_filename = feature_variant_covariate_filename)
    #import pdb; pdb.set_trace()
    if(feature_list == None or len(feature_list) == 0):
        print ('No features to be tested.')
        sys.exit()
    
    #Open output files
    qtl_loader_utils.ensure_dir(output_dir)
    if not selectionStart is None :
        output_writer = qtl_output.hdf5_writer(output_dir + 'qtl_results_{}_{}_{}.h5'.format(chromosome, selectionStart, selectionEnd))
    else :
        output_writer = qtl_output.hdf5_writer(output_dir + 'qtl_results_{}.h5'.format(chromosome))
    if(write_permutations):
        if not selectionStart is None :
            permutation_writer = qtl_output.hdf5_permutations_writer(output_dir + 'perm_results_{}_{}_{}.h5'.format(chromosome, selectionStart, selectionEnd), n_perm)
        else :
            permutation_writer = qtl_output.hdf5_permutations_writer(output_dir + 'perm_results_{}.h5'.format(chromosome), n_perm)
    
    #Arrays to store indices of snps tested and pass and fail QC SNPs for features without missingness.
    tested_snp_idxs = []
    pass_qc_snps_all = []
    fail_qc_snps_all = []
    fail_qc_features = []
    alpha_params = []
    beta_params = []
    n_samples = []
    n_e_samples = []
    currentFeatureNumber = 0
    snpQcInfoMain = None
    for feature_id in feature_list:
        snpQcInfo = None
        currentFeatureNumber+= 1
        if (len(phenotype_df.loc[feature_id,:]) < minimum_test_samples):
            print("Feature: "+feature_id+" not tested: not enough samples do QTL test.")
            continue
        data_written = False
        snpQuery = utils.do_snp_selection(feature_id, complete_annotation_df, bim, cis_mode, window_size, skipAutosomeFiltering)
        snp_cov_df = None
        if(feature_variant_covariate_df is not None):
            if(feature_id in feature_variant_covariate_df['feature'].values):
                covariateSnp = feature_variant_covariate_df['snp_id'].values[feature_variant_covariate_df['feature'] == feature_id]
                if(any(i in bim['snp'].values for i in covariateSnp)):
                    snpQuery_cov = bim.loc[bim['snp'].map(lambda x: x in list(covariateSnp)),:]
                    snp_cov_df_t = pd.DataFrame(data = bed[snpQuery_cov['i'].values,:].compute().transpose(), index = fam.index, columns = snpQuery_cov['snp'],)
                    snp_cov_df = pd.DataFrame(fill_NaN.fit_transform(snp_cov_df_t))
                    snp_cov_df.index = snp_cov_df_t.index
                    snp_cov_df.columns = snp_cov_df_t.columns
                    snp_cov_df_t = None

        if (len(snpQuery) != 0) and (snp_filter_df is not None):
            snpQuery = snpQuery.loc[snpQuery['snp'].map(lambda x: x in list(map(str, snp_filter_df.index)))]

        if (len(snpQuery) != 0) and (snp_feature_filter_df is not None):
            snpQuery = snpQuery.loc[snpQuery['snp'].map(lambda x: x in list(snp_feature_filter_df['snp_id'].loc[snp_feature_filter_df['feature'] == feature_id]))]
        
        if len(snpQuery) != 0:
            phenotype_ds = phenotype_df.loc[feature_id]
            contains_missing_samples = any(~np.isfinite(phenotype_ds))

            if(contains_missing_samples):
                print ('Feature: ' + feature_id + ' contains missing data.')
                phenotype_ds.dropna(inplace = True)

            '''select indices for relevant individuals in genotype matrix
            These are not unique. NOT to be used to access phenotype/covariates data
            '''
            individual_ids = sample2individual_df.loc[phenotype_ds.index, 'iid'].values
            sample2individual_feature = sample2individual_df.loc[phenotype_ds.index]
            
            if(contains_missing_samples):
                tmp_unique_individuals = genetically_unique_individuals
                genetically_unique_individuals = utils.get_unique_genetic_samples(kinship_df.loc[individual_ids, individual_ids], relatedness_score);
            
            if phenotype_ds.empty or (len(genetically_unique_individuals) < minimum_test_samples) :
                print("Feature: " + feature_id + " not tested: not enough samples do QTL test.")
                fail_qc_features.append(feature_id)
                if contains_missing_samples:
                    genetically_unique_individuals = tmp_unique_individuals
                continue
            elif np.var(phenotype_ds.values) == 0:
                print("Feature: " + feature_id + " has no variance in selected individuals.")
                fail_qc_features.append(feature_id)
                if contains_missing_samples:
                    genetically_unique_individuals = tmp_unique_individuals
                continue

            #If no missing samples we can use the previous SNP Qc information before actually loading data.
            #This allowes for more efficient blocking and retreaving of data
            if not contains_missing_samples:
               snpQuery = snpQuery.loc[snpQuery['snp'].map(lambda x: x not in list(map(str, fail_qc_snps_all)))]
            #import pdb; pdb.set_trace()
            n_sample_t = (phenotype_ds.size)
            n_e_samples_t = (len(genetically_unique_individuals))
            print ('For feature: ' + str(currentFeatureNumber) + '/' +str(len(feature_list)) + ' (' + feature_id + '): ' + str(snpQuery.shape[0]) + ' SNPs need to be tested.\n Please stand by.')
            if(n_perm != 0):
                bestPermutationPval = np.ones((n_perm), dtype = np.float)
            for snpGroup in utils.chunker(snpQuery, blocksize):
                snp_idxs = snpGroup['i'].values
                snp_names = snpGroup['snp'].values

                tested_snp_idxs.extend(snp_idxs)

                #subset genotype matrix, we cannot subselect at the same time, do in two steps.
                snp_df = pd.DataFrame(data = bed[snp_idxs,:].compute().transpose(), index = fam.index, columns = snp_names)
                snp_df = snp_df.loc[individual_ids,:]
                
                #SNP QC.
                #Now we do more proper QC on non-identical samples.
                #However, we do not use it when checking for missingness.
                #That could be extended but gives alot of overhead.
                if not contains_missing_samples:
                    #remove snps from snp_df if they fail QC
                    snp_df = snp_df.loc[:,snp_df.columns[~snp_df.columns.isin(fail_qc_snps_all)]]
                    if snp_df.shape[1] == 0:
                        continue
                    snps_to_test_df = snp_df.loc[:,snp_df.columns[~snp_df.columns.isin(pass_qc_snps_all)]]
                    if snps_to_test_df.shape[1] > 0:
                        #Only do QC on relevant SNPs. join pre-QCed list and new QCed list.
                        if kinship_df is not None:
                            passed_snp_names, failed_snp_names, call_rate, maf, hweP = do_snp_qc(snps_to_test_df.iloc[np.unique(snps_to_test_df.index, return_index = 1)[1]].loc[genetically_unique_individuals,:], 
                                                                                                min_call_rate, min_maf, min_hwe_P)
                        else:
                            passed_snp_names, failed_snp_names, call_rate, maf, hweP = do_snp_qc(snps_to_test_df, min_call_rate, min_maf, min_hwe_P)
                        snps_to_test_df = None
                        #append snp_names and failed_snp_names
                        pass_qc_snps_all.extend(passed_snp_names)
                        fail_qc_snps_all.extend(failed_snp_names)
                    snp_df = snp_df.loc[:,snp_df.columns[snp_df.columns.isin(pass_qc_snps_all)]]
                else:
                    #Do snp QC for relevant section.
                    if kinship_df is not None:
                        passed_snp_names, failed_snp_names, call_rate, maf, hweP = do_snp_qc(snp_df.iloc[np.unique(snp_df.index, return_index = 1)[1]].loc[genetically_unique_individuals,:], 
                                                                                            min_call_rate, min_maf, min_hwe_P) 
                    else:
                        passed_snp_names, failed_snp_names, call_rate, maf, hweP = do_snp_qc(snp_df, min_call_rate, min_maf, min_hwe_P)
                    snp_df = snp_df.loc[:,snp_df.columns[snp_df.columns.isin(passed_snp_names)]]
                #print('step 0')
                if len(snp_df.columns) == 0:
                    continue
                
                snpQcInfo_t = None
                if call_rate is not None:
                    snpQcInfo_t = call_rate.transpose()
                    if maf is not None:
                        snpQcInfo_t = snpQcInfo_t.merge(maf.transpose(), how = 'outer')
                        if hweP is not None:
                            snpQcInfo_t = snpQcInfo_t.merge(hweP.transpose(), how = 'outer')
                call_rate = None
                maf = None
                hweP = None
                if snpQcInfo is None and snpQcInfo_t is not None:
                    snpQcInfo = snpQcInfo_t
                elif snpQcInfo_t is not None:
                    snpQcInfo = pd.concat([snpQcInfo, snpQcInfo_t], axis = 1)
                
                #We could make use of relatedness when imputing.
                snp_matrix_DF = pd.DataFrame(fill_NaN.fit_transform(snp_df), index = snp_df.index, columns = snp_df.columns)
                snp_df = None

#                test if the environments, covariates, kinship, snp and phenotype are in the same order
                if ((all(snp_matrix_DF.index == kinship_df.loc[individual_ids,individual_ids].index) if kinship_df is not None else True) &\
                     (all(phenotype_ds.index == covariate_df.loc[sample2individual_feature['sample'],:].index)if covariate_df is not None else True)&\
                     all(phenotype_ds.index == environment_df.loc[sample2individual_feature['sample'],:].index)&\
                     all(snp_matrix_DF.index == sample2individual_feature.loc[phenotype_ds.index]['iid'])):
                    '''
                    if all lines are in order put in arrays the correct genotype and phenotype
                    x=a if cond1 else b <---> equivalent to if cond1: x=a else x=b;                 better readability of the code
                    '''
                    kinship_mat = kinship_df.loc[individual_ids,individual_ids].values if kinship_df is not None else None
                    env =  environment_df.loc[sample2individual_feature['sample'],:].values
                    if covariate_df is not None:
                        cov_matrix =  covariate_df.loc[sample2individual_feature['sample'],:].values
                    else:
                        cov_matrix = None
                    if(snp_cov_df is not None and cov_matrix is not None):
                        snp_cov_df_tmp = snp_cov_df.loc[individual_ids,:]
                        snp_cov_df_tmp.index = sample2individual_feature['sample']
                        cov_matrix = np.concatenate((cov_matrix,snp_cov_df_tmp.values), 1)
                    elif snp_cov_df is not None :
                        snp_cov_df_tmp = snp_cov_df.loc[individual_ids,:]
                        snp_cov_df_tmp.index = sample2individual_feature['sample']
                        cov_matrix = snp_cov_df_tmp.values
                        #cov_matrix = np.concatenate((np.ones(snp_cov_df_tmp.shape[0]).reshape(np.ones(snp_cov_df_tmp.shape[0]).shape[0],1),snp_cov_df_tmp.values),1)

                    phenotype = utils.force_normal_distribution(phenotype_ds.values, method = gaussianize_method) if gaussianize_method is not None else phenotype_ds.values
                else:
                    print ('There is an issue in mapping phenotypes and genotypes')
                    sys.exit()
                #For limix 1.1 we need to switch to lm our selfs if there is no K.
#                return[snp_matrix_DF,phenotype, kinship_mat,cov_matrix]
                #import pdb; pdb.set_trace()
                try: 
                    if(association_mode == True):
                        res_structlmm = struct_lmm(snp_matrix_DF, np.atleast_2d(phenotype.T).T, np.atleast_2d(env.T).T, covs = np.asarray(cov_matrix, float), no_interaction_test = True)
                    else: # i.e. interaction mode
                        res_structlmm = struct_lmm(snp_matrix_DF, np.atleast_2d(phenotype.T).T, np.atleast_2d(env.T).T, covs = np.asarray(cov_matrix, float), no_association_test = True)
                        # only one of the two analyses can be performed, so copy pvalues from interaction run as simply 'pv'
                        res_structlmm['pv'] = res_structlmm['pv_int']
                except: 
                    print (feature_id)
                    print ('Struct-LMM failed')
                    sys.exit()
                #import pdb; pdb.set_trace()
                if(n_perm!=0):
                    pValueBuffer = []
                    totalSnpsToBeTested = (snp_matrix_DF.shape[1]*n_perm)
                    permutationStepSize = np.floor(totalSnpsToBeTested/blocksize)
                    if(permutationStepSize == 0):
                        permutationStepSize = 1
                    perm = 0;
                    if(write_permutations):
                        perm_df = pd.DataFrame(index = range(len(snp_matrix_DF.columns)), columns = ['snp_id'] + ['permutation_' + str(x) for x in range(n_perm)])
                        perm_df['snp_id'] = snp_matrix_DF.columns
                    for currentNperm in utils.chunker(list(range(1, n_perm+1)), permutationStepSize):
                        if kinship_df is not None:
                            temp = utils.get_shuffeld_genotypes_preserving_kinship(genetically_unique_individuals, relatedness_score, snp_matrix_DF,kinship_df.loc[individual_ids,individual_ids], len(currentNperm))
                        else :
                            temp = utils.get_shuffeld_genotypes(snp_matrix_DF, len(currentNperm))
                        #import pdb; pdb.set_trace()
                        # make into "geno_" dataframe for running structLMM 
                        temp = pd.DataFrame(data = temp)
                        if(association_mode == True):
                            res_structlmm_perms = struct_lmm(temp, np.atleast_2d(phenotype.T).T, np.atleast_2d(env.T).T, covs = np.asarray(cov_matrix, float), no_interaction_test = True)
                        else:
                            res_structlmm_perms = struct_lmm(temp, np.atleast_2d(phenotype.T).T, np.atleast_2d(env.T).T, covs = np.asarray(cov_matrix, float), no_association_test = True)
                            res_structlmm_perms['pv'] = res_structlmm_perms['pv_int'] 
                        pValueBuffer.extend(np.asarray(res_structlmm_perms['pv']))
                    if(not(len(pValueBuffer) == totalSnpsToBeTested)):
                        #print(len(pValueBuffer))
                        #print(pValueBuffer)
                        #print(totalSnpsToBeTested)
                        print('Error in blocking logic for permutations.')
                        sys.exit()
                    for relevantOutput in utils.chunker(pValueBuffer, snp_matrix_DF.shape[1]) :
                        if(write_permutations):
                            perm_df['permutation_' + str(perm)] = relevantOutput
                        if(bestPermutationPval[perm] > min(relevantOutput)):
                            bestPermutationPval[perm] = min(relevantOutput)
                        perm+=1
                #print('step 3')
                #add these results to qtl_results
                temp_df = pd.DataFrame(index = range(len(snp_matrix_DF.columns)), columns = ['feature_id','snp_id','p_value','beta','beta_se','n_samples','empirical_feature_p_value'])
                temp_df['snp_id'] = snp_matrix_DF.columns
                temp_df['feature_id'] = feature_id
                temp_df['beta'] = 0
                temp_df['p_value'] = np.asarray(res_structlmm['pv'])
                temp_df['beta_se'] = -1
                #insert default dummy value
                temp_df['empirical_feature_p_value'] = -1.0
                if not temp_df.empty :
                    data_written = True
                    output_writer.add_result_df(temp_df)
                    if(write_permutations):
                        permutation_writer.add_permutation_results_df(perm_df, feature_id)
            if contains_missing_samples:
                snpQcInfo = snpQcInfo.transpose()
                snpQcInfo.to_csv(output_dir + '/snp_qc_metrics_naContaining_feature_{}.txt'.format(feature_id), sep = '\t')
                genetically_unique_individuals = tmp_unique_individuals
            else:
                if (snpQcInfo is not None and snpQcInfoMain is not None):
                    cols_to_use = snpQcInfo.columns.difference(snpQcInfoMain.columns)
                    snpQcInfoMain = pd.concat([snpQcInfoMain, snpQcInfo[cols_to_use]], axis = 1)
                elif snpQcInfo is not None :
                    snpQcInfoMain = snpQcInfo.copy(deep = True)
                #print('step 4')
            #This we need to change in the written file.
        if(n_perm > 1 and data_written):
            #updated_permuted_p_in_hdf5(bestPermutationPval, feature_id);
            alpha_para, beta_para = output_writer.apply_pval_correction(feature_id, bestPermutationPval, cis_mode)
            alpha_params.append(alpha_para)
            beta_params.append(beta_para)
        if(data_written) :
            n_samples.append(phenotype_ds.size)
            n_e_samples.append(len(genetically_unique_individuals))
        if not data_written :
            fail_qc_features.append(feature_id)
        #print('step 5')
    output_writer.close()
    if(write_permutations):
        permutation_writer.close()

    #gather unique indexes of tested snps
    tested_snp_idxs = list(set(tested_snp_idxs))
    #write annotation and snp data to file
    snp_df = pd.DataFrame()
    snp_df['snp_id'] = bim['snp']
    snp_df['chromosome'] = bim['chrom']
    snp_df['position'] = bim['pos']
    snp_df['assessed_allele'] = bim['a1']
    snp_df = snp_df.ix[tested_snp_idxs,:]
    snp_df = snp_df.drop_duplicates()
    snp_df.index = snp_df['snp_id']
    if snpQcInfoMain is not None :
        snp_df = snp_df.transpose()
        snpQcInfoMain = snpQcInfoMain.transpose()
        snpQcInfoMain['snp_id'] = snpQcInfoMain.index
        snpQcInfoMain =  snpQcInfoMain.drop_duplicates().transpose()
        #print(snp_df)
        #print(snpQcInfoMain)
        snp_df = snp_df.merge(snpQcInfoMain, how = 'outer') 
        snp_df = snp_df.transpose()
        snp_df.columns = ['snp_id', 'chromosome', 'position', 'assessed_allele', 'call_rate', 'maf', 'hwe_pvalue']
    feature_list = [x for x in feature_list if x not in fail_qc_features]
    
    annotation_df = annotation_df.loc[feature_list,:]
    annotation_df['n_samples'] = n_samples
    annotation_df['n_e_samples'] = n_e_samples
    if(n_perm>1):
        annotation_df['alpha_param'] = alpha_params
        annotation_df['beta_param'] = beta_params
    if not selectionStart is None :
        snp_df.to_csv(output_dir + '/snp_metadata_{}_{}_{}.txt'.format(chromosome, selectionStart, selectionEnd), sep = '\t', index = False)
        annotation_df.to_csv(output_dir + '/feature_metadata_{}_{}_{}.txt'.format(chromosome, selectionStart, selectionEnd), sep = '\t')
    else :
        snp_df.to_csv(output_dir + '/snp_metadata_{}.txt'.format(chromosome), sep = '\t', index = False)
        annotation_df.to_csv(output_dir + '/feature_metadata_{}.txt'.format(chromosome), sep = '\t')


if __name__=='__main__':
    args = qtl_parse_args.get_struct_args()
    plink  = args.plink
    bgen = args.bgen
    anno_file = args.annotation_file
    extended_anno_file = args.extended_annotation_file
    pheno_file = args.phenotype_file
    env_file = args.environment_file #structLMM envs
    output_dir = args.output_directory
    window_size = args.window
    genetic_range = args.genomic_range
    covariates_file = args.covariates_file
    kinship_file = args.kinship_file
    samplemap_file = args.sample_mapping_file
    min_maf = args.minor_allel_frequency
    min_hwe_P = args.hardy_weinberg_equilibrium
    min_call_rate = args.call_rate
    block_size = args.block_size
    n_perm = int(args.number_of_permutations)
    snps_filename = args.variant_filter
    snp_feature_filename = args.feature_variant_filter
    feature_variant_covariate_filename = args.feature_variant_covariate
    random_seed = args.seed
    feature_filename = args.feature_filter
    relatedness_score = args.relatedness_score
    minimum_test_samples = args.minimum_test_samples
    gaussianize = args.gaussianize_method
    # can only perform cis or trans analysis
    cis = args.cis
    trans = args.trans
    # can only perform association or interaction structLMM
    association_mode = args.association_mode
    interaction_mode = args.interaction_mode    

    write_permutations = args.write_permutations
    includeAllChromsomes = args.no_chromosome_filter

    if ((plink is None) and (bgen is None)):
        raise ValueError("No genotypes provided. Either specify a path to a binary plink genotype file or a bgen file.")
    if ((plink is not None) and (bgen is not None)):
        raise ValueError("Only one genotype file can be provided at once, not both plink and bgen")

    if (bgen is not None) :
        plinkGenotype = False
        geno_prefix = bgen
        raise ValueError("Not supported")
    else:
        plinkGenotype = True
        geno_prefix = plink

    if (cis and trans):
        raise ValueError("cis and trans cannot be specified simultaneously")
    elif (not cis and not trans):
        raise ValueError("At least one run mode (-c / -t) is needed.")

    if (association_mode and interaction_mode):
        raise ValueError("association and interaction analyses cannot be specified simultaneously")
    elif(not association_mode and not interaction_mode):
        raise ValueError("At least one run mode (-am / -im) is needed.")   

    if (random_seed is None):
        random_seed = np.random.randint(40000)

    if(n_perm == 0 and write_permutations):
        write_permutations = False

    if(n_perm > 1 and n_perm < 10):
        n_perm = 10
        print("Defaults to 10 permutations, if permutations are only used for calibration please give in 1.")
    run_struct_QTL_analysis(pheno_file, anno_file, env_file, geno_prefix, plinkGenotype, output_dir, int(window_size),
                     min_maf = float(min_maf), min_hwe_P = float(min_hwe_P), min_call_rate = float(min_call_rate), association_mode = association_mode, 
                     cis_mode = cis, blocksize = int(block_size), skipAutosomeFiltering = includeAllChromsomes, gaussianize_method = gaussianize, 
                     seed = int(random_seed), minimum_test_samples = int(minimum_test_samples),  n_perm = int(n_perm), write_permutations = write_permutations, 
                     relatedness_score = float(relatedness_score), feature_variant_covariate_filename = feature_variant_covariate_filename,
                     snps_filename = snps_filename, feature_filename = feature_filename, snp_feature_filename = snp_feature_filename, 
                     genetic_range = genetic_range, covariates_filename = covariates_file, kinship_filename = kinship_file, sample_mapping_filename = samplemap_file, 
                     extended_anno_filename = extended_anno_file)
