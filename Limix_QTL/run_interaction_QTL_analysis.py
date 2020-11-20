from __future__ import division
##External base packages.
import time
import glob
import os
import gc
import pdb
import sys
##External packages.
import pandas as pd
import numpy as np
from sklearn.preprocessing import Imputer
from numpy_sugar.linalg import economic_qs, economic_svd
from limix.stats import effsizes_se, lrt_pvalues
from glimix_core.lmm import LMM
from bgen_reader import read_bgen
#Internal code.
import qtl_output
import qtl_loader_utils
import qtl_parse_args
import qtl_utilities as utils
from qtl_snp_qc import do_snp_qc

#V0.1.4
def run_interaction_QTL_analysis(pheno_filename, anno_filename, geno_prefix, plinkGenotype, output_dir, interaction_term, window_size=250000, min_maf=0.05, min_hwe_P=0.001, min_call_rate=0.95, blocksize=1000,
                     cis_mode=True, skipAutosomeFiltering = False, gaussianize_method=None, minimum_test_samples= 10, seed=np.random.randint(40000), n_perm=0, write_permutations = False, relatedness_score=0.95, feature_variant_covariate_filename = None, snps_filename=None, feature_filename=None, snp_feature_filename=None, genetic_range='all',
                     covariates_filename=None, kinship_filename=None, sample_mapping_filename=None, extended_anno_filename=None, regressCovariatesUpfront = False, regres_snp_from_env=False):
    fill_NaN = Imputer(missing_values=np.nan, strategy='mean', axis=0, copy=False)
    print('Running interaction QTL analysis.')
    lik = 'normal'
    minimumProbabilityStep=0.1
    '''Core function to take input and run QTL tests on a given chromosome.'''
    if relatedness_score is not None:
        relatedness_score = float(relatedness_score)
    [phenotype_df, kinship_df, covariate_df, sample2individual_df,complete_annotation_df, annotation_df, snp_filter_df, snp_feature_filter_df, geneticaly_unique_individuals, minimum_test_samples, feature_list, bim, fam, bed, bgen, chromosome, selectionStart, selectionEnd, feature_variant_covariate_df]=\
    utils.run_QTL_analysis_load_intersect_phenotype_covariates_kinship_sample_mapping(pheno_filename=pheno_filename, anno_filename=anno_filename, geno_prefix=geno_prefix, plinkGenotype=plinkGenotype, cis_mode=cis_mode, skipAutosomeFiltering = skipAutosomeFiltering,
                      minimum_test_samples= minimum_test_samples,  relatedness_score=relatedness_score, snps_filename=snps_filename, feature_filename=feature_filename, snp_feature_filename=snp_feature_filename, selection=genetic_range,
                     covariates_filename=covariates_filename, kinship_filename=kinship_filename, sample_mapping_filename=sample_mapping_filename, extended_anno_filename=extended_anno_filename, feature_variant_covariate_filename=feature_variant_covariate_filename)
    
    mixed = kinship_df is not None
    if (kinship_df is None) or (relatedness_score is None) : 
        geneticaly_unique_individuals = sample2individual_df['iid'].values
    QS = None
    if(feature_list==None or len(feature_list)==0):
        print ('No features to be tested.')
        sys.exit()
    
    #Open output files
    qtl_loader_utils.ensure_dir(output_dir)
    if not selectionStart is None :
        output_writer = qtl_output.hdf5_writer(output_dir+'/qtl_results_{}_{}_{}.h5'.format(chromosome,selectionStart,selectionEnd))
    else :
        output_writer = qtl_output.hdf5_writer(output_dir+'/qtl_results_{}.h5'.format(chromosome))
    if(write_permutations):
        if not selectionStart is None :
            permutation_writer = qtl_output.hdf5_permutations_writer(output_dir+'/perm_results_{}_{}_{}.h5'.format(chromosome,selectionStart,selectionEnd),n_perm)
        else :
            permutation_writer = qtl_output.hdf5_permutations_writer(output_dir+'/perm_results_{}.h5'.format(chromosome),n_perm)

    if(covariate_df is None or not (interaction_term in covariate_df.columns)):
        print ('Interaction term is not found in the covariates')
        print((interaction_term))
        sys.exit()
    
    #Arrays to store indices of snps tested and pass and fail QC SNPs for features without missingness.
    tested_snp_ids = []
    pass_qc_snps_all = []
    fail_qc_snps_all = []
    fail_qc_features = []
    alpha_params = []
    beta_params = []
    n_samples = []
    n_e_samples = []
    na_containing_features=0
    currentFeatureNumber = 0
    snpQcInfoMain = None
    for feature_id in feature_list:
        snpQcInfo = None
        currentFeatureNumber+= 1
        if (len(phenotype_df.loc[feature_id,:]))<minimum_test_samples:
            print("Feature: "+feature_id+" not tested not enough samples do QTL test.")
            fail_qc_features.append(feature_id)
            geneticaly_unique_individuals = tmp_unique_individuals
            continue
        data_written = False
        contains_missing_samples = False
        snpQuery = utils.do_snp_selection(feature_id, complete_annotation_df, bim, cis_mode, window_size, skipAutosomeFiltering)
        snp_cov_df = None
        if(feature_variant_covariate_df is not None):
            if(feature_id in feature_variant_covariate_df['feature'].values):
                covariateSnp = feature_variant_covariate_df['snp_id'].values[feature_variant_covariate_df['feature']==feature_id]
                if(any(i in  bim['snp'].values for i in covariateSnp)):
                    snpQuery_cov = bim.loc[bim['snp'].map(lambda x: x in list(covariateSnp)),:]
                    if(plinkGenotype):
                        snp_cov_df = pd.DataFrame(data=bed[snpQuery_cov['i'].values,:].compute().transpose(),index=fam.index,columns=snpQuery_cov['snp'],)
                    else:
                        ##Here we make some assumptions on the SNPs. They are expected to be ploidy 2!
                        ##Also we don't use a minimal quality to assure a value is present for all samples.
                        print('Warning, during the regression of SNPs we assume ploidy 2.')
                        snp_cov_df_t = pd.DataFrame(columns=fam.index)
                        rowNumber = 0
                        for snpId in snpQuery_cov['i'] :
                            geno = bgen["genotype"][snpId].compute()
                            if(geno["phased"]):
                                snp_df_dosage_t = geno["probs"][:,[0,2]].sum(1).astype(float)
                                snp_df_dosage_t[(np.amax(geno["probs"][:,:2],1)+np.amax(geno["probs"][:,2:4],1))<(1+minimumProbabilityStep)] = float('NaN')
                            else :
                                snp_df_dosage_t = (geno["probs"][:,0]* 2)+geno["probs"][:,1]
                                snp_df_dosage_t[np.amax(geno["probs"][:,:3],1)<((1/3)+minimumProbabilityStep)] = float('NaN')
                            snp_df_dosage_t = pd.Series(snp_df_dosage_t, index= fam.index)
                            snp_df_dosage_t.name = snpId
                            snp_cov_df_t = snp_cov_df_t.append(snp_df_dosage_t)
                            rowNumber = rowNumber +1
                        snp_cov_df_t = snp_cov_df_t.transpose()
        
        if (len(snpQuery) != 0) and (snp_filter_df is not None):
            toSelect = set(snp_filter_df.index).intersection(set(snpQuery['snp']))
            snpQuery = snpQuery.loc[snpQuery['snp'].isin(toSelect)]
        
        if (len(snpQuery) != 0) and (snp_feature_filter_df is not None):
            toSelect = set(np.unique(snp_feature_filter_df['snp_id'].loc[snp_feature_filter_df['feature']==feature_id])).intersection(set(snpQuery['snp']))
            snpQuery = snpQuery.loc[snpQuery['snp'].isin(toSelect)]
        if len(snpQuery) == 0:
            print("Feature: "+feature_id+" not tested. No SNPS passed QC for phenotype.")
            fail_qc_features.append(feature_id)
            continue
        else:
            phenotype_ds = phenotype_df.loc[feature_id]
            contains_missing_samples = any(~np.isfinite(phenotype_ds))
            if(contains_missing_samples):
                print ('Feature: ' + feature_id + ' contains missing data.')
                phenotype_ds.dropna(inplace=True)
                na_containing_features = na_containing_features+1
            '''select indices for relevant individuals in genotype matrix
            These are not unique. NOT to be used to access phenotype/covariates data
            '''
            individual_ids = sample2individual_df.loc[phenotype_ds.index,'iid'].values
            sample2individual_feature= sample2individual_df.loc[phenotype_ds.index]
            
            if(contains_missing_samples):
                tmp_unique_individuals = geneticaly_unique_individuals
                if (kinship_df is not None) and (relatedness_score is not None):
                    geneticaly_unique_individuals = utils.get_unique_genetic_samples(kinship_df.loc[individual_ids,individual_ids], relatedness_score);
                else:
                    geneticaly_unique_individuals = individual_ids
            else:
                #If no missing samples we can use the previous SNP Qc information before actually loading data.
                #This allows for more efficient blocking and retrieving of data
                snpQuery = snpQuery.loc[snpQuery['snp'].map(lambda x: x not in list(map(str, fail_qc_snps_all)))]
            
            if phenotype_ds.empty or len(geneticaly_unique_individuals)<minimum_test_samples :
                print("Feature: "+feature_id+" not tested not enough samples do QTL test.")
                fail_qc_features.append(feature_id)
                if contains_missing_samples:
                    geneticaly_unique_individuals = tmp_unique_individuals
                continue
            elif np.var(phenotype_ds.values) == 0:
                print("Feature: "+feature_id+" has no variance in selected individuals.")
                fail_qc_features.append(feature_id)
                if contains_missing_samples:
                    geneticaly_unique_individuals = tmp_unique_individuals
                continue
            
            print ('For feature: ' +str(currentFeatureNumber)+ '/'+str(len(feature_list))+ ' (' + feature_id + '): ' + str(snpQuery.shape[0]) + ' SNPs need to be tested.\n Please stand by.')
            
            if(n_perm!=0):
                bestPermutationPval = np.ones((n_perm), dtype=np.float)
            
            #Here we need to start preparing the LMM, can use the fam for sample IDS in SNP matrix.
            #test if the covariates, kinship, snp and phenotype are in the same order
            if ((all(kinship_df.loc[individual_ids,individual_ids].index==sample2individual_feature.loc[phenotype_ds.index]['iid']) if kinship_df is not None else True) &\
                 (all(phenotype_ds.index==covariate_df.loc[sample2individual_feature['sample'],:].index)if covariate_df is not None else True)):
                '''
                if all lines are in order put in arrays the correct genotype and phenotype
                x=a if cond1 else b <---> equivalent to if cond1: x=a else x=b;                 better readability of the code
                 '''
                if kinship_df is not None:
                    kinship_mat = kinship_df.loc[individual_ids,individual_ids].values
                    kinship_mat = kinship_mat.astype(float)
                    ##GOWER normalization of Kinship matrix.
                    kinship_mat *= (kinship_mat.shape[0] - 1) / (kinship_mat.trace() - kinship_mat.mean(0).sum())
                    ## This needs to go with the subselection stuff.
                    if(QS is None and not contains_missing_samples):
                        QS = economic_qs(kinship_mat)
                    elif (contains_missing_samples):
                        QS_tmp = QS
                        QS = economic_qs(kinship_mat)
                if kinship_df is None:
                    K = np.eye(len(phenotype_ds.index))
                    if(QS is None and not contains_missing_samples):
                        QS = economic_qs(K)
                    elif (contains_missing_samples):
                        QS_tmp = QS
                        QS = economic_qs(K)
                
                cov_matrix =  covariate_df.loc[sample2individual_feature['sample'],:] if covariate_df is not None else None
                inter = cov_matrix.loc[:,interaction_term]
                cov_matrix = cov_matrix.values
                
                if snp_cov_df is not None:
                    snp_cov_df_tmp = snp_cov_df.loc[individual_ids,:]
                    snp_cov_df_tmp.index=sample2individual_feature['sample']
                    snp_cov_df = pd.DataFrame(fill_NaN.fit_transform(snp_cov_df_tmp))
                    snp_cov_df.index=snp_cov_df_tmp.index
                    snp_cov_df.columns=snp_cov_df_tmp.columns
                    cov_matrix = np.concatenate((cov_matrix,snp_cov_df.values),1)
                    snp_cov_df_tmp = None
                    snp_cov_df = None
                cov_matrix = cov_matrix.astype(float)
            else:
                print ('There is an issue in mapping phenotypes vs covariates and/or kinship')
                sys.exit()
            
            phenotype = utils.force_normal_distribution(phenotype_ds.values,method=gaussianize_method) if gaussianize_method is not None else phenotype_ds.values
            
            phenotype = phenotype.astype(float)
            

            countChunker = 0
            for snpGroup in utils.chunker(snpQuery, blocksize):
                countChunker=countChunker+1
                #print(countChunker)
                #Fix seed at the start of the first chunker so all permutations are based on the same random first split.
                np.random.seed(seed)
                #print(snpGroup)
                snp_idxs = snpGroup['i'].values
                snp_names = snpGroup['snp'].values
                
                tested_snp_ids.extend(snp_names)
                #subset genotype matrix, we cannot subselect at the same time, do in two steps.
                if(plinkGenotype):
                    snp_df = pd.DataFrame(data=bed[snp_idxs,:].compute().transpose(),index=fam.index,columns=snp_names)
                else :
                    snp_df_dosage = pd.DataFrame(np.nan,index=fam.index, columns = snp_names)
                    snp_df = pd.DataFrame(np.nan,index=fam.index, columns = snp_names)
                    rowNumber = 0
                    for snpId in snp_idxs :
                        geno = bgen["genotype"][snpId].compute()
                        if (geno["ploidy"].min()>1 & geno["ploidy"].max()<3) :
                            if(geno["phased"]):
                                snp_df_dosage_t = geno["probs"][:,[0,2]].sum(1).astype(float)
                                snp_df_t = (np.abs(np.argmax(geno["probs"][:,:2], axis=1)-1)+np.abs(np.argmax(geno["probs"][:,2:4], axis=1)-1)).astype(float)
                                naId = (np.amax(geno["probs"][:,:2],1)+np.amax(geno["probs"][:,2:4],1))<(1+minimumProbabilityStep)
                                snp_df_dosage_t[naId] = float('NaN')
                                snp_df_t[naId] = float('NaN')
                            else :
                                snp_df_dosage_t = (geno["probs"][:,0]* 2)+geno["probs"][:,1]
                                snp_df_t = (np.abs(np.argmax(geno["probs"][:,:3], axis=1)-2)).astype(float)
                                naId = np.amax(geno["probs"][:,:3],1)<((1/3)+minimumProbabilityStep)
                                snp_df_dosage_t[naId] = float('NaN')
                                snp_df_t[naId] = float('NaN')
                            snp_df_dosage.loc[:,snp_names[rowNumber]] = snp_df_dosage_t
                            snp_df.loc[:,snp_names[rowNumber]] = snp_df_t
                        rowNumber = rowNumber +1
                    snp_df_dosage = snp_df_dosage.loc[individual_ids,:]
                    
                snp_df = snp_df.loc[individual_ids,:]
                
                snp_df = snp_df.loc[:,np.unique(snp_df.columns)[np.unique(snp_df.columns,return_counts=1)[1]==1]]
                #SNP QC.
                if not contains_missing_samples:
                    #remove SNPs from snp_df if they have previously failed QC
                    snp_df = snp_df.loc[:,snp_df.columns[~snp_df.columns.isin(fail_qc_snps_all)]]
                    if snp_df.shape[1] == 0:
                        continue
                    snps_to_test_df = snp_df.loc[:,snp_df.columns[~snp_df.columns.isin(pass_qc_snps_all)]]
                    if snps_to_test_df.shape[1] > 0:
                        #Only do QC on relevant SNPs. join pre-QCed list and new QCed list.
                        if kinship_df is not None:
                            passed_snp_names,failed_snp_names,call_rate,maf,hweP = do_snp_qc(snps_to_test_df.iloc[np.unique(snps_to_test_df.index,return_index=1)[1]].loc[geneticaly_unique_individuals,:], min_call_rate, min_maf, min_hwe_P)
                        else:
                            passed_snp_names,failed_snp_names,call_rate,maf,hweP = do_snp_qc(snps_to_test_df, min_call_rate, min_maf, min_hwe_P)
                        snps_to_test_df = None
                        #append snp_names and failed_snp_names
                        pass_qc_snps_all.extend(passed_snp_names)
                        fail_qc_snps_all.extend(failed_snp_names)
                    snp_df = snp_df.loc[:,snp_df.columns[snp_df.columns.isin(pass_qc_snps_all)]]
                else:
                    #Do snp QC for relevant section.
                    #Get relevant slice from: phenotype_ds
                    if kinship_df is not None:
                        passed_snp_names,failed_snp_names,call_rate,maf,hweP = do_snp_qc(snp_df.iloc[np.unique(snp_df.index,return_index=1)[1]].loc[geneticaly_unique_individuals,:], min_call_rate, min_maf, min_hwe_P) 
                    else:
                        passed_snp_names,failed_snp_names,call_rate,maf,hweP = do_snp_qc(snp_df, min_call_rate, min_maf, min_hwe_P)
                    snp_df = snp_df.loc[:,snp_df.columns[snp_df.columns.isin(passed_snp_names)]]
                snpQcInfo_t = None
                if call_rate is not None:
                    snpQcInfo_t = call_rate
                    if maf is not None:
                        snpQcInfo_t = pd.concat([snpQcInfo_t,maf.reindex(snpQcInfo_t.index)],axis=1)
                        if hweP is not None:
                            snpQcInfo_t = pd.concat([snpQcInfo_t,hweP.reindex(snpQcInfo_t.index)],axis=1)
                call_rate = None
                maf = None
                hweP = None
                if snpQcInfo is None and snpQcInfo_t is not None:
                    snpQcInfo = snpQcInfo_t
                elif snpQcInfo_t is not None:
                    snpQcInfo = pd.concat([snpQcInfo, snpQcInfo_t], axis=0, sort = False)
                #First process SNPQc than check if we can continue.
                if len(snp_df.columns) == 0:
                    continue
                elif (not plinkGenotype):
                    snp_df_dosage= snp_df_dosage.loc[:,np.unique(snp_df.columns)]
                
                #We could make use of relatedness when imputing.  And impute only based on genetically unique individuals.
                snp_df = pd.DataFrame(fill_NaN.fit_transform(snp_df),index=snp_df.index,columns=snp_df.columns)
                if (not plinkGenotype):
                    snp_df_dosage = pd.DataFrame(fill_NaN.fit_transform(snp_df_dosage),index=snp_df_dosage.index,columns=snp_df_dosage.columns)
                ##No more snp_matrix_DF > snp_df
#                test if the covariates, kinship, snp and phenotype are in the same order
                if (len(snp_df.index) != len(sample2individual_feature.loc[phenotype_ds.index]['iid']) or not all(snp_df.index==sample2individual_feature.loc[phenotype_ds.index]['iid'])):
                    print ('There is an issue in mapping phenotypes and genotypes')
                    sys.exit()
                #print(snp_df)
                #pdb.set_trace()
                for snp_selection in range(snp_df.shape[1]):
                    #print(snp_selection)
                    #pdb.set_trace()
                    snpForTest = snp_df.loc[:,snp_df.columns[snp_selection]].copy(deep=True)
                    if (not plinkGenotype):
                        snpForTest = snp_df_dosage.loc[:,snp_df_dosage.columns[snp_selection]].copy(deep=True)
                    cov_matrix_snp = np.column_stack((cov_matrix, snpForTest))
                    #Add snp to covariate matrix.
                    lmm = LMM(phenotype, cov_matrix_snp, QS)
                    if not mixed:
                        lmm.delta = 1
                        lmm.fix('delta')
                    #Prepare null model.
                    lmm.fit(verbose=False)
                    if regressCovariatesUpfront:
                        phenotype_corrected = phenotype-cov_matrix_snp[:,1:].dot(lmm.beta[1:])
                        cov_matrix_corrected = cov_matrix_snp[:,0]
                        lmm = LMM(phenotype_corrected, cov_matrix_corrected, QS)
                        lmm.fit(verbose=False)
                    #fit new null model.
                    null_lml = lmm.lml()
                    flmm = lmm.get_fast_scanner()
                    #pdb.set_trace()
                    if(regres_snp_from_env):
                        inter = utils.regressOut(inter,np.concatenate(([snpForTest], [np.ones_like(inter.values)]),axis=0).T)
                    G = np.atleast_2d((snpForTest.values * inter.values).T).T
                    G = G.astype(float)
                    G_index = snp_df.columns[snp_selection]
                    
                    alt_lmls, effsizes = flmm.fast_scan(G, verbose=False)
                    var_pvalues = lrt_pvalues(null_lml, alt_lmls)
                    var_effsizes_se = effsizes_se(effsizes, var_pvalues)
                    #pdb.set_trace()
                    #add these results to qtl_results
                    temp_df = pd.DataFrame(index = range(1),columns=['feature_id','snp_id','p_value','beta','beta_se','empirical_feature_p_value'])
                    temp_df['snp_id'] = G_index
                    temp_df['feature_id'] = feature_id
                    temp_df['beta'] = np.asarray(effsizes)
                    temp_df['p_value'] = np.asarray(var_pvalues)
                    temp_df['beta_se'] = np.asarray(var_effsizes_se)
                    #insert default dummy value
                    temp_df['empirical_feature_p_value'] = -1.0
                    
                    if(n_perm!=0):
                        snpForTest = snpForTest.to_frame(name=snp_df.columns[snp_selection])
                        pValueBuffer = []
                        totalSnpsToBeTested = (G.shape[1]*n_perm)
                        permutationStepSize = np.floor(n_perm/(totalSnpsToBeTested/blocksize))
                        if(permutationStepSize>n_perm):
                            permutationStepSize=n_perm
                        elif(permutationStepSize<1):
                            permutationStepSize=1
                        #pdb.set_trace()
                        if(write_permutations):
                            perm_df = pd.DataFrame(index = range(1),columns=['snp_id'] + ['permutation_'+str(x) for x in range(n_perm)])
                            perm_df['snp_id'] = G_index
                        for currentNperm in utils.chunker(list(range(1, n_perm+1)), permutationStepSize):
                            if kinship_df is not None:
                                temp = utils.get_shuffeld_genotypes_preserving_kinship(geneticaly_unique_individuals, relatedness_score, snpForTest, kinship_df.loc[individual_ids,individual_ids], len(currentNperm))
                            else :
                                temp = utils.get_shuffeld_genotypes(snpForTest, len(currentNperm))
                            
                            temp = temp.astype(float)
                            for i in range(0,temp.shape[1]):
                                temp[:,i] = temp[:,i] * inter.values
                            #pdb.set_trace()
                            alt_lmls_p, effsizes_p = flmm.fast_scan(temp, verbose=False)
                            var_pvalues_p = lrt_pvalues(null_lml, alt_lmls_p)
                            pValueBuffer.extend(np.asarray(var_pvalues_p))
                        if(not(len(pValueBuffer)==totalSnpsToBeTested)):
                            #print(len(pValueBuffer))
                            #print(pValueBuffer)
                            #print(totalSnpsToBeTested)
                            print('Error in blocking logic for permutations.')
                            sys.exit()
                        perm = 0
                        for relevantOutput in utils.chunker(pValueBuffer,G.shape[1]) :
                            if(write_permutations):
                                perm_df['permutation_'+str(perm)] = relevantOutput
                            if(bestPermutationPval[perm] > min(relevantOutput)):
                                bestPermutationPval[perm] = min(relevantOutput)
                            perm = perm+1
                            #print(relevantOutput)
                            #print('permutation_'+str(perm))
                
                    if not temp_df.empty :
                        data_written = True
                        output_writer.add_result_df(temp_df)
                        if(write_permutations):
                           permutation_writer.add_permutation_results_df(perm_df,feature_id)
            #This we need to change in the written file.
        if(n_perm>1 and data_written):
            #updated_permuted_p_in_hdf5(bestPermutationPval, feature_id);
            alpha_para, beta_para = output_writer.apply_pval_correction(feature_id,bestPermutationPval, cis_mode)
            #np.savetxt(output_dir+"/Permutation.pValues."+feature_id+".txt",bestPermutationPval)
            alpha_params.append(alpha_para)
            beta_params.append(beta_para)
        if not data_written :
            fail_qc_features.append(feature_id)
        else:
            n_samples.append(phenotype_ds.size)
            n_e_samples.append(len(geneticaly_unique_individuals))
        if contains_missing_samples:
            QS = QS_tmp
            geneticaly_unique_individuals = tmp_unique_individuals
            del QS_tmp
            del tmp_unique_individuals
            if snpQcInfo is not None:
                snpQcInfo.index.name = "snp_id"
                snpQcInfo.to_csv(output_dir+'/snp_qc_metrics_naContaining_feature_{}.txt'.format(feature_id),sep='\t')
        else:
            if (snpQcInfo is not None and snpQcInfoMain is not None):
                snpQcInfoMain = pd.concat([snpQcInfoMain, snpQcInfo], axis=0, sort=False)
            elif snpQcInfo is not None :
                snpQcInfoMain = snpQcInfo.copy(deep=True)
        #if snpQcInfo is not None:
            #snpQcInfo2 = snpQcInfo.copy().transpose()
            #snpQcInfo2.to_csv(output_dir+'/snp_qc_metrics_feature_{}.txt'.format(feature_id),sep='\t')
        #print('step 5')
        gc.collect()
    output_writer.close()
    
    if(write_permutations):
        permutation_writer.close()
    fail_qc_features = np.unique(fail_qc_features)
    if((len(feature_list)-len(fail_qc_features))==0):
        time.sleep(15)
        #Safety timer to make sure the file is unlocked.
        print("Trying to remove the h5 file. Nothing has been tested.")
        print(output_dir+'qtl_results_{}_{}_{}.h5'.format(chromosome,selectionStart,selectionEnd))
        if not selectionStart is None :
            os.remove(output_dir+'qtl_results_{}_{}_{}.h5'.format(chromosome,selectionStart,selectionEnd))
        else :
            os.remove(output_dir+'qtl_results_{}.h5'.format(chromosome))
        sys.exit()
    #gather unique indexes of tested SNPs
    tested_snp_ids = list(set(tested_snp_ids))
    #write annotation and snp data to file
    snp_df = pd.DataFrame()
    snp_df['snp_id'] = bim['snp']
    snp_df['chromosome'] = bim['chrom']
    snp_df['position'] = bim['pos']
    snp_df['assessed_allele'] = bim['a1']
    snp_df.index = snp_df['snp_id']
    snp_df = snp_df.drop_duplicates()
    snp_df = snp_df.reindex(tested_snp_ids)
    snp_df = snp_df.drop_duplicates()
    
    if snpQcInfoMain is not None :
        snpQcInfoMain['index']=snpQcInfoMain.index
        snpQcInfoMain =  snpQcInfoMain.drop_duplicates()
        del snpQcInfoMain['index']
        snp_df = pd.concat([snp_df, snpQcInfoMain.reindex(snp_df.index)], axis=1)
        if(snp_df.shape[1]==5):
            snp_df.columns = ['snp_id','chromosome','position','assessed_allele','call_rate']
        elif(snp_df.shape[1]==6):
            snp_df.columns = ['snp_id','chromosome','position','assessed_allele','call_rate','maf']
        else :
            snp_df.columns = ['snp_id','chromosome','position','assessed_allele','call_rate','maf','hwe_p']
    
    feature_list = list(set(feature_list) - set(fail_qc_features))
    annotation_df = annotation_df.reindex(feature_list)
    annotation_df['n_samples'] = n_samples
    annotation_df['n_e_samples'] = n_e_samples

    if(n_perm>1):
        annotation_df['alpha_param'] = alpha_params
        annotation_df['beta_param'] = beta_params
    
    if not selectionStart is None :
        snp_df.to_csv(output_dir+'/snp_metadata_{}_{}_{}.txt'.format(chromosome,selectionStart,selectionEnd),sep='\t',index=False)
        annotation_df.to_csv(output_dir+'/feature_metadata_{}_{}_{}.txt'.format(chromosome,selectionStart,selectionEnd),sep='\t')
    else :
        snp_df.to_csv(output_dir+'/snp_metadata_{}.txt'.format(chromosome),sep='\t',index=False)
        annotation_df.to_csv(output_dir+'/feature_metadata_{}.txt'.format(chromosome),sep='\t')

        

if __name__=='__main__':
    args = qtl_parse_args.get_interaction_args()
    plink  = args.plink
    bgen = args.bgen
    anno_file = args.annotation_file
    extended_anno_file = args.extended_annotation_file
    pheno_file = args.phenotype_file
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
    cis = args.cis
    trans = args.trans
    interaction_term = args.interaction_term
    write_permutations = args.write_permutations
    includeAllChromsomes = args.no_chromosome_filter
    regress_out_snp_from_env = args.regress_snp_interaction
    regressBefore = args.regress_covariates

    if ((plink is None) and (bgen is None)):
        raise ValueError("No genotypes provided. Either specify a path to a binary plink genotype file or a bgen file.")
    if ((plink is not None) and (bgen is not None)):
        raise ValueError("Only one genotype file can be provided at once, not both plink and bgen")

    if (bgen is not None) :
        plinkGenotype=False
        geno_prefix = bgen
    else:
        plinkGenotype=True
        geno_prefix = plink

    if (cis and trans):
        raise ValueError("cis and trans cannot be specified simultaneously")
    elif (not cis and not trans):
        raise ValueError("At least one run mode (-c / -t) is needed.")
    if (random_seed is None):
        random_seed = np.random.randint(40000)

    if(n_perm==0 and write_permutations):
        write_permutations=False
    if(n_perm==1):
        print("Warning: With only 1 permutation P-value correction is not performed.")
    if(n_perm<50):
        print("Warning: With less than 50 permutations P-values correction is not very accurate.")

    run_interaction_QTL_analysis(pheno_file, anno_file,geno_prefix, plinkGenotype, output_dir, interaction_term, int(window_size),
                     min_maf=float(min_maf), min_hwe_P=float(min_hwe_P), min_call_rate=float(min_call_rate), blocksize=int(block_size),
                     cis_mode=cis, skipAutosomeFiltering= includeAllChromsomes, gaussianize_method = gaussianize, minimum_test_samples= int(minimum_test_samples), seed=int(random_seed), 
                     n_perm=int(n_perm), write_permutations = write_permutations, relatedness_score=relatedness_score, feature_variant_covariate_filename = feature_variant_covariate_filename,
                     snps_filename=snps_filename, feature_filename=feature_filename, snp_feature_filename=snp_feature_filename, genetic_range=genetic_range, covariates_filename=covariates_file,
                     kinship_filename=kinship_file, sample_mapping_filename=samplemap_file, extended_anno_filename=extended_anno_file, regres_snp_from_env = regress_out_snp_from_env, regressCovariatesUpfront = regressBefore)