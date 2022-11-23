# -*- coding: utf-8 -*-
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
import math
from sklearn.impute import SimpleImputer
from glimix_core.lmm import LMM
import dask.array as da
#Internal code.
import qtl_output
import qtl_loader_utils
import qtl_parse_args
import qtl_utilities as utils
from qtl_snp_qc import do_snp_qc
import time
#V0.2.0

def run_PrsQtl_analysis(pheno_filename, anno_filename, prsFile, output_dir, min_call_rate=0.95, blocksize=1000,
                     skipAutosomeFiltering = False, gaussianize_method=None, minimum_test_samples= 10, seed=np.random.randint(40000), n_perm=0, write_permutations = False, 
                     write_zscore = False, write_feature_top_permutations = False, relatedness_score=None, feature_variant_covariate_filename = None,
                     snps_filename=None, feature_filename=None, snp_feature_filename=None, genetic_range='all', covariates_filename=None, randomeff_filename=None,
                     sample_mapping_filename=None, regressCovariatesUpfront = False, lr_random_effect = False, debugger=False):
    if debugger:
        fun_start = time.time()
    
    ##Ensure output directory
    qtl_loader_utils.ensure_dir(output_dir)
    print('Running GRS QT analysis.')
    #Define start objects (p1)
    tot_time = 0
    idx = 0

    fill_NaN = SimpleImputer(missing_values=np.nan, strategy='mean', copy=False)
                                     
    lik = 'normal'
    minimumProbabilityStep=0.1

    '''Core function to take input and run QTL tests on a given chromosome.'''
    if relatedness_score is not None:
        relatedness_score = float(relatedness_score)
    # Intersect files together to list the amount of samples with enough files
    if debugger:
        print(relatedness_score)
        fun_start = time.time()
    [phenotype_df, kinship_df, randomeff_df, covariate_df, sample2individual_df, annotation_df, snp_filter_df, snp_feature_filter_df, geneticaly_unique_individuals, minimum_test_samples, feature_list, risk_df, chromosome, selectionStart, selectionEnd, feature_variant_covariate_df]=\
    utils.run_PrsQtl_analysis_load_intersect_phenotype_covariates_kinship_sample_mapping(pheno_filename=pheno_filename, anno_filename=anno_filename, prsFile=prsFile, skipAutosomeFiltering = skipAutosomeFiltering,
                      minimum_test_samples= minimum_test_samples,  relatedness_score=relatedness_score, snps_filename=snps_filename, feature_filename=feature_filename, snp_feature_filename=snp_feature_filename, selection=genetic_range,
                     covariates_filename=covariates_filename, randomeff_filename=randomeff_filename, sample_mapping_filename=sample_mapping_filename, feature_variant_covariate_filename=feature_variant_covariate_filename)
                
                             
                                                                      

    mixed = kinship_df is not None
    if (kinship_df is None) or (relatedness_score is None) : 
        geneticaly_unique_individuals = np.unique(sample2individual_df['iid'].values)
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

    #Arrays to store indices of snps tested and pass and fail QC SNPs for features without missingness.
    tested_snp_names = []
    fail_qc_features = []
    alpha_params = []
    beta_params = []
    n_samples = []
    n_e_samples = []
    na_containing_features=0
    currentFeatureNumber = 0
    snpQcInfoMain = None
    random_eff_param = []
    log = {}
    Sigma = {}
    Sigma_qs = {}
    randomeff_mix = False
    
    minRho = 0.1
    #minRho = 0
    maxRho = 0.9
    #maxRho = 1
    stepSize=0.1
    #stepSize=0.05
    
    precission =  (len(str(stepSize))-2)
    rho1 = np.round(np.arange(minRho,np.round(maxRho+stepSize,precission),stepSize),precission)
    
    #Here we check if there is a middle value in rho1 to start from.
    if len(rho1)%2==0 :
        rho2 = [None] * (len(rho1)+1)
        halfValue = int(np.round(len(rho1)/2,0))
        rho2[0:halfValue+1] = rho1[0:halfValue+1]
        rho2[(halfValue+2):(len(rho1)+1)] = rho1[(halfValue+1):(len(rho1)+1)]
        rho2[(halfValue+1)] = ((rho1[halfValue] + rho1[halfValue+1])/2)
        rho1 = np.round(rho2,precission)
    for feature_id in feature_list:
        gc.collect()
        snpQcInfo = None
        currentFeatureNumber+= 1
        if (len(phenotype_df.loc[feature_id,:]))<minimum_test_samples:
            print("Feature: "+feature_id+" not tested not enough samples do QTL test (n="+str(len(phenotype_df.loc[feature_id,:]))+").")
            fail_qc_features.append(feature_id)
            geneticaly_unique_individuals = tmp_unique_individuals
            continue
        data_written = False
        contains_missing_samples = False
        snpQuery = risk_df.index.values
        snp_cov_df = None
        
        if(feature_variant_covariate_df is not None):
            if(feature_id in feature_variant_covariate_df['feature_id'].values):
                covariateSnp = feature_variant_covariate_df['snp_id'].values[feature_variant_covariate_df['feature_id']==feature_id]
                if(any(i in  risk_df.index.values for i in covariateSnp)):
                    snp_cov_df = risk_df.loc[risk_df.index.map(lambda x: x in list(covariateSnp)),:].transpose()
        
        if (len(snpQuery) != 0) and (snp_filter_df is not None):
            snpQuery = list(set(snp_filter_df.index).intersection(set(snpQuery)))
        
        if (len(snpQuery) != 0) and (snp_feature_filter_df is not None):
            snpQuery = list(set(np.unique(snp_feature_filter_df['snp_id'].loc[snp_feature_filter_df['feature_id']==feature_id])).intersection(set(snpQuery)))
        
        if len(snpQuery) == 0:
            print("Feature: "+feature_id+" not tested. No SNPS passed QC for phenotype.")
            fail_qc_features.append(feature_id)
            continue
        else:
            phenotype_ds = phenotype_df.loc[feature_id]
            contains_missing_samples = any(~np.isfinite(phenotype_ds))
            if(contains_missing_samples):
                #import pdb; pdb.set_trace()
                print ('Feature: ' + feature_id + ' contains missing data.')
                phenotype_ds.dropna(inplace=True)
                na_containing_features = na_containing_features+1
            '''select indices for relevant individuals in genotype matrix
            These are not unique. NOT to be used to access phenotype/covariates data
            '''
            individual_ids = sample2individual_df.loc[phenotype_ds.index,'iid'].values
            sample2individual_feature= sample2individual_df.loc[phenotype_ds.index]
            
            if contains_missing_samples:
                tmp_unique_individuals = geneticaly_unique_individuals
                if (kinship_df is not None) and (relatedness_score is not None):
                    geneticaly_unique_individuals = utils.get_unique_genetic_samples(kinship_df.loc[individual_ids,individual_ids], relatedness_score);
                else :
                    geneticaly_unique_individuals = individual_ids
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
            
            print ('For feature: ' +str(currentFeatureNumber)+ '/'+str(len(feature_list))+ ' (' + feature_id + '): ' + str(len(snpQuery)) + ' risk scores will be tested.\n Please stand by.')
            if(n_perm!=0):
                bestPermutationPval = np.ones((n_perm), dtype=float)
            
            #Here we need to start preparing the LMM, can use the fam for sample IDS in SNP matrix.
            #test if the covariates, kinship, snp and phenotype are in the same order            
            if ((all(kinship_df.loc[individual_ids,individual_ids].index==sample2individual_feature.loc[phenotype_ds.index]['iid']) if kinship_df is not None else True) &\
                 (all(phenotype_ds.index==randomeff_df.loc[sample2individual_feature['sample'],:].index) if randomeff_df is not None else True) &\
                 (all(phenotype_ds.index==covariate_df.loc[sample2individual_feature['sample'],:].index) if covariate_df is not None else True)):
                '''
                if all lines are in order put in arrays the correct genotype and phenotype
                x=a if cond1 else b <---> equivalent to if cond1: x=a else x=b;                 better readability of the code
                 '''
                if kinship_df is not None and randomeff_df is None:
                    if( not lr_random_effect):
                        kinship_mat = kinship_df.loc[individual_ids,individual_ids].values
                        kinship_mat = kinship_mat.astype(float)
                        ##GOWER normalization of Kinship matrix.
                        kinship_mat *= (kinship_mat.shape[0] - 1) / (kinship_mat.trace() - kinship_mat.mean(0).sum())
                        ## This needs to go with the subselection stuff.
                        if(QS is None or contains_missing_samples):
                            QS = utils.economic_qs(kinship_mat)
                    elif(lr_random_effect):
                        kinship_mat = kinship_df.loc[individual_ids,:].values
                        kinship_mat = kinship_mat.astype(float)
                        if(QS is None or contains_missing_samples):
                            QS = utils.economic_qs_linear(kinship_mat, return_q1=False)
                # combining the two matrices
                if kinship_df is not None and randomeff_df is not None:
                    #Here we need to match names and make sure that the order is the same and the right samples get mixed.
                    if( not lr_random_effect):
                        randomeff_mix = True
                        if(not Sigma_qs or contains_missing_samples):
                            kinship_mat = kinship_df.loc[individual_ids,individual_ids].values
                            kinship_mat = kinship_mat.astype(float)
                            randomeff_mat = randomeff_df.loc[sample2individual_feature['sample'],sample2individual_feature['sample']].values
                            randomeff_mat = randomeff_mat.astype(float)
                            
                            if (not (kinship_mat.shape[0] == randomeff_mat.shape[0] and  kinship_mat.shape[1] == randomeff_mat.shape[1])):
                                print ('There is an issue in mapping between the second random effect term and the main association information.')
                                sys.exit()
                            
                            for rho in rho1:
                                Sigma[rho] = rho * kinship_mat + (1 - rho) * randomeff_mat
                                ##GOWER normalization of Kinship matrix.
                                Sigma[rho] *= (Sigma[rho].shape[0] - 1) / (Sigma[rho].trace() - Sigma[rho].mean(0).sum())
                                Sigma_qs[rho] = utils.economic_qs(Sigma[rho])
                    elif(lr_random_effect):
                        randomeff_mix = True
                        if(not Sigma_qs or contains_missing_samples):
                            ##We only take the genetically unique indivudals to reflect the kinship here.
                            kinship_mat = kinship_df.loc[individual_ids,geneticaly_unique_individuals].values
                            kinship_mat = kinship_mat.astype(float)
                            randomeff_mat = randomeff_df.loc[sample2individual_feature['sample'],:].values
                            randomeff_mat = randomeff_mat.astype(float)
                            
                            if (kinship_mat.shape[0] != randomeff_mat.shape[0]):
                                print ('There is an issue in mapping between the second random effect term and the main association information.')
                                sys.exit()
                            
                            for rho in rho1:
                                ##Stck the two matrices together [per rho value]
                                Sigma_qs[rho] = utils.economic_qs_linear(np.concatenate([np.sqrt(rho) * kinship_mat] + [np.sqrt(1 - rho) * randomeff_mat], axis=1), return_q1=False)

                    else:
                        print ('There is an issue in mapping between the second random effect term and the main association information.')
                        sys.exit()

                ##This cant happen!
                # if kinship_df is None and randomeff_df is not None: 
                #     randomeff_mat = randomeff_df.loc[individual_ids,individual_ids].values
                #     randomeff_mat = randomeff_mat.astype(float)

                #     ##GOWER normalization of Kinship matrix.
                #     randomeff_mat *= (randomeff_mat.shape[0] - 1) / (randomeff_mat.trace() - randomeff_mat.mean(0).sum())
                #     ## This needs to go with the subselection stuff.
                #     if(QS is None and not contains_missing_samples):
                #         QS = utils.economic_qs(randomeff_mat)
                #     elif (contains_missing_samples):
                #         QS = utils.economic_qs(randomeff_mat)

                # creating a fake QS if none random effect is present or use the read depth one
                if kinship_df is None:
                    if randomeff_df is None:
                        K = np.eye(len(phenotype_ds.index))
                        if(QS is None and not contains_missing_samples):
                            QS = utils.economic_qs(K)
                        elif (contains_missing_samples):
                            QS = utils.economic_qs(K)
                    else:
                        if(QS is None and not contains_missing_samples):
                            QS = utils.economic_qs(randomeff_df)
                        elif (contains_missing_samples):
                            QS = utils.economic_qs(randomeff_df)
                cov_matrix =  covariate_df.loc[sample2individual_feature['sample'],:].values if covariate_df is not None else None
                if covariate_df is None:
                    cov_matrix = np.ones((len(individual_ids), 1))
                #pdb.set_trace()
                if snp_cov_df is not None:
                    snp_cov_df_tmp = snp_cov_df.loc[individual_ids,:]
                    snp_cov_df = pd.DataFrame(fill_NaN.fit_transform(snp_cov_df_tmp))
                    snp_cov_df.index=sample2individual_feature['sample']
                    snp_cov_df.columns=snp_cov_df_tmp.columns
                    cov_matrix = np.concatenate((cov_matrix,snp_cov_df.values),1)
                    snp_cov_df_tmp = None
                    snp_cov_df = None
                cov_matrix = cov_matrix.astype(float)
            else:
                print ('There is an issue in mapping phenotypes vs covariates and/or kinship and/or second random effect term.')
                sys.exit()
            
            phenotype = utils.force_normal_distribution(phenotype_ds.values,method=gaussianize_method) if gaussianize_method is not None else phenotype_ds.values
            
            phenotype = phenotype.astype(float)
            ###########################################################################################################################################################
            
            #Prepare LMM
            ##########################################################################################################################################################
            # Computing Null Model 
            if debugger:
                fun_start = time.time()
            if randomeff_mix:
                #mixingParameters = utils.rhoTest(best=None, phenotype = phenotype, cov_matrix=cov_matrix, Sigma_qs=Sigma_qs, mixed=mixed, lastMove=None, rhoArray = rho1, verbose = True)
                mixingParameters = utils.rhoTestBF(best=None, phenotype = phenotype, cov_matrix=cov_matrix, Sigma_qs=Sigma_qs, mixed=mixed, lastMove=None, rhoArray = rho1, verbose = False)
                
                lmm = mixingParameters["lmm"]
                feature_best_rho = mixingParameters["rho"]
                log[(feature_id)].append(feature_best_rho)
                
                if debugger:
                    if mixingParameters["rho"]!=0:
                        print("Random effect has influence, mixing parameter: "+str(feature_best_rho))
                    else :
                        print("Only kinship has effect.")
                    
            else:
                lmm = LMM(phenotype, cov_matrix, QS)
                if not mixed:
                    lmm.delta = 1
                    lmm.fix('delta')
                lmm.fit(verbose=False)
            if debugger:
                fun_end = time.time()
                print(" Computing Null model took {}".format(fun_end-fun_start))
            ##########################################################################################################################################################
            
            if regressCovariatesUpfront:
                #pdb.set_trace();
                phenotype_corrected = phenotype-cov_matrix[:,1:].dot(lmm.beta[1:])
                cov_matrix_corrected = cov_matrix[:,0]
                if randomeff_mix:
                    lmm = LMM(phenotype_corrected, cov_matrix_corrected, Sigma_qs[mixingParameters["rho"]])
                else:
                    lmm = LMM(phenotype_corrected, cov_matrix_corrected, QS)
                lmm.fit(verbose=False)
            
            null_lml = lmm.lml()
            flmm = lmm.get_fast_scanner()
            countChunker = 0
            #pdb.set_trace();
            for snpGroup in utils.chunker(snpQuery, blocksize):
                countChunker=countChunker+1
                #print(countChunker)
                #Fix seed at the start of the first chunker so all permutations are based on the same random first split.
                np.random.seed(seed)
                
                snp_names = snpGroup
                
                tested_snp_names.extend(snp_names)
                snp_matrix_DF = risk_df.loc[snp_names,individual_ids].transpose()
                ##GRS var QC
                snp_matrix_DF = snp_matrix_DF.loc[:,snp_matrix_DF.isna().sum(axis=0)!=snp_matrix_DF.shape[0],]
                snp_matrix_DF = snp_matrix_DF.loc[:,(np.nanstd(snp_matrix_DF,axis=0)>0)]
#               test if the covariates, kinship, snp and phenotype are in the same order
                if (len(snp_matrix_DF.index) != len(sample2individual_feature.loc[phenotype_ds.index]['iid']) or not all(snp_matrix_DF.index==sample2individual_feature.loc[phenotype_ds.index]['iid'])):
                    print ('There is an issue in mapping phenotypes and genotypes')
                    sys.exit()
                #Impute missingness
                #pdb.set_trace()
                call_rate = 1-snp_matrix_DF.isnull().sum()/len(snp_matrix_DF.index)
                if snpQcInfo is None and call_rate is not None:
                    snpQcInfo = call_rate
                elif call_rate is not None:
                    snpQcInfo = pd.concat([snpQcInfo, call_rate], axis=0)
                selection = call_rate > min_call_rate
                snp_matrix_DF = snp_matrix_DF.loc[:,list(snp_matrix_DF.columns[selection])]
                if  snp_matrix_DF.shape[1]==0:
                    continue
                
                snp_matrix_DF = pd.DataFrame(fill_NaN.fit_transform(snp_matrix_DF),index=snp_matrix_DF.index,columns=snp_matrix_DF.columns)
                #
                
                G = snp_matrix_DF.values
                G = G.astype(float)
                G_index = snp_matrix_DF.columns
                
                scannerOut = flmm.fast_scan(G, verbose=False)
                alt_lmls = scannerOut['lml']
                effsizes = scannerOut['effsizes1']
                var_effsizes_se = scannerOut['effsizes1_se']
                var_pvalues = utils.lrt_pvalues(null_lml, alt_lmls)
                
                #add these results to qtl_results
                temp_df = pd.DataFrame(index = range(len(G_index)),columns=['feature_id','snp_id','p_value','beta','beta_se','empirical_feature_p_value'])
                temp_df['snp_id'] = G_index
                temp_df['feature_id'] = feature_id.replace("/","-")
                temp_df['beta'] = np.asarray(effsizes)
                temp_df['p_value'] = np.asarray(var_pvalues)
                temp_df['beta_se'] = np.asarray(var_effsizes_se)
                #insert default dummy value
                temp_df['empirical_feature_p_value'] = -1.0
                
                if(n_perm!=0):
                    pValueBuffer = []
                    if(write_zscore):
                        zScoreBuffer = []
                    totalSnpsToBeTested = (G.shape[1]*n_perm)
                    permutationStepSize = np.floor(n_perm/(totalSnpsToBeTested/blocksize))
                    if(permutationStepSize>n_perm):
                        permutationStepSize=n_perm
                    elif(permutationStepSize<1):
                        permutationStepSize=1
                    
                    if(write_permutations):
                        perm_df = pd.DataFrame(index = range(len(G_index)),columns=['snp_id'] + ['permutation_'+str(x) for x in range(n_perm)])
                        perm_df['snp_id'] = G_index
                    for currentNperm in utils.chunker(list(range(1, n_perm+1)), permutationStepSize):
                        if (kinship_df is not None) and (relatedness_score is not None):
                            temp = utils.get_shuffeld_genotypes_preserving_kinship(geneticaly_unique_individuals, relatedness_score, snp_matrix_DF,kinship_df.loc[individual_ids,individual_ids], len(currentNperm))
                        else :
                            temp = utils.get_shuffeld_genotypes(snp_matrix_DF, len(currentNperm))
                        temp = temp.astype(float)
                        scannerOut = flmm.fast_scan(temp, verbose=False)
                        alt_lmls_p = scannerOut['lml']
                        var_pvalues_p = utils.lrt_pvalues(null_lml, alt_lmls_p)
                        if(write_zscore):
                            zScoreBuffer.extend(np.asarray(scannerOut['effsizes1']/scannerOut['effsizes1_se']))
                        pValueBuffer.extend(np.asarray(var_pvalues_p))
                    if(not(len(pValueBuffer)==totalSnpsToBeTested)):
                        print(len(pValueBuffer))
                        print(pValueBuffer)
                        print(totalSnpsToBeTested)
                        print('Error in blocking logic for permutations.')
                        sys.exit()
                    perm = 0
                    for relevantOutput in utils.chunker(pValueBuffer,G.shape[1]) :
                        if(write_permutations):
                            if(write_zscore):
                                sPos = 0 + G.shape[1]*perm
                                ePos = sPos + G.shape[1]
                                perm_df['permutation_'+str(perm)] = zScoreBuffer[sPos:ePos]
                            else :
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
        if not data_written :
            fail_qc_features.append(feature_id)
        else:
            n_samples.append(phenotype_ds.size)
            n_e_samples.append(len(geneticaly_unique_individuals))
            if(n_perm>1):
                #updated_permuted_p_in_hdf5(bestPermutationPval, feature_id);
                alpha_para, beta_para = output_writer.apply_pval_correction(feature_id.replace("/","-"),bestPermutationPval,False)
                if write_feature_top_permutations:
                    np.savetxt(output_dir+"/Permutation.pValues."+feature_id.replace("/","-")+".txt",bestPermutationPval)
                alpha_params.append(alpha_para)
                beta_params.append(beta_para)
            if randomeff_mix :
                random_eff_param.append(feature_best_rho)
            
        if contains_missing_samples:
            QS = None
            Sigma_qs = None
            geneticaly_unique_individuals = tmp_unique_individuals
            snpQcInfo = snpQcInfo.to_frame(name="call_rate")
            snpQcInfo.index.name = "snp_id"
            snpQcInfo.to_csv(output_dir+'/snp_qc_metrics_naContaining_feature_{}.txt'.format(feature_id),sep='\t')
            del tmp_unique_individuals
        else:
            if (snpQcInfo is not None and snpQcInfoMain is not None):
                snpQcInfoMain = pd.concat([snpQcInfoMain, snpQcInfo], axis=0)
            elif snpQcInfo is not None :
                snpQcInfoMain = snpQcInfo.copy(deep=True)
            #print('step 5')
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
    #gather unique indexes of tested snps
    #write annotation and snp data to file
    snp_df = pd.DataFrame()
    snp_df['snp_id'] = np.unique(tested_snp_names)
    snp_df.index = np.unique(tested_snp_names)
    snp_df['chromosome'] = "NA"
    snp_df['position'] = "NA"
    snp_df['assessed_allele'] = "NA"
    if (snpQcInfoMain is not None):
        snpQcInfoMain = snpQcInfoMain.to_frame(name="call_rate")
        snpQcInfoMain['index']=snpQcInfoMain.index
        snpQcInfoMain = snpQcInfoMain.drop_duplicates()
        del snpQcInfoMain['index']
        snp_df = pd.concat([snp_df, snpQcInfoMain.reindex(snp_df.index)], axis=1)
    
    feature_list = list(set(feature_list) - set(fail_qc_features))
    annotation_df = annotation_df.reindex(feature_list)
    annotation_df['n_samples'] = n_samples
    annotation_df['n_e_samples'] = n_e_samples
    if(n_perm>1):
        annotation_df['alpha_param'] = alpha_params
        annotation_df['beta_param'] = beta_params
    if randomeff_mix:
        annotation_df['rho'] = random_eff_param
    if not selectionStart is None :
        snp_df.to_csv(output_dir+'/snp_metadata_{}_{}_{}.txt'.format(chromosome,selectionStart,selectionEnd),sep='\t',index=False)
        annotation_df.to_csv(output_dir+'/feature_metadata_{}_{}_{}.txt'.format(chromosome,selectionStart,selectionEnd),sep='\t')
    else :
        snp_df.to_csv(output_dir+'/snp_metadata_{}.txt'.format(chromosome),sep='\t',index=False)
        annotation_df.to_csv(output_dir+'/feature_metadata_{}.txt'.format(chromosome),sep='\t')

if __name__=='__main__':
    args = qtl_parse_args.get_grsQtl_args()
    grsFile  = args.genetic_risk_scores
    anno_file = args.annotation_file
    pheno_file = args.phenotype_file
    output_dir = args.output_directory
    genetic_range = args.genomic_range
    covariates_file = args.covariates_file
    randeff_file = args.randomeff_files
    samplemap_file = args.sample_mapping_file
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
    write_permutations = args.write_permutations
    write_zscore = args.write_zscores
    includeAllChromsomes = args.no_chromosome_filter
    regressBefore = args.regress_covariates
    write_feature_top_permutations = args.write_feature_top_permutations
    lr_random_effect = args.low_rank_random_effect
    debugger = args.debugger
    
    if (random_seed is None):
        random_seed = np.random.randint(40000)
    if(write_zscore):
        #Need to make sure we are also outputing (z-score writing turns the normal writer into a z-score writer.)
        write_permutations = write_zscore
    if(n_perm==0 and write_permutations):
        write_permutations=False
    
    if(n_perm==1):
        print("Warning: With only 1 permutation P-value correction is not performed.")
    if(n_perm<50):
        print("Warning: With less than 50 permutations P-values correction is not very accurate.")
    run_PrsQtl_analysis(pheno_file, anno_file, grsFile, output_dir, min_call_rate=float(min_call_rate), blocksize=int(block_size), 
                        skipAutosomeFiltering= includeAllChromsomes, gaussianize_method = gaussianize,  minimum_test_samples= int(minimum_test_samples), seed=int(random_seed), 
                        n_perm=int(n_perm), write_permutations = write_permutations, write_zscore = write_zscore, write_feature_top_permutations = write_feature_top_permutations, relatedness_score=relatedness_score, feature_variant_covariate_filename = feature_variant_covariate_filename, 
                        snps_filename=snps_filename, feature_filename=feature_filename, snp_feature_filename=snp_feature_filename, genetic_range=genetic_range, covariates_filename=covariates_file, 
                        randomeff_filename=randeff_file, sample_mapping_filename=samplemap_file, regressCovariatesUpfront = regressBefore, lr_random_effect = lr_random_effect, debugger= debugger)