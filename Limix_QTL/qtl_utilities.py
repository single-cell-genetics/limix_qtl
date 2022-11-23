import os
import sys
import gc
import numpy as np
import pandas as pd
import math
import scipy.stats as scst
import scipy as sp
import scipy.linalg as la
from scipy.stats import chi2
import qtl_loader_utils
import pdb
from glimix_core.lmm import LMM
from glimix_core.glmm._glmm import GLMM
from numpy.linalg import eigh, svd, pinv, solve, norm
import dask.array as da

def run_QTL_analysis_load_intersect_phenotype_covariates_kinship_sample_mapping(pheno_filename, anno_filename, geno_prefix,
        plinkGenotype, minimum_test_samples= 10, relatedness_score=None, cis_mode=True, skipAutosomeFiltering = False, snps_filename=None,
        feature_filename=None, snp_feature_filename=None, selection='all', covariates_filename=None, randomeff_filename=None,
        sample_mapping_filename=None, extended_anno_filename=None, feature_variant_covariate_filename=None):


    # pheno_filename = "/Users/chaaya/dhonveli_dkfz/hipsci_pipeline/geuvadis_CEU_test_data/Expression/Geuvadis_CEU_YRI_Expr.txt.gz"
    # anno_filename = "/Users/chaaya/dhonveli_dkfz/hipsci_pipeline/geuvadis_CEU_test_data/Expression/Geuvadis_CEU_Annot.txt"
    # geno_prefix = "/Users/chaaya/dhonveli_dkfz/limix_qtl/limix_qtl-master/Limix_QTL/test_data/Genotypes/Geuvadis"
    # plinkGenotype = "/Users/chaaya/dhonveli_dkfz/limix_qtl/limix_qtl-master/Limix_QTL/test_data/Genotypes/Geuvadis"
    # minimum_test_samples = 10
    # relatedness_score = 0.95
    # cis_mode = True
    # skipAutosomeFiltering = False
    # snps_filename = None
    # feature_filename = None
    # snp_feature_filename = None
    # selection = 'all'
    # covariates_filename = "/Users/chaaya/dhonveli_dkfz/limix_qtl/limix_qtl-master/Limix_QTL/test_data/Expression/Geuvadis_CEU_YRI_covariates.txt"
    # randomeff_filename = "/Users/chaaya/dhonveli_dkfz/limix_qtl/limix_qtl-master/Limix_QTL/test_data/Genotypes/Geuvadis_chr1_kinship.normalized.txt,/Users/chaaya/dhonveli_dkfz/limix_qtl/limix_qtl-master/Limix_QTL/test_data/Genotypes/Geuvadis_readdepth.txt"
    # sample_mapping_filename = "/Users/chaaya/dhonveli_dkfz/limix_qtl/limix_qtl-master/Limix_QTL/test_data/Geuvadis_CEU_gte.txt"
    # extended_anno_filename = None
    # feature_variant_covariate_filename = None
    # output_dir = "/Users/chaaya/dhonveli_dkfz/limix_qtl/limix_qtl-master/Limix_QTL/test_data/Output2/"
    # window_size = 250000
    # min_maf = 0.05
    # min_hwe_P = 0.001
    # min_call_rate = 0.95
    # blocksize = 1000
    # gaussianize_method = None
    # genetic_range = "all"
    # seed = np.random.randint(40000)
    # n_perm = 0
    # write_permutations = False
    # regressCovariatesUpfront = False
    # write_feature_top_permutations = False

    
    # selection based on coordinates
    selectionStart = None
    selectionEnd = None
    if(":" in selection):
        parts = selection.split(":")
        if("-" not in parts[1]):
            print("No correct sub selection.")
            print("Given in: "+selection)
            print("Expected format: (chr number):(start location)-(stop location)")
            sys.exit()
        chromosome = parts[0]
        if("-" in parts[1]):
            parts2 = parts[1].split("-")
            selectionStart = int(parts2[0])
            selectionEnd = int(parts2[1])
    else :
        chromosome=selection

    ''' function to take input and intersect sample and genotype.'''
    #Load input data files & filter for relevant data
    #Load input data filesf

    # loading phenotype and annotation files
    phenotype_df = qtl_loader_utils.get_phenotype_df(pheno_filename)
    annotation_df = qtl_loader_utils.get_annotation_df(anno_filename)

    phenotype_df.columns = phenotype_df.columns.astype("str")
    phenotype_df.index = phenotype_df.index.astype("str")
    annotation_df.columns = annotation_df.columns.astype("str")
    annotation_df.index = annotation_df.index.astype("str")

    # loading genotype, file type is handeled in the loader utils.
    bim,fam,bed,bgen = qtl_loader_utils.get_genotype_data(geno_prefix, plinkGenotype)

    # converting chromsome names
    annotation_df.replace(['X', 'Y', 'XY', 'MT'], ['23', '24', '25', '26'],inplace=True)
    if chromosome=='X' :
        chromosome = '23'
    elif chromosome=='Y':
        chromosome = '24'
    elif chromosome=='XY':
        chromosome='25'
    elif chromosome=='MT':
        chromosome='26'

    print("Intersecting data.")
    if(annotation_df.shape[0] != annotation_df.groupby(annotation_df.index).first().shape[0]):
        print("Only one location per feature supported. If multiple locations are needed please look at: --extended_anno_file")
        sys.exit()

    ##Make sure that there is only one entry per feature id!.

    sample2individual_df = qtl_loader_utils.get_samplemapping_df(sample_mapping_filename,list(phenotype_df.columns),'sample')
    sample2individual_df.index = sample2individual_df.index.astype('str')
    sample2individual_df = sample2individual_df.astype('str')
    sample2individual_df['sample']=sample2individual_df.index
    sample2individual_df = sample2individual_df.drop_duplicates();

    ##Filter first the linking files!
    #Subset linking to relevant genotypes.
    orgSize = sample2individual_df.shape[0]
    sample2individual_df = sample2individual_df.loc[sample2individual_df['iid'].map(lambda x: x in list(map(str, fam.index))),:]
    diff = orgSize- sample2individual_df.shape[0]
    orgSize = sample2individual_df.shape[0]
    print("Dropped: "+str(diff)+" samples because they are not present in the genotype file.")

    #Subset linking to relevant phenotypes.
    sample2individual_df = sample2individual_df.loc[np.intersect1d(sample2individual_df.index,phenotype_df.columns),:]
    diff = orgSize- sample2individual_df.shape[0]
    orgSize = sample2individual_df.shape[0]
    print("Dropped: "+str(diff)+" samples because they are not present in the phenotype file.")

    #Subset linking vs kinship.
    kinship_df = None
    readdepth_df = None
    if randomeff_filename is not None:
        kinship_df,readdepth_df = qtl_loader_utils.get_randeff_df(randomeff_filename)

    if kinship_df is not None:
        #Filter from individual2sample_df & sample2individual_df since we don't want to filter from the genotypes.
        sample2individual_df = sample2individual_df[sample2individual_df['iid'].map(lambda x: x in list(map(str, kinship_df.index)))]
        diff = orgSize- sample2individual_df.shape[0]
        orgSize = sample2individual_df.shape[0]
        print("Dropped: "+str(diff)+" samples because they are not present in the kinship file.")

    if readdepth_df is not None:
        #This needs to come from the covariate site not the genotype side!
        #Filter from individual2sample_df & sample2individual_df since we don't want to filter from the genotypes.
        sample2individual_df = sample2individual_df[sample2individual_df['sample'].map(lambda x: x in list(map(str, readdepth_df.index)))]
        diff = orgSize- sample2individual_df.shape[0]
        orgSize = sample2individual_df.shape[0]
        print("Dropped: "+str(diff)+" samples because they are not present in the second random effect file.")

    #Subset linking vs covariates.
    covariate_df = qtl_loader_utils.get_covariate_df(covariates_filename)

    if covariate_df is not None:
        if np.nansum(covariate_df==1,0).max()<covariate_df.shape[0]: covariate_df.insert(0, 'ones', np.ones(covariate_df.shape[0]))
        sample2individual_df = sample2individual_df.loc[list(set(sample2individual_df.index) & set(covariate_df.index)),:]
        diff = orgSize- sample2individual_df.shape[0]
        orgSize = sample2individual_df.shape[0]
        print("Dropped: "+str(diff)+" samples because they are not present in the covariate file.")
    ###

    print("Number of samples with genotype & phenotype data: " + str(sample2individual_df.shape[0]))
    if(sample2individual_df.shape[0]<minimum_test_samples):
        print("Not enough samples with both genotype & phenotype data.")
        sys.exit()
        
    ##Filter now the actual data!
    #Filter phenotype data based on the linking files.
    phenotype_df = phenotype_df.loc[list(set(phenotype_df.index)&set(annotation_df.index)),sample2individual_df.index.values]

    #Filter kinship data based on the linking files.
    genetically_unique_individuals = None
    if kinship_df is not None:
        kinship_df = kinship_df.loc[np.intersect1d(kinship_df.index,sample2individual_df['iid']),np.intersect1d(kinship_df.index,sample2individual_df['iid'])]
    if (kinship_df is not None) and (relatedness_score is not None):
        genetically_unique_individuals = get_unique_genetic_samples(kinship_df, relatedness_score);
    
    #Filter covariate data based on the linking files.

    snp_feature_filter_df= qtl_loader_utils.get_snp_feature_df(snp_feature_filename)

    try:
        feature_filter_df = qtl_loader_utils.get_snp_df(feature_filename)
    except:
        if feature_filename  is not None:
            feature_filter_df=pd.DataFrame(index=feature_filename)
    #Do filtering on features.
    if feature_filter_df is not None:
        phenotype_df = phenotype_df.loc[feature_filter_df.index,:]
        ##Filtering on features to test.
    if snp_feature_filter_df is not None:
        lst3 = set(phenotype_df.index).intersection(np.unique(snp_feature_filter_df['feature_id']))
        phenotype_df = phenotype_df.loc[lst3,:]
        ##Filtering on features  to test from the combined feature snp filter.

    if ((not cis_mode) and len(set(bim['chrom']))<22) :
        print("Warning, running a trans-analysis on snp data from less than 22 chromosomes.\nTo merge data later the permutation P-values need to be written out.")

    if(cis_mode):
        #Remove features from the annotation that are on chromosomes which are not present anyway.
        annotation_df = annotation_df[np.in1d(annotation_df['chromosome'],list(set(bim['chrom'])))]

    #Prepare to filter on snps.
    snp_filter_df = qtl_loader_utils.get_snp_df(snps_filename)
    if snp_filter_df is not None:
        toSelect = set(snp_filter_df.index).intersection(set(bim['snp']))
        bim = bim.loc[bim['snp'].isin(toSelect)]
        ##Filtering on SNPs to test from the snp filter.

    if snp_feature_filter_df is not None:
        toSelect = set(np.unique(snp_feature_filter_df['snp_id'])).intersection(set(bim['snp']))
        bim = bim.loc[bim['snp'].isin(toSelect)]
        ##Filtering on features  to test from the combined feature snp filter.

    #Filtering for sites on non allosomes.
    if not skipAutosomeFiltering :
        annotation_df = annotation_df[annotation_df['chromosome'].map(lambda x: x in list(map(str, range(1, 23))))]

    #Determine features to be tested
    if chromosome=='all':
        feature_list = list(set(annotation_df.index)&set(phenotype_df.index))
    else:
        if not selectionStart is None :
            lowest = min([selectionStart,selectionEnd])
            highest = max([selectionStart,selectionEnd])
            annotation_df['mean'] = ((annotation_df["start"] + annotation_df["end"])/2)
            feature_list = list(set(annotation_df.iloc[(annotation_df['chromosome'].values==chromosome) & (annotation_df['mean'].values>=lowest) & (annotation_df["mean"].values<highest)].index.values)&set(phenotype_df.index))
            del annotation_df['mean']
        else :
            feature_list = list(set(annotation_df[annotation_df['chromosome']==chromosome].index)&set(phenotype_df.index))
    #Drop not used feature information.
    phenotype_df = phenotype_df.loc[feature_list,:]
    gc.collect()
    print("Number of features to be tested: " + str(len(feature_list)))
    print("Total number of variants to be considered, before variante QC and feature intersection: " + str(bim.shape[0]))

    if(phenotype_df.shape[1]<minimum_test_samples):
        print("Not enough samples with both genotype & phenotype data, for current number of covariates.")
        sys.exit()

    if extended_anno_filename is not None:
        complete_annotation_df = pd.read_csv(extended_anno_filename,sep='\t',index_col=0)
        annotation_df['index']=annotation_df.index
        complete_annotation_df['index']=complete_annotation_df.index
        complete_annotation_df = pd.concat([annotation_df,complete_annotation_df]).drop_duplicates()
        del complete_annotation_df['index']
    else:
        complete_annotation_df = annotation_df

    feature_variant_covariate_df = qtl_loader_utils.get_snp_feature_df(feature_variant_covariate_filename)
    
    if len(bim["snp"].values) > len(set(bim["snp"].values)):
        print("Warning duplicated SNP ids (After filtering if applied).")
        print("Removing variants observed twice.")
        snpC = bim["snp"].value_counts()
        snpC = snpC.index[np.where(snpC==1)].values
        bim = bim.loc[bim['snp'].isin(snpC)]
    
    return [phenotype_df, kinship_df,readdepth_df, covariate_df, sample2individual_df, complete_annotation_df, annotation_df, snp_filter_df,
        snp_feature_filter_df, genetically_unique_individuals, minimum_test_samples, feature_list, bim, fam, bed, bgen, chromosome,
        selectionStart, selectionEnd, feature_variant_covariate_df]

def run_PrsQtl_analysis_load_intersect_phenotype_covariates_kinship_sample_mapping\
        (pheno_filename, anno_filename, prsFile, minimum_test_samples= 10, relatedness_score=0.95, skipAutosomeFiltering = False, snps_filename=None,
         feature_filename=None, snp_feature_filename=None, selection='all', covariates_filename=None, randomeff_filename=None, sample_mapping_filename=None, feature_variant_covariate_filename=None):

    selectionStart = None
    selectionEnd = None
    if(":" in selection):
        parts = selection.split(":")
        if("-" not in parts[1]):
            print("No correct sub selection.")
            print("Given in: "+selection)
            print("Expected format: (chr number):(start location)-(stop location)")
            sys.exit()
        chromosome = parts[0]
        if("-" in parts[1]):
            parts2 = parts[1].split("-")
            selectionStart = int(parts2[0])
            selectionEnd = int(parts2[1])
    else :
        chromosome=selection

    ''' function to take input and intersect sample and genotype.'''
    #Load input data files & filter for relevant data
    #Load input data filesf
    #import pdb; pdb.set_trace();
    phenotype_df = qtl_loader_utils.get_phenotype_df(pheno_filename)
    annotation_df = qtl_loader_utils.get_annotation_df(anno_filename)

    phenotype_df.columns = phenotype_df.columns.astype("str")
    phenotype_df.index = phenotype_df.index.astype("str")
    annotation_df.columns = annotation_df.columns.astype("str")
    annotation_df.index = annotation_df.index.astype("str")

    if(annotation_df.shape[0] != annotation_df.groupby(annotation_df.index).first().shape[0]):
        print("Only one location per feature supported. If multiple locations are needed please look at: --extended_anno_file")
        sys.exit()

    #Determine features to be tested
    if chromosome!='all':
        if not selectionStart is None :
            lowest = min([selectionStart,selectionEnd])
            highest = max([selectionStart,selectionEnd])
            annotation_df['mean'] = ((annotation_df["start"] + annotation_df["end"])/2)
            feature_list = list(set(annotation_df.iloc[(annotation_df['chromosome'].values==chromosome) & (annotation_df['mean'].values>=lowest) & (annotation_df["mean"].values<highest)].index.values))
            annotation_df = annotation_df.loc[feature_list,]
            del annotation_df['mean']
        else :
            feature_list = list(annotation_df[annotation_df['chromosome']==chromosome].index)
            annotation_df = annotation_df.loc[feature_list,]

    #To be able to read variants from a large file we change the loading here.
    #First we subset the genes to the chunk and get the relevant SNPs based on that.

    snp_feature_filter_df= qtl_loader_utils.get_snp_feature_df(snp_feature_filename)
    feature_filter_df = qtl_loader_utils.get_snp_df(feature_filename)
    snp_filter_df = qtl_loader_utils.get_snp_df(snps_filename)
    feature_variant_covariate_df = qtl_loader_utils.get_snp_feature_df(feature_variant_covariate_filename)

    #import pdb; pdb.set_trace()
    #Do filtering on variants and features first stage.
    if snp_feature_filter_df is not None:
        if feature_filter_df is not None:
            toSelect = set(feature_filter_df.index.values).intersection(set(annotation_df.index.values))
            annotation_df = annotation_df.loc[toSelect,]
        toSelect = list(set(snp_feature_filter_df['feature_id'].values).intersection(set(annotation_df.index.values)))
        snp_feature_filter_df = snp_feature_filter_df.loc[snp_feature_filter_df['feature_id'].isin(toSelect)]
        relSnps = snp_feature_filter_df['snp_id'].values

        if snp_filter_df is not None:
            relSnps = set(snp_filter_df.index).intersection(set(relSnps))
        if feature_variant_covariate_df is not None:
            feature_variant_covariate_df = feature_variant_covariate_df.loc[feature_variant_covariate_df['feature_id'].isin(toSelect)]
            relSnps = np.union1d(relSnps, feature_variant_covariate_df["snp_id"].values)

        relSnps = np.unique(relSnps)
        risk_df = qtl_loader_utils.get_grs_subset_df(prsFile, relSnps)

        if risk_df is None:
            print("No variants selected during SNP reading.")
            sys.exit()
        risk_df = risk_df.assign(SnpId=risk_df.index.values)
        risk_df = risk_df.drop_duplicates(keep='first')
        risk_df = risk_df.drop(['SnpId'], axis='columns')
        risk_df = risk_df.loc[risk_df.isnull().sum(axis=1)!=risk_df.shape[1],]
    elif snp_filter_df is not None:
        relSnps = snp_filter_df.index

        if feature_variant_covariate_df is not None:
            feature_variant_covariate_df = feature_variant_covariate_df.loc[feature_variant_covariate_df['feature_id'].isin(toSelect)]
            relSnps = np.union1d(relSnps, feature_variant_covariate_df["snp_id"].values)

        relSnps = np.unique(relSnps)
        risk_df = qtl_loader_utils.get_grs_subset_df(prsFile, relSnps)
        if risk_df is None:
            print("No variants selected during SNP reading.")
            sys.exit()
        risk_df = risk_df.assign(SnpId=risk_df.index.values)
        risk_df = risk_df.drop_duplicates(keep='first')
        risk_df = risk_df.drop(['SnpId'], axis='columns')
        risk_df = risk_df.loc[risk_df.isnull().sum(axis=1)!=risk_df.shape[1],]
    else :
        risk_df = qtl_loader_utils.get_phenotype_df(prsFile)
    print("Intersecting data.")
    risk_df =  risk_df.astype(float)
    #pdb.set_trace();
    ##Make sure that there is only one entry per feature id!.

    sample2individual_df = qtl_loader_utils.get_samplemapping_df(sample_mapping_filename,list(phenotype_df.columns),'sample')
    sample2individual_df.index = sample2individual_df.index.astype('str')
    sample2individual_df = sample2individual_df.astype('str')
    sample2individual_df['sample']=sample2individual_df.index
    sample2individual_df = sample2individual_df.drop_duplicates();
    ##Filter first the linking files!
    #Subset linking to relevant genotypes.
    orgSize = sample2individual_df.shape[0]
    sample2individual_df = sample2individual_df.loc[sample2individual_df['iid'].map(lambda x: x in list(map(str, risk_df.columns))),:]
    diff = orgSize- sample2individual_df.shape[0]
    orgSize = sample2individual_df.shape[0]
    print("Dropped: "+str(diff)+" samples because they are not present in the genotype file.")

    #Subset linking to relevant phenotypes.
    sample2individual_df = sample2individual_df.loc[np.intersect1d(sample2individual_df.index,phenotype_df.columns),:]
    diff = orgSize- sample2individual_df.shape[0]
    orgSize = sample2individual_df.shape[0]
    print("Dropped: "+str(diff)+" samples because they are not present in the phenotype file.")
    #Subset linking vs kinship.
    # extract filename from randomeffects filename
    kinship_df = None
    readdepth_df = None
    if randomeff_filename is not None:
        kinship_df,readdepth_df = qtl_loader_utils.get_randeff_df(randomeff_filename)

    if kinship_df is not None:
        #Filter from individual2sample_df & sample2individual_df since we don't want to filter from the genotypes.
        sample2individual_df = sample2individual_df[sample2individual_df['iid'].map(lambda x: x in list(map(str, kinship_df.index)))]
        diff = orgSize- sample2individual_df.shape[0]
        orgSize = sample2individual_df.shape[0]
        print("Dropped: "+str(diff)+" samples because they are not present in the kinship file.")

    if readdepth_df is not None:
        #This needs to come from the covariate site not the genotype side!
        #Filter from individual2sample_df & sample2individual_df since we don't want to filter from the genotypes.
        sample2individual_df = sample2individual_df[sample2individual_df['sample'].map(lambda x: x in list(map(str, readdepth_df.index)))]
        diff = orgSize- sample2individual_df.shape[0]
        orgSize = sample2individual_df.shape[0]
        print("Dropped: "+str(diff)+" samples because they are not present in the second random effect file.")

    #Subset linking vs covariates.
    covariate_df = qtl_loader_utils.get_covariate_df(covariates_filename)

    if covariate_df is not None:
        if np.nansum(covariate_df==1,0).max()<covariate_df.shape[0]: covariate_df.insert(0, 'ones',np.ones(covariate_df.shape[0]))
        sample2individual_df = sample2individual_df.loc[list(set(sample2individual_df.index) & set(covariate_df.index)),:]
        diff = orgSize- sample2individual_df.shape[0]
        orgSize = sample2individual_df.shape[0]
        print("Dropped: "+str(diff)+" samples because they are not present in the covariate file.")

    ###
    print("Number of samples with genotype & phenotype data: " + str(sample2individual_df.shape[0]))
    if(sample2individual_df.shape[0]<minimum_test_samples):
        print("Not enough samples with both genotype & phenotype data.")
        sys.exit()
    #import pdb; pdb.set_trace()
    ##Filter now the actual data!
    #Filter phenotype data based on the linking files.
    phenotype_df = phenotype_df.loc[list(set(phenotype_df.index)&set(annotation_df.index)),sample2individual_df.index.values]

    #Filter kinship data based on the linking files.
    genetically_unique_individuals = None
    if kinship_df is not None:
        kinship_df = kinship_df.loc[np.intersect1d(kinship_df.index,sample2individual_df['iid']),np.intersect1d(kinship_df.index,sample2individual_df['iid'])]
    if kinship_df is not None and (relatedness_score is not None):
        genetically_unique_individuals = get_unique_genetic_samples(kinship_df, relatedness_score);

    #Filter covariate data based on the linking files.

    #Do filtering on features.
    if feature_filter_df is not None:
        toSelect = set(feature_filter_df.index.values).intersection(set(phenotype_df.index.values))
        phenotype_df = phenotype_df.loc[toSelect,:]
        ##Filtering on features to test.
    if snp_feature_filter_df is not None:
        toSelect = set(snp_feature_filter_df['feature_id'].values).intersection(set(phenotype_df.index.values))
        phenotype_df = phenotype_df.loc[toSelect,:]
        if feature_filter_df is not None:
            snp_feature_filter_df = snp_feature_filter_df.loc[snp_feature_filter_df['feature_id'].isin(toSelect)]
        ##Filtering on features  to test from the combined feature snp filter.

    #Prepare to filter on SNPs.
    if snp_filter_df is not None:
        toSelect = set(snp_filter_df.index).intersection(set(risk_df.index.values))
        risk_df=risk_df.loc[toSelect,:]
        ##Filtering on SNPs to test from the snp filter.

    if snp_feature_filter_df is not None:
        toSelect = set(np.unique(snp_feature_filter_df['snp_id'])).intersection(set(risk_df.index.values))
        risk_df=risk_df.loc[toSelect,:]
        ##Filtering on features to test from the combined feature snp filter.

    #Filtering for sites on non allosomes.
    if not skipAutosomeFiltering :
        annotation_df = annotation_df[annotation_df['chromosome'].map(lambda x: x in list(map(str, range(1, 23))))]

    feature_list = list(set(annotation_df.index)&set(phenotype_df.index))
    print("Number of features to be tested: " + str(len(feature_list)))
    print("Total number of variants to be considered, before variante QC and feature intersection: " + str(risk_df.shape[0]))

    if(phenotype_df.shape[1]<minimum_test_samples):
        print("Not enough samples with both genotype & phenotype data, for current number of covariates.")
        sys.exit()

    return [phenotype_df, kinship_df, readdepth_df, covariate_df, sample2individual_df, annotation_df, snp_filter_df, snp_feature_filter_df, genetically_unique_individuals, minimum_test_samples, feature_list, risk_df, chromosome, selectionStart, selectionEnd, feature_variant_covariate_df]

def merge_QTL_results(results_dir):
    '''Merge QTL results for individual chromosomes into a combined, indexed
    hdf5 file.'''
    qtl_results_files = sorted(glob.glob(results_dir+'qtl_results_*.txt'))

    hdf5_outfile = qtl_output.hdf5_writer(results_dir+'qtl_results.h5')

    for filename in qtl_results_files:
        df = pd.read_csv(filename,sep='\t')
        hdf5_outfile.add_result_df(df)

    hdf5_outfile.close()

def chunker(seq, size):
    return (seq[pos:pos + np.int(size)] for pos in range(0, len(seq), np.int(size)))

def get_unique_genetic_samples(kinship_df, relatedness_score):
#    tril returns the lower triungular.
#    if two lines are > identity then  kinship_df>=relatedness_score should have an offdiagonal 1.
#    if there is one 1  in the tril then it means that  in the upper triul there was a line withe identical genotype
    return (kinship_df.index[(np.tril(kinship_df>=relatedness_score,k=-1)).sum(1)==0])

def force_normal_distribution(phenotype, method='gaussnorm', reference=None):
    _doc='rank transform x into ref/ gaussian;keep the range; keep ties'

    if method=='log':
        return np.log(1+phenotype)

    if method=='log_standardize':
        temp=np.log(1+phenotype)
        return (temp-np.nanmean(temp))/np.nanstd(temp)

    if method=='arcsin':
        return np.arcsin(np.sqrt(phenotype))

    if method=='arcsin_standardize':
        temp=np.arcsin(np.sqrt(phenotype))
        return (temp-np.nanmean(temp))/np.nanstd(temp)

    if method=='standardize':
        return (phenotype-np.nanmean(phenotype))/np.nanstd(phenotype)

    indextoupdate = np.isfinite(phenotype)
    y1 = phenotype[indextoupdate]
    yuni,yindex=np.unique(y1, return_inverse=True)
    phenotypenorm=phenotype.copy()

    if method =='gaussnorm':

        sref = scst.norm.isf(np.linspace(0.001, 0.999,num=yuni.shape[0])[::-1])
        phenotypenorm[indextoupdate]=sref[yindex]
        return phenotypenorm

    elif method=='ranknorm':
        try:
            xref1=np.unique(reference[np.isfinite(reference)])
            sref=np.sort(xref1)[np.linspace(0,xref1.shape[0]-0.001, num=y1.shape[0]).astype(int)]
        except:
            print ('reference missing. provide reference to force_normal_distribution or choose gaussnorm')
            return 1
        phenotypenorm[indextoupdate]=sref[np.argsort(np.argsort(y1))]
        return phenotypenorm

    elif method=='ranknorm_duplicates':
        try:
            xref1=np.unique(reference[np.isfinite(reference)])### unique values from reference
            sref=np.sort(xref1)[np.linspace(0,xref1.shape[0]-0.001, num=yuni.shape[0]).astype(int)]
        except:
            print ('reference missing. provide reference to force_normal_distribution or choose gaussnorm')
            return 1
        phenotypenorm[indextoupdate]=sref[yindex]
        return phenotypenorm

    else:
        print ('methods are: log, log_standardize, standardize, gaussnorm, ranknorm, ranknorm_duplicates, arcsin, arcsin_standardize')

#get_shuffeld_genotypes_preserving_kinship(genetically_unique_individuals, relatedness_score, snp_matrix_DF,kinship_df.loc[individual_ids,individual_ids], n_perm)
def get_shuffeld_genotypes_preserving_kinship(genetically_unique_individuals, relatedness_score, snp_matrix_DF,kinship_df1,n_perm):

    # take only one line for replicates (those with the same name)
    boolean_selection = ~snp_matrix_DF.index.duplicated()
    temp = snp_matrix_DF.loc[boolean_selection,:]

    boolean_selection = ~kinship_df1.index.duplicated()
    kinship_df1 = kinship_df1.loc[boolean_selection, boolean_selection]

    # subset snp matrix to genetically_unique_individuals
    u_snp_matrix = temp.loc[genetically_unique_individuals,:]

    '''has replicates but not same lines form donor (np.setdiff1d(individual_ids,genetically_unique_individuals))'''
    #Shuffle and reinflate
    locationBuffer = np.zeros(snp_matrix_DF.shape[0], dtype=np.int)
    #Prepare location search for permuted snp_matrix_df.
    index_samples = np.arange(u_snp_matrix.shape[0])
    for index,current_name in enumerate(genetically_unique_individuals):
        # find all samples that are related to current_name, or are current_name itself.
        kinship_row = kinship_df1.loc[current_name]
        selection = np.logical_or(kinship_row>=relatedness_score, kinship_row.index==current_name)
        locationBuffer[np.where(selection)] = index 
    snp_matrix_copy = np.zeros((snp_matrix_DF.shape[0],snp_matrix_DF.shape[1]*n_perm))
    counter = 0
    end = (snp_matrix_DF.shape[1])
    for perm_id in range(0,n_perm) :
        np.random.shuffle(index_samples)
        temp_u = u_snp_matrix.values[index_samples,:]
        snp_matrix_copy[:,counter:end] = temp_u[locationBuffer,:]
        counter+= snp_matrix_DF.shape[1]
        end+= snp_matrix_DF.shape[1]
    return(snp_matrix_copy)

def get_shuffeld_genotypes(snp_matrix_DF,n_perm):
    snp_matrix_copy = np.zeros((snp_matrix_DF.shape[0],snp_matrix_DF.shape[1]*n_perm))
    counter = 0
    end = (snp_matrix_DF.shape[1])

    index_samples = np.arange(snp_matrix_DF.shape[0])
    for perm_id in range(0,n_perm) :
        np.random.shuffle(index_samples)
        snp_matrix_copy[:,counter:end] = snp_matrix_DF.values[index_samples,:]
        counter+= snp_matrix_DF.shape[1]
        end+= snp_matrix_DF.shape[1]
    return(snp_matrix_copy)

def qtl_plot(snp_matrix_DF, phenotype,K=None, M=None,LMM=None,snp=None,show_reg_cov=True):
    if LMM is None:
        LMM = limix.qtl.qtl_test_lmm(snp_matrix_DF.values, phenotype,K=K,M=M,verbose=False)

    if snp is None:
        snp=snp_matrix_DF.values[:,np.argmin(LMM.variant_effsizes)]

    cov=LMM.null_covariate_effsizes;betacov=np.array([cov[key] for key in cov.keys()])


    temp=pd.DataFrame(data=np.hstack([np.vstack([snp,phenotype,np.zeros(phenotype.shape[0]).astype(bool)]),\
                                          np.vstack([snp,phenotype-np.dot(M,betacov[:,None])[:,0],np.ones(phenotype.shape[0]).astype(bool)])]).T,\
                                            columns=['Variant','Phenotype donor','Covariates regressed'])

    temp['Covariates regressed']=temp['Covariates regressed'].astype(bool)

    indexunique=np.unique(np.array([l.split('-')[1] for l in snp_matrix_DF.index]),return_index=1)[1]
    ax = sb.boxplot(y=temp['Phenotype donor'],x=temp['Variant'])
    for patch in ax.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .03))

#    sb.swarmplot(y=temp['Phenotype donor'],x=temp['Variant'],color='grey')
    if show_reg_cov:
        index1=np.hstack([indexunique,indexunique+phenotype.shape[0]])
        sb.swarmplot(y=temp['Phenotype donor'].iloc[index1],x=temp['Variant'].iloc[index1],  hue=temp['Covariates regressed'].iloc[index1], split=True)
    else:
        index1= indexunique
        sb.swarmplot(y=temp['Phenotype donor'].iloc[index1],x=temp['Variant'].iloc[index1])

def do_snp_selection(feature_id, annotation_df, bim, cis_mode, window_size, skipAutosomeFiltering = False):
    annotation_sub_df = annotation_df.loc[[feature_id],:]
    list_of_snp_dfs = []
    snpQuery = None
    for feature_id,annotation_ds in annotation_sub_df.iterrows():
        chrom = str(annotation_ds.loc['chromosome'])
        start = annotation_ds.loc['start']
        end = annotation_ds.loc['end']
        # make robust to features selfpecified back-to-front
        lowest = min([start,end])
        highest = max([start,end])
        if (cis_mode) :
            # for cis, we sequentially add snps that fall within each region
            snpQuery = bim.query("chrom == '%s' & pos > %d & pos < %d" % (chrom, lowest-window_size, highest+window_size))
            list_of_snp_dfs.append(snpQuery)
        else :
            # for trans, we sequentially exclude snps that fall within each region
            if snpQuery is None :
                #Initial search for trans window.
                snpQuery = bim.query("(chrom == '%s' & (pos < %d | pos > %d))|chrom != '%s'" % (chrom, lowest-window_size, highest+window_size,chrom))
            else :
                #Here we start intersecting regions we want to exclude.
                snpQuery = snpQuery.query("(chrom == '%s' & (pos < %d | pos > %d))|chrom != '%s'" % (chrom, lowest-window_size, highest+window_size,chrom))

    if (cis_mode):
        selected_snp_df = pd.concat(list_of_snp_dfs).drop_duplicates()
    else:
        selected_snp_df = snpQuery
    if not skipAutosomeFiltering :
        # filtering for sites on non allosomes.
        selected_snp_df = selected_snp_df.loc[selected_snp_df['chrom'].map(lambda x: x in list(map(str, range(1, 23)))).values]

    return selected_snp_df

def reduce_snp(snp_df, threshold=0.975):
    ''' input a snp df  samples(rows) x snps( columns)'''
    ''' returns a df with columns: 'lead_snp_id' the pruned snp names and snp_id'''
    allsnps=snp_df.columns
    cor=abs(np.corrcoef(snp_df.T))>=threshold
    cor1=pd.DataFrame(data=np.triu(cor,k=1),index=allsnps,columns=allsnps).sum(1)
    cor2=pd.DataFrame(data=np.triu(cor,k=1),index=allsnps,columns=allsnps)
    uniquesnps=cor1[cor1==0].index.values
    duplicatedsnps=cor1[cor1>0].index.values
    '''the unique ones'''
    rez=pd.DataFrame(data=uniquesnps,index=uniquesnps,columns=['snp_id'])
    ''' the ones that are in LD.. if a,c and a,b but not b,c returns a,b)'''
    rez2=pd.DataFrame(data=duplicatedsnps,index=allsnps[np.argmax(cor2.loc[duplicatedsnps].values*cor2.loc[duplicatedsnps].values.sum(0),1)],columns=['snp_id'])
    rez=pd.concat([rez,rez2])
    rez['lead_snp_id']=rez.indexs
    return(rez)

def regressOut(Y, X):
    """
    regresses out X from Y
    """
    Xd = la.pinv(X)
    Y_out = Y-X.dot(Xd.dot(Y))
    return Y_out
    


##Smart optimization.
def rhoTest(best, phenotype, cov_matrix, Sigma_qs, mixed, lastMove, rhoArray, verbose = False):
    if best is None:
        # initialize rho opitimization
        best = {}
        best["rho_mid"] = int(np.ceil(len(rhoArray)/2))
        best["rho_left"] = best["rho_mid"]-1
        best["rho_right"] = best["rho_mid"]+1
        lmm_left = LMM(phenotype, cov_matrix, Sigma_qs[rhoArray[best["rho_left"]]])
        lmm_mid = LMM(phenotype, cov_matrix, Sigma_qs[rhoArray[best["rho_mid"]]])
        lmm_right = LMM(phenotype, cov_matrix, Sigma_qs[rhoArray[best["rho_right"]]])
        if not mixed:
            lmm_left.delta = 1
            lmm_left.fix('delta')
            lmm_mid.delta = 1
            lmm_mid.fix('delta')
            lmm_right.delta = 1
            lmm_right.fix('delta')
        lmm_left.fit(verbose=False)
        lmm_mid.fit(verbose=False)
        lmm_right.fit(verbose=False)
        best["lmm_left"] = lmm_left
        best["lml_left"] = lmm_left.lml()
        best["lmm_mid"] = lmm_mid
        best["lml_mid"] = lmm_mid.lml()
        best["lmm_right"] = lmm_right
        best["lml_right"] = lmm_right.lml()
        if verbose:
            print(rhoArray[best["rho_left"]])
            print(best["lml_left"])
            print(rhoArray[best["rho_mid"]])
            print(best["lml_mid"])
            print(rhoArray[best["rho_right"]])
            print(best["lml_right"])
    move=None
    ##Bigger is better lml
    if(best["rho_left"] == 0 or best["rho_right"] == (len(rhoArray)-1)):
        ##Check what is the best of the three.
        return(returnBestRho(best,rhoArray))
    elif best["lml_left"] >= best["lml_mid"] and best["lml_left"] > best["lml_right"] :
        #Need to do a step to the left.
        move = "left"
        ##Recycle old fits.
        best["rho_right"] = best["rho_mid"]
        best["lmm_right"] = best["lmm_mid"]
        best["lml_right"] = best["lml_mid"]
        best["rho_mid"] = best["rho_left"]
        best["lml_mid"] = best["lml_left"]
        best["lmm_mid"] = best["lmm_left"]
        
        ##Compute new.
        best["rho_left"] = best["rho_left"]-1
        lmm_left = LMM(phenotype, cov_matrix, Sigma_qs[rhoArray[best["rho_left"]]])
        if not mixed:
            lmm_left.delta = 1
            lmm_left.fix('delta')
        lmm_left.fit(verbose=False)
        best["lmm_left"] = lmm_left
        best["lml_left"] = lmm_left.lml()
        
    elif best["lml_right"] >= best["lml_mid"] and  best["lml_right"] > best["lml_left"] :
        #Need to do a step to the right.
        move = "right"
        ##Recycle.
        best["rho_left"] = best["rho_mid"]
        best["lmm_left"] = best["lmm_mid"]
        best["lml_left"] = best["lml_mid"]
        best["rho_mid"] = best["rho_right"]
        best["lmm_mid"] = best["lmm_right"]
        best["lml_mid"] = best["lml_right"]
        
        ##Compute new.
        best["rho_right"] = best["rho_right"]+1
        lmm_right = LMM(phenotype, cov_matrix, Sigma_qs[rhoArray[best["rho_right"]]])
        if not mixed:
            lmm_right.delta = 1
            lmm_right.fix('delta')
        lmm_right.fit(verbose=False)
        best["lmm_right"] = lmm_right
        best["lml_right"] = lmm_right.lml()
    #Here we start dealing with the special cases.
    elif ((best["lml_left"] == best["lml_mid"] and best["lml_mid"] == best["lml_right"]) or (best["lml_left"] == best["lml_right"] and best["lml_mid"] < best["lml_right"])):
        ##We are not at an extreme yet. And might nog be done optimizing.
        ##We are either at a flat spot (condition 1)
        ##Or we are symetrically surounding a local minima (moving in either directions is better (and the same))
        
        ## Solution fit extra values around it, to see if we can get away from the flatspot / minima.
        # Fit extra on the right side.
        lmm_right = LMM(phenotype, cov_matrix, Sigma_qs[rhoArray[best["rho_right"]+1]])
        if not mixed:
            lmm_right.delta = 1
            lmm_right.fix('delta')
        lmm_right.fit(verbose=False)
        #Fit one to the left
        lmm_left = LMM(phenotype, cov_matrix, Sigma_qs[rhoArray[best["rho_left"]-1]])
        if not mixed:
            lmm_left.delta = 1
            lmm_left.fix('delta')
        lmm_left.fit(verbose=False)
        #Determine which way to step.
        
        if (lmm_left.lml() > lmm_right.lml()) :
            #Left is the way to go. Furhter left is larger than mid and larger than further right.
            move = "left"
            #moving to the left.
            best["rho_right"] = best["rho_mid"]
            best["lmm_right"] = best["lmm_mid"]
            best["lml_right"] = best["lml_mid"]
            best["rho_mid"] = best["rho_left"]
            best["lmm_mid"] = best["lmm_left"]
            best["lml_mid"] = best["lml_left"]
            
            ##Geting already computed new values.
            best["rho_left"] = best["rho_left"]-1
            best["lmm_left"] = lmm_left
            best["lml_left"] = lmm_left.lml()
        elif (lmm_right.lml() > lmm_left.lml()) :
            #Right is the way to go. Furhter right is larger than mid and larger than further left.
            move = "right"
            #moving to the right.
            best["rho_left"] = best["rho_mid"]
            best["lmm_left"] = best["lmm_mid"]
            best["lml_left"] = best["lml_mid"]
            best["rho_mid"] = best["rho_right"]
            best["lmm_mid"] = best["lmm_right"]
            best["lml_mid"] = best["lml_right"]
            ##Geting already computed new values.
            best["rho_right"] = best["rho_right"]+1
            best["lmm_right"] = lmm_right
            best["lml_right"] = lmm_right.lml()
        else:
            ##Mid might be a minima or plateau and both move in the same lml value or stay flat.
            ##Defaulting back to testing all values.
            return(rhoTestBF(best, phenotype, cov_matrix, Sigma_qs, mixed, lastMove, rhoArray, verbose))
    else:
        ##We are already surrounding an optimum. (might be local).
        ##Check what is the best of the three.
        return(returnBestRho(best, rhoArray))
    
    if(move=="left" and lastMove=="right") or (move=="right" and lastMove=="left"):
        ##Moving back and forth, need to check global optimal (in range). 
        return(rhoTestBF(best, phenotype, cov_matrix, Sigma_qs, mixed, lastMove, rhoArray, verbose))
        
    if verbose:
        print(rhoArray[best["rho_left"]])
        print(best["lml_left"])
        print(rhoArray[best["rho_mid"]])
        print(best["lml_mid"])
        print(rhoArray[best["rho_right"]])
        print(best["lml_right"])
        print(move)
        print()
    
    return rhoTest(best, phenotype, cov_matrix, Sigma_qs, mixed, move, rhoArray, verbose)

##rho brute force.
def rhoTestBF(best, phenotype, cov_matrix, Sigma_qs, mixed, lastMove, rhoArray, verbose = False):
    maxLml = -sys.float_info.max
    mixingParameters = {}
    bestLmm = None
    loc = []
    posBuffer = 0
    
    for i in rhoArray:
        lmm = LMM(phenotype, cov_matrix, Sigma_qs[i])
        if not mixed:
            lmm.delta = 1
            lmm.fix('delta')
        lmm.fit(verbose=False)
        if verbose:
            print(i)
            print(lmm.lml())
        
        if(lmm.lml()>maxLml):
            maxLml = lmm.lml()
            loc = [posBuffer]
            bestLmm = lmm
        elif(lmm.lml()==maxLml):
            loc.append(posBuffer)
        posBuffer +=1
    
    if(len(loc)==1):
        #We have one optimum.
        if verbose:
                print(rhoArray[loc[0]])
        mixingParameters["rho"] = rhoArray[loc[0]]
        mixingParameters["lmm"] =  bestLmm
    else :
        #pdb.set_trace()
        #There are multiple best values, picking the one closest to the mid-point.
        minPoint = (np.ceil(len(rhoArray)/2)-1)
        locArray = np.abs(np.asarray(loc)-4)
        i = rhoArray[loc[np.random.choice(np.where(locArray==min(locArray))[0])]]
        lmm = LMM(phenotype, cov_matrix, Sigma_qs[i])
        lmm.fit(verbose=False)
        if verbose:
            print("Picked value closest to mid point.")
            print(i)
        mixingParameters["rho"] = i
        mixingParameters["lmm"] =  lmm
    return(mixingParameters)

def returnBestRho(best, rhoArray):
    mixingParameters = {}
    if (best["lml_left"] > best["lml_mid"] and best["lml_left"] > best["lml_right"]):
        #Left is best.
        mixingParameters["rho"] = rhoArray[best["rho_left"]]
        mixingParameters["lmm"] =  best["lmm_left"]
    elif (best["lml_left"] < best["lml_right"] and best["lml_mid"] < best["lml_right"]):
        #Right is best.
        mixingParameters["rho"] = rhoArray[best["rho_right"]]
        mixingParameters["lmm"] = best["lmm_right"]
    elif (best["lml_left"] <= best["lml_mid"] and best["lml_mid"] >= best["lml_right"]):
        #midle is best, or all are equal.
        mixingParameters["rho"] = rhoArray[best["rho_mid"]]
        mixingParameters["lmm"] = best["lmm_mid"]
    elif (best["lml_left"] == best["lml_right"] and best["lml_mid"] < best["lml_right"] and best["lml_mid"] < best["lml_left"]):
        #This is the unlikely case that the middle is a minima and steps away from this behave the same.
        print("Stuck in minima, and can't get out of it easily. ")
        print("(Continuning with the mid point: "+str(best["rho_mid"])+")")
        mixingParameters["rho"] = rhoArray[best["rho_mid"]]
        mixingParameters["lmm"] = best["lmm_mid"]
    else:
        print(rhoArray[best["rho_left"]])
        print(best["lml_left"])
        print("")
        print(rhoArray[best["rho_mid"]])
        print(best["lml_mid"])
        print("")
        print(rhoArray[best["rho_right"]])
        print(best["lml_right"])
        print("")
        print("Broken logic")
    return(mixingParameters)

def economic_qs(K, epsilon=np.sqrt(np.finfo(float).eps)):
    r"""Economic eigen decomposition for symmetric matrices.
    A symmetric matrix ``K`` can be decomposed in
    :math:`\mathrm Q_0 \mathrm S_0 \mathrm Q_0^\intercal + \mathrm Q_1\
    \mathrm S_1 \mathrm Q_1^ \intercal`, where :math:`\mathrm S_1` is a zero
    matrix with size determined by ``K``'s rank deficiency.
    Args:
        K (array_like): Symmetric matrix.
        epsilon (float): Eigen value threshold. Default is
                         ``sqrt(finfo(float).eps)``.
    Returns:
        tuple: ``((Q0, Q1), S0)``.
    """
    
    (S, Q) = eigh(K)
    
    nok = abs(max(Q[0].min(), Q[0].max(), key=abs)) < epsilon
    nok = nok and abs(max(K.min(), K.max(), key=abs)) >= epsilon
    if nok:
        from scipy.linalg import eigh as sp_eigh
        (S, Q) = sp_eigh(K)
    
    ok = S >= epsilon
    nok = np.logical_not(ok)
    S0 = S[ok]
    Q0 = Q[:, ok]
    Q1 = Q[:, nok]
    return ((Q0, Q1), S0)

def economic_qs_linear(G, return_q1=True):
    """
    Economic eigen decomposition for a symmetric matrix 𝙺=𝙶𝙶ᵀ.
    Let us define ::
        𝙺 = [𝚀₀  𝚀₁] [𝚂₀  𝟎] [𝚀₀ᵀ]
                     [ 𝟎  𝟎] [𝚀₁ᵀ]
    where the eigenvectors are the columns of [𝚀₀  𝚀₁] and the positive
    eigenvalues are the diagonal elements of 𝚂₀.
    Args:
        G (array_like): Matrix.
        return_q1 (bool): Return 𝚀₁ matrix. Defaults to ``True``.
    Returns:
        tuple: ((𝚀₀, 𝚀₁), 𝚂₀).
    """
    
    if not isinstance(G, da.Array):
        G = np.asarray(G, float)
    
    if not return_q1:
        return _economic_qs_linear_noq1(G)
    
    if G.shape[0] > G.shape[1]:
        (Q, Ssq, _) = svd(G, full_matrices=True)
        S0 = Ssq ** 2
        rank = len(S0)
        Q0, Q1 = Q[:, :rank], Q[:, rank:]
        return ((Q0, Q1), S0)

    return economic_qs(G.dot(G.T))

def _economic_qs_linear_noq1(G):
    if G.shape[0] > G.shape[1]:
        (Q0, Ssq, _) = svd(G, full_matrices=False)
        S0 = Ssq ** 2
        return ((Q0,), S0)
    QS = economic_qs(G.dot(G.T))
    return ((QS[0][0],), QS[1])
    
def lrt_pvalues(null_lml, alt_lmls, dof=1):
    """
    Compute p-values from likelihood ratios.
    These are likelihood ratio test p-values.
    Parameters
    ----------
    null_lml : float
        Log of the marginal likelihood under the null hypothesis.
    alt_lmls : array_like
        Log of the marginal likelihoods under the alternative hypotheses.
    dof : int
        Degrees of freedom.
    Returns
    -------
    pvalues : ndarray
        P-values.
    """
    
    from numpy import clip, inf
    super_tiny = np.finfo(float).tiny
    tiny = np.finfo(float).eps

    lrs = clip(-2 * null_lml + 2 * np.asarray(alt_lmls, float), super_tiny, inf)
    pv = chi2(df=dof).sf(lrs)
    return clip(pv, super_tiny, 1 - tiny)

def ddot(L, R, left=None, out=None):
    r"""Dot product of a matrix and a diagonal one.
    Args:
        L (array_like): Left matrix.
        R (array_like): Right matrix.
        out (:class:`numpy.ndarray`, optional): copy result to.
    Returns:
        :class:`numpy.ndarray`: Resulting matrix.
    """
    L = np.asarray(L, float)
    R = np.asarray(R, float)
    if left is None:
        ok = min(L.ndim, R.ndim) == 1 and max(L.ndim, R.ndim) == 2
        if not ok:
            msg = "Wrong array layout. One array should have"
            msg += " ndim=1 and the other one ndim=2."
            raise ValueError(msg)
        left = L.ndim == 1
    if left:
        if out is None:
            out = np.copy(R)
        L = L.reshape(list(L.shape) + [1] * (R.ndim - 1))
        return np.multiply(L, R, out=out)
    else:
        if out is None:
            out = np.copy(L)
        return np.multiply(L, R, out=out)

def sum2diag(A, D, out=None):
    r"""Add values ``D`` to the diagonal of matrix ``A``.
    Args:
        A (array_like): Left-hand side.
        D (array_like or float): Values to add.
        out (:class:`numpy.ndarray`, optional): copy result to.
    Returns:
        :class:`numpy.ndarray`: Resulting matrix.
    """
    A = np.asarray(A, float)
    D = np.asarray(D, float)
    if out is None:
        out = np.copy(A)
    else:
        np.copyto(out, A)
    np.einsum("ii->i", out)[:] += D
    return out

def glmm_posteriori_covariance_safe_decomposition(objct, feature):
    """Covariance of the estimated posteriori."""
    # We want to compute the posterior covariance:
    #   (K⁻¹ + T⁻¹)⁻¹ = (QS⁻¹Qᵗ + T⁻¹)⁻¹ = T - TQ(S + QᵗTQ)⁻¹QᵗT
    
    K = GLMM.covariance(objct)
    S, Q = eigh(K)
    T = objct._ep._posterior.tau
    
    M = Q.T @ ddot(T, Q, left=True)
    L = ddot(T, Q, left=True)
    R = L.T
    posterior_covariance = -sum2diag(L @ solve(sum2diag(M, S), R), -T)
    
    return posterior_covariance
    
##OLD
#def glmm_posteriori_covariance_safe_decomposition(objct, feature):
#    # We want to compute the posterior covariance:
#    # 
#    #   (K⁻¹ + T⁻¹)⁻¹ = (QS⁻¹Qᵗ + T⁻¹)⁻¹ = T - TQ(S + QᵗTQ)⁻¹QᵗT
#    
#    r"""Covariance of the estimated posteriori."""
#    
#    K = GLMM.covariance(objct)
#    K = K + (1e-3 * np.identity(K.shape[0]))
#    
#    tau = objct._ep._posterior.tau
#    tau = tau + 1e-3
#    
#    if not np.all(np.isfinite(1 / tau)):
#        raise Exception("Tau is way to small")
#        pdb.set_trace();
#    if not np.all(np.isfinite(K)):
#        raise Exception("What? K is weird man...");
#        pdb.set_trace();
#    
#    np.savetxt("./K."+feature+".txt",K)
#    np.savetxt("./tau."+feature+".txt",tau)
#    
#    return pinv(pinv(K + np.diag(1 / tau)))
