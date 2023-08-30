from pandas_plink import read_plink
from bgen_reader import read_bgen
import pandas as pd
import numpy as np
import os

#V0.1.1

##loader functions
def ensure_dir(file_path):
    '''Check if directory exists for output, and create it if it doesn't.'''
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)

##Read snp filterlist
def get_snp_df(snps_filename):
    if snps_filename:
        snp_filter_df = pd.read_csv(snps_filename,sep='\t',index_col=0)
    else:
        snp_filter_df = None
    return snp_filter_df

def get_randeff_df(randomeff_filename):
    kinship_filename = False
    readdepth_filename = False
    if ',' in randomeff_filename :
        kinship_filename, readdepth_filename = randomeff_filename.split(",")
    else :
        kinship_filename = randomeff_filename
    if kinship_filename:
        kinship_df = pd.read_csv(kinship_filename,sep='\t',index_col=0)
        kinship_df.index = kinship_df.index.astype(str)
        kinship_df.columns = kinship_df.columns.astype(str)
    else:
        kinship_df = None
    if readdepth_filename:
        readdepth_df = pd.read_csv(readdepth_filename,sep='\t',index_col=0)
        readdepth_df.index = readdepth_df.index.astype(str)
        readdepth_df.columns = readdepth_df.columns.astype(str)
    else:
        readdepth_df = None
    return kinship_df,readdepth_df

def get_samplemapping_df(sample_mapping_filename,sample_labels,key_from):
    assert(key_from in ['iid','sample'])
    if sample_mapping_filename:
        mapping_df = pd.read_csv(sample_mapping_filename,sep='\t',header=None,names=['iid','sample'], dtype={'iid': str, 'sample': str})
        mapping_df.set_index(key_from,inplace=True)
    else:
        #assume the mapping is the identity mapping
        identifiers = sample_labels
        mapping_df = pd.DataFrame(data=np.vstack([identifiers,identifiers]).T,index=identifiers,columns=['iid','sample'], dtype={'iid': str, 'sample': str})
    return mapping_df

def get_snp_feature_df(snp_feature_filename):
    if snp_feature_filename:
        snp_feature_df = pd.read_csv(snp_feature_filename,sep='\t')
    else :
        snp_feature_df = None
    return snp_feature_df

def get_covariate_df(covariates_filename):
    if covariates_filename:
        covariate_df = pd.read_csv(covariates_filename,sep='\t',index_col=0)
        covariate_df.index = covariate_df.index.astype(str)
        covariate_df.columns = covariate_df.columns.astype(str)
    else:
        covariate_df = None
    return covariate_df

def get_genotype_data(geno_prefix, plinkGenotype):
    if(plinkGenotype):
        bim,fam,bed = read_plink(geno_prefix,verbose=False)
        fam.set_index('iid',inplace=True)
        bgen=None
    else :
        bgen = read_bgen(geno_prefix+'.bgen', verbose=False)
        bed=None
        fam =bgen['samples']
        fam = fam.to_frame("iid")
        fam.set_index('iid',inplace=True)
        fam.index = fam.index.astype(str)
        
        bim = bgen['variants'].compute()
        bim = bim.assign(i = range(bim.shape[0]))
        bim['id'] = bim['rsid']
        bim = bim.rename(index=str, columns={"id": "snp"})
        bim['a1'] = bim['allele_ids'].str.split(",", expand=True)[0]
        bim.index = bim["snp"].astype(str).values
        bim.index.name = "candidate"
        
        ##Fix chromosome ids
        bim['chrom'].replace('^chr','',regex=True,inplace=True)
        bim['chrom'].replace(['X', 'Y', 'XY', 'MT'], ['23', '24', '25', '26'],inplace=True)
        ##Remove non-biallelic & non-ploidy 2 (to be sure). (These can't happen in binary plink files).
        print("Warning, the current software only supports biallelic SNPs and ploidy 2")
        bim.loc[np.logical_and(bim['nalleles']<3,bim['nalleles']>0),:]
    
    return bim,fam,bed,bgen

def get_annotation_df(anno_filename):
    annotation_col_dtypes = {'feature_id':str,
                         'gene_id':str,
                         'gene_name':str,
                         'chromosome':str,
                         'start':np.int64,
                         'end':np.int64,
                         'strand':str}
    annotation_df = pd.read_csv(anno_filename,sep='\t',index_col=0,dtype=annotation_col_dtypes)
    return annotation_df

def get_env_df(env_filename):
    env_df = pd.read_csv(env_filename,sep='\t',index_col=0)
    env_df.index = env_df.index.astype(str)
    env_df.columns = env_df.columns.astype(str)
    return env_df
    

def get_phenotype_df(pheno_filename):
    pheno_df = pd.read_csv(pheno_filename,sep='\t',index_col=0, na_values=['.'])
    pheno_df.index = pheno_df.index.astype(str)
    pheno_df.columns = pheno_df.columns.astype(str)
    return pheno_df

def get_grs_subset_df(grs_filename, relSnps):
    iter_csv = pd.read_csv(grs_filename, chunksize=1000, sep='\t',index_col=0, na_values=['.'])
    risk_df = None
    for chunk in iter_csv:
        if any(chunk.index.isin(relSnps)):
            risk_df = pd.concat([risk_df,chunk.reindex(labels=relSnps,axis ='index')])
    return risk_df

def get_top_qtl_results(top_qtl_results_filename):
    if top_qtl_results_filename:
        top_qtl_results_df = pd.read_csv(top_qtl_results_filename,sep='\t',index_col=0)
    else:
        top_qtl_results_df = None
    return top_qtl_results_df
