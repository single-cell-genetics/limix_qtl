import limix
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

def get_kinship_df(kinship_filename):
    if kinship_filename:
        kinship_df = pd.read_csv(kinship_filename,sep='\t',index_col=0)
    else:
        kinship_df = None
    return kinship_df

def get_samplemapping_df(sample_mapping_filename,sample_labels,key_from):
    assert(key_from in ['iid','sample'])
    if sample_mapping_filename:
        mapping_df = pd.read_csv(sample_mapping_filename,sep='\t',header=None,names=['iid','sample'])
        mapping_df.set_index(key_from,inplace=True)
    else:
        #assume the mapping is the identity mapping
        identifiers = sample_labels
        mapping_df = pd.DataFrame(data=np.vstack([identifiers,identifiers]).T,index=identifiers,columns=['iid','sample'])
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
    else:
        covariate_df = None
    return covariate_df

def get_genotype_data(geno_prefix):
    bim,fam,bed = limix.io.plink.read(geno_prefix,verbose=False)
    fam.set_index('iid',inplace=True)
    return bim,fam,bed

def get_annotation_df(anno_filename):
    annotation_col_dtypes = {'feature_id':np.object,
                         'gene_id':np.object,
                         'gene_name':np.object,
                         'chromosome':np.object,
                         'start':np.int64,
                         'end':np.int64,
                         'strand':np.object}
    annotation_df = pd.read_csv(anno_filename,sep='\t',index_col=0,dtype=annotation_col_dtypes)
    return annotation_df

def get_env_df(env_filename):
    return pd.read_csv(env_filename,sep='\t',index_col=0)

def get_phenotype_df(pheno_filename):
    return pd.read_csv(pheno_filename,sep='\t',index_col=0, na_values=['.'])

def get_grs_subset_df(grs_filename, relSnps):
    iter_csv = pd.read_csv(grs_filename, chunksize=1000, sep='\t',index_col=0, na_values=['.'])
    risk_df = None
    for chunk in iter_csv:
        if any(chunk.index.isin(relSnps)):
            risk_df = pd.concat([risk_df,chunk.reindex(labels=relSnps,axis ='index')])
    return risk_df
