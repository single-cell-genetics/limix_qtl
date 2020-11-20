import glob
import os.path
import pandas as pd
import argparse
import numpy as np
import scipy.stats
from scipy.stats import beta

#V0.1.2
BETA_SHAPE1_MIN = 0.1
BETA_SHAPE1_MAX = 10
BETA_SHAPE2_MIN_CIS = 1
BETA_SHAPE2_MIN_TRANS = 5
BETA_SHAPE2_MAX_CIS = 1000000
BETA_SHAPE2_MAX_TRANS = 100000000

def estimate_beta_function_paras(top_pvalues_perm):
    mean = np.mean(top_pvalues_perm)
    variance = np.var(top_pvalues_perm)
    alpha_para = mean * (mean * (1 - mean ) / variance - 1)
    beta_para = alpha_para * (1 / mean - 1)
    return alpha_para,beta_para

def correction_function_fdr(pValue, top_pvalues_perm, nPerm):
    fdrPval = len(np.where(top_pvalues_perm<pValue)[0])/nPerm
    if(fdrPval>1.0):
        fdrPval =1.0
    return fdrPval

def add_global_fdr_measures(QTL_Dir, OutputDir, relevantGenes, qtl_results_file="top_qtl_results_all.txt"):
    if QTL_Dir[-1:] == "/" :
        QTL_Dir = QTL_Dir[:-1]
    if OutputDir[-1:] == "/" :
        OutputDir = OutputDir[:-1]
    
    if relevantGenes is not None :
        genesToParse = pd.read_csv(relevantGenes, header=None)[0].values
        toRead = set(QTL_Dir+"/Permutation.pValues."+genesToParse+".txt")
    
    permutationInformtionToProcess = (glob.glob(QTL_Dir+"/Permutation.pValues.*.txt"))
    
    if relevantGenes is not None :
        permutationInformtionToProcess = set(permutationInformtionToProcess).intersection(toRead)
    pValueBuffer = []
    genesTested = 0
    for file in permutationInformtionToProcess :
        #print(file)
        pValueBuffer.extend(np.loadtxt(file))
        genesTested +=1
    nPerm = len(pValueBuffer)/genesTested
    
    pValueBuffer=np.float_(pValueBuffer)
    alpha_para, beta_para = estimate_beta_function_paras(pValueBuffer)
    beta_dist_mm = scipy.stats.beta(alpha_para,beta_para)
    correction_function_beta = lambda x: beta_dist_mm.cdf(x)
    
    qtlResults = pd.read_table(QTL_Dir+"/"+qtl_results_file,sep='\t')
    
    if relevantGenes is not None :
        qtlResults = qtlResults.loc[qtlResults['feature_id'].isin(genesToParse)]
    
    qtlResults['empirical_global_p_value'] = correction_function_beta(qtlResults["p_value"])
    
    fdrBuffer = []
    for p in qtlResults["p_value"] :
        fdrBuffer.append(correction_function_fdr(p , pValueBuffer, nPerm))
    qtlResults['emperical_global_fdr'] = fdrBuffer
    qtlResults.to_csv(QTL_Dir+"/top_qtl_results_all_global_FDR_info.txt",sep='\t',index=False)

def parse_args():
    parser = argparse.ArgumentParser(description='Run global gene and snp wide correction on QTLs.')
    parser.add_argument('--input_dir','-id',required=True)
    parser.add_argument('--ouput_dir','-od',required=True)
    parser.add_argument('--gene_selection','-gs',required=False, default=None)
    parser.add_argument('--qtl_filename','-qf',required=False, default=None)
    args = parser.parse_args()
    return args

if __name__=='__main__':
    args = parse_args()
    inputDir  = args.input_dir
    outputDir = args.ouput_dir
    relevantGenes = args.gene_selection
    qtlFileName = args.qtl_filename
    
    if qtlFileName is not None :
        add_global_fdr_measures(inputDir, outputDir, relevantGenes, qtlFileName)
    else :
        add_global_fdr_measures(inputDir, outputDir, relevantGenes)



