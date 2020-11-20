import numpy as np
import scipy.stats
from scipy.stats import beta
#from joblib import Parallel

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

def define_correction_function(top_pvalues_perm, cis_mode):
    #Always try to use the MLE estimator, new default to 10 permutations.
    #If the MLE estimator fails we go back to the cruder estimation of the beta distribution.
    offset = (np.finfo(np.double).tiny*100)
    ##Replace zero's value with smallest number not 0.
    top_pvalues_perm[top_pvalues_perm == 0] = offset
    ##Replace highest value with highest number not 1.
    top_pvalues_perm[top_pvalues_perm == 1] = 1-offset
    try :
        alpha_para,beta_para,loc,fscale =  beta.fit(top_pvalues_perm,floc=0,fscale=1)
    except (scipy.stats._continuous_distns.FitSolverError):
        alpha_para,beta_para = estimate_beta_function_paras(top_pvalues_perm)
    except (scipy.stats._continuous_distns.FitDataError):
        alpha_para,beta_para = estimate_beta_function_paras(top_pvalues_perm)
    if(cis_mode):
        if(alpha_para<BETA_SHAPE1_MIN or alpha_para>BETA_SHAPE1_MAX or alpha_para<BETA_SHAPE2_MIN_CIS or alpha_para>BETA_SHAPE2_MAX_CIS):
            alpha_para,beta_para = estimate_beta_function_paras(top_pvalues_perm)
            ### If pvalues become more significant after multiple testing correction we put them back to the orignal test Pvalue in a seperate step.
    else :
        if(alpha_para<BETA_SHAPE1_MIN or alpha_para>BETA_SHAPE1_MAX or alpha_para<BETA_SHAPE2_MIN_TRANS or alpha_para>BETA_SHAPE2_MAX_TRANS):
            alpha_para,beta_para = estimate_beta_function_paras(top_pvalues_perm)
    
    beta_dist = scipy.stats.beta(alpha_para,beta_para)
    correction_function = lambda x: beta_dist.cdf(x)
    #Would be good to replace 0 with minimal double value of python.
    return [correction_function, alpha_para, beta_para]

def calculate_corrected_pvalues(top_pvalues_perm,nominal_pvalues):
    correction_function = define_correction_function(top_pvalues_perm)
    #apply correction to nominal p_values - potentially slow
    corrected_pvalues = np.array([correction_function(x) for x in nominal_pvalues])
    #corrected_pvalues = Parallel(n_jobs=-2)(np.array([correction_function(x) for x in nominal_pvalues]))
    return corrected_pvalues
