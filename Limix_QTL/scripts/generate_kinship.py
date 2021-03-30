import numpy as np
import math
import scipy

def generate_kinship(genotypes):
    kchr = genotypes
    #standardise
    kchr -= kchr.mean(axis=0)
    kchr /= kchr.std(axis=0)
    
    kinship = scipy.dot(kchr, kchr.T)
    return kinship

