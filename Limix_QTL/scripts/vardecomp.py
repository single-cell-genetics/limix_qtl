import sys
import os
import argparse
import cPickle
import numpy as np
import pandas as pd
import scipy as sp
import limix
from limix.varDecomp import VarianceDecomposition
# import limix.modules.varianceDecomposition as var
import limix.modules.qtl as qtl
import limix.io.data as data
import limix.io.genotype_reader as gr
import limix.io.phenotype_reader as phr
import limix.io.data_util as data_util
import limix.utils.preprocess as preprocess
sys.path.append('/nfs/users/nfs_n/nc10/bin/LIMIX-eQTL/pysrc')
from CFG.settings import *
import limix.stats.fdr as fdr
from include.data import QtlData
from include.utils import smartDumpDictHdf5
from include.utils import getLambda
from limix.utils.preprocess import gaussianize
import pylab as pl
import h5py
import copy
import warnings
import time

def match_samples(*sampleIDs):
    sampleID_common = sampleIDs[0]
    for sampleID in sampleIDs: 
        sampleID_common = sp.intersect1d(sampleID, sampleID_common)
    idxs = []
    for sampleID in sampleIDs:
        _idxs = sp.array([sp.where(sampleID==sample)[0][0] for sample in sampleID_common])
        assert (sampleID[_idxs]==sampleID_common).all()
        idxs.append(_idxs)
    return idxs

# def is_pos_def(x):
#     return np.all(np.linalg.eigvals(x) > 0)
def isPSD(A, tol=1e-6):
    E,V = sp.linalg.eigh(A)
    return np.all(E > -tol)

def get_args():

    parser = argparse.ArgumentParser(description='Test HierGWAS on MD genotypes and classifier phenotypes')
    parser.add_argument('-geno_dir','--geno_dir', required=True)
    parser.add_argument('-pheno_file','--pheno_file', required=False)
    parser.add_argument('-mt_file','--mt_file', required=True)
    parser.add_argument('-out_file','--out_file', required=True)
    parser.add_argument('-cis_window_kb','--cis_window_kb', required=True)
    args = parser.parse_args()

    geno_dir = args.geno_dir
    pheno_file = args.pheno_file
    out_file = args.out_file
    cis_window_kb = args.cis_window_kb
    mt_file=args.mt_file

    return geno_dir, pheno_file, out_file, cis_window_kb, mt_file

if __name__=='__main__':

    geno_dir, pheno_file, out_file, cis_window_kb, mt_file= get_args()

    cis_window = 1000 * cis_window_kb

    MT = gr.genotype_reader_tables(mt_file)
    Mpos = MT.getPos()
    Msnps = MT.getGenotypes()
    Msnps-= Msnps.mean(0)
    Msnps/= Msnps.std(0)
    M = sp.dot(Msnps, Msnps.T)    
    Mids=MT.sample_ID

    # MTrare = gr.genotype_reader_tables(mt_rare_file)
    # MTrarepos = MTrare.getPos()
    # MTraresnps = MTrare.getGenotypes()
    # MTrareids=MTrare.sample_ID

    chr=1
    chrgeno=gr.genotype_reader_tables(geno_dir+"/chrom"+str(chr)+".h5")
    chrsnps = chrgeno.getGenotypes()
    Kchr = chrsnps
    Kchr-= chrsnps.mean(0)
    Kchr/= Kchr.std(0)
    Kchr = sp.dot(Kchr, Kchr.T)
    Kids= chrgeno.sample_ID

    chromosomes=range(2,23)
    for chr in chromosomes:
        chrgeno=gr.genotype_reader_tables(geno_dir+"/chrom"+str(chr)+".h5")
        chrsnps = chrgeno.getGenotypes()
        kchr = chrsnps
        kchr-= chrsnps.mean(0)
        kchr/= kchr.std(0)
        kchr = sp.dot(kchr, kchr.T)
        Kchr+= kchr 

    # read data and split dataset into jobs
    data = QtlData(geno_dir, pheno_file)
    all_ps = data.pID.copy() 
    n_ps = all_ps.shape[0]
    Icv = sp.floor(1 * sp.arange(n_ps) / n_ps)
    ps = all_ps[Icv==0]

    # match samples
    idxs_p, idxs_g, idxs_m, idxs_k = match_samples(data.sampleID_p, data.sampleID_g, Mids, Kids)

    # # import rrm
    # frrm = h5py.File(rrm_file, 'r')
    # K = frrm['RRM'][:][idxs_g][:, idxs_g]
    # frrm.close()
    Kall = Kchr[:][idxs_k][:, idxs_k]

    Kallstd=np.array(Kall)
    Kallstd/= Kallstd.diagonal().mean()
    Kallstd+=1e-4*sp.eye(Kallstd.shape[0])
    Kids=[idxs_k]

    M = M[:][idxs_m][:, idxs_m]
    Mids=Mids[idxs_m]
    M/= M.diagonal().mean()
    Km=M
    Km+=1e-4*sp.eye(Km.shape[0])
    Msnps=Msnps[idxs_m]

    with open (out_file + ".vardecomp.txt", "a") as decomp_file:
        with open (out_file + ".vardecomp_perm.txt", "a") as perm_file:
            with open (out_file + ".lmm.txt", "a") as lmm_file:
                for gene in ps:
                    print '.. Analyzing gene %s' % gene
                    t0 = time.time()
                    Y = data.getPhenotypes(gene)
                    Y = gaussianize(Y[idxs_p])
                    p_info = data.getGeneInfo(gene)
                    try:
                        Xc = data.getGenotypes(gene,cis_window=cis_window)
                        Xc = Xc[idxs_g].astype(float)
                        Xc-= Xc.mean(0)
                        Xc/= Xc.std(0)
                        Kc = sp.dot(Xc, Xc.T)
                        Kt = Kall-Kc
                        Kc /= Kc.diagonal().mean()
                        Kt /= Kt.diagonal().mean()
                        Kc+=1e-4*sp.eye(Kc.shape[0])
                        Kt+=1e-4*sp.eye(Kt.shape[0])
                        # vc1 = var.VarianceDecomposition(Y)
                        ## get null
                        vcnull = VarianceDecomposition(Y)
                        vcnull.addFixedEffect()
                        vcnull.addRandomEffect(is_noise=True)
                        vcnull.optimize()
                        varnull=vcnull.getVarianceComps()
                        lmnull=vcnull.getLML()
                        ## get cis trans effects 
                        vc0 = VarianceDecomposition(Y)
                        vc0.addFixedEffect()
                        vc0.addRandomEffect(K=Kallstd)
                        vc0.addRandomEffect(is_noise=True)
                        vc0.optimize()
                        var0=vc0.getVarianceComps()
                        lm0=vcnull.getLML()-vc0.getLML()
                        ## get cis trans effects 
                        vc1 = VarianceDecomposition(Y)
                        vc1.addFixedEffect()
                        vc1.addRandomEffect(K=Kc)
                        vc1.addRandomEffect(K=Kt)
                        vc1.addRandomEffect(is_noise=True)
                        vc1.optimize()
                        var1=vc1.getVarianceComps()
                        lm1=vcnull.getLML()-vc1.getLML()
                        ## get cis trans mt effects 
                        vc2 = VarianceDecomposition(Y)
                        vc2.addFixedEffect()
                        vc2.addRandomEffect(K=Kc)
                        vc2.addRandomEffect(K=Kt)
                        vc2.addRandomEffect(K=Km)
                        vc2.addRandomEffect(is_noise=True)
                        vc2.optimize()
                        var2=vc2.getVarianceComps()
                        lm2=vc1.getLML()-vc2.getLML()

                        var1/=var1.sum(1)[:,sp.newaxis]
                        var2/=var2.sum(1)[:,sp.newaxis]
                        cis, trans, noise = list(var1[0,:])
                        cismt, transmt, mt, noisemt = list(var2[0,:])
                        lms = "\t".join(map(str, [lmnull, lm0, lm1, lm2]))
                        vcs = "\t".join(map(str, [cis, trans, noise, cismt, transmt, mt, noisemt]))
                        oline="\t".join([gene, lms, vcs])
                        decomp_file.write(oline + "\n")

                        for i in range(0,50):
                            idxs_perm = sp.random.permutation(M.shape[0])
                            Mperm = np.array(M[:][idxs_perm][:, idxs_perm])
                            Mperm = Mperm[idxs_perm]
                            Mperm-= Mperm.mean(0)
                            Mperm/= Mperm.std(0)
                            Mperm = sp.dot(Mperm, Mperm.T)   
                            Mperm/= Mperm.diagonal().mean()
                            Kmperm=Mperm
                            Kmperm+=1e-4*sp.eye(Kmperm.shape[0])
                            vcperm = VarianceDecomposition(Y)
                            vcperm.addFixedEffect()
                            vcperm.addRandomEffect(K=Kc)
                            vcperm.addRandomEffect(K=Kt)
                            vcperm.addRandomEffect(K=Mperm)
                            vcperm.addRandomEffect(is_noise=True)
                            vcperm.optimize()
                            permlm1=vc1.getLML()-vcperm.getLML()
                            Kallperm = np.array(Kchr[:][idxs_perm][:, idxs_perm])
                            Kallperm = Kallperm[idxs_perm]
                            Kallperm-= Kallperm.mean(0)
                            Kallperm/= Kallperm.std(0)
                            Kallperm = sp.dot(Kallperm, Kallperm.T)   
                            Kallperm/= Kallperm.diagonal().mean()
                            Kallperm+=1e-4*sp.eye(Kallperm.shape[0])
                            vcperm = VarianceDecomposition(Y)
                            vcperm.addFixedEffect()
                            vcperm.addRandomEffect(K=Kallperm)
                            vcperm.addRandomEffect(is_noise=True)
                            vcperm.optimize()
                            permlm0=vcnull.getLML()-vcperm.getLML()
                            perm_file.write("\t".join(map(str,[permlm0,permlm1]))+"\n")
                        ## get trans PCs 
                        S_R, U_R = sp.linalg.eigh(Kc)
                        F1 = U_R[:, ::-1][:, :10] 
                        # add an intercept term
                        F1 = sp.concatenate([F1, sp.ones((F1.shape[0], 1))], 1)
                        test="lrt"                  #specify type of statistical test
                        lmm0 = qtl.test_lmm(snps=Msnps,pheno=Y,K=Kallstd,covs=F1,test=test)
                        pvalues = lmm0.getPv()       # 1xS vector of p-values (S=X.shape[1])
                        betas = lmm0.getBetaSNP()    # 1xS vector of effect sizes (S=X.shape[1])
                        ses=lmm0.beta_ste # 1xS vector of effect sizes standard errors (S=X.shape[1]
                        RV=Mpos
                        RV["pvaluesCisPCs"]=pvalues.T
                        RV["betasCisPCs"]=betas.T
                        RV["sesCisPCs"]=ses.T
                        RV["gene"]=gene

                        test="lrt"                  #specify type of statistical test
                        lmm2 = qtl.test_lmm(snps=Msnps,pheno=Y,K=Kallstd,covs=sp.ones((F1.shape[0], 1)),test=test)
                        pvalues = lmm2.getPv()       # 1xS vector of p-values (S=X.shape[1])
                        betas = lmm2.getBetaSNP()    # 1xS vector of effect sizes (S=X.shape[1])
                        ses=lmm2.beta_ste # 1xS vector of effect sizes standard errors (S=X.shape[1]
                        RV["pvaluesK"]=pvalues.T
                        RV["betasK"]=betas.T
                        RV["sesK"]=ses.T

                        RV.to_csv(lmm_file, header=False)
                    except Exception:
                        pass

# wdir="/lustre/scratch115/projects/blood_mt/GTEX/WGS"
# mt_file="/lustre/scratch115/projects/blood_mt/GTEX/WGS/mtplinknew/all.hom.meanmaf5.h5"
# geno_dir="/lustre/scratch115/projects/blood_mt/GTEX/WGS/plink/h5"
# cis_window_kb=1000
# nperm=100
# pheno_file="/lustre/scratch115/projects/blood_mt/GTEX/WGS/pheno/peer/Cells_Transformed_fibroblasts/chr11.h5"

