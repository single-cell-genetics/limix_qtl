import sys
import os
sys.path.append(os.getcwd())
import h5py
import glob
import os.path
import numpy as np
import pandas as pd
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import dask.array
import glob
from numpy_sugar.linalg import economic_qs, economic_svd
from glimix_core.lmm import LMM
import qtl_utilities as utils
import qtl_loader_utils



def qtl_plots(row,phenotype_ds,phenotype_corrected_ds,plinkGenotype,
              bim,fam,bed,annotation_df,sample2individual_df):
    # extracting SNP index
    print("Processing SNP: "+row[0])
    snp_idx = bim[bim.index == row[0]]["i"]

    # extracting snp-genotypes association
    if (plinkGenotype):
        snp_df = pd.DataFrame(data=bed[snp_idx.values, :].compute().transpose(), index=fam.index, columns=[row[0]])
    else:
        snp_df_dosage = pd.DataFrame(np.nan, index=fam.index, columns=row.index)
        snp_df = pd.DataFrame(np.nan, index=fam.index, columns=snp_names)
        geno = bgen["genotype"][snpId].compute()
        if (geno["ploidy"].min() > 1 & geno["ploidy"].max() < 3):
            if (geno["phased"]):
                snp_df_dosage_t = geno["probs"][:, [0, 2]].sum(1).astype(float)
                snp_df_t = (np.abs(np.argmax(geno["probs"][:, :2], axis=1) - 1) + np.abs(
                    np.argmax(geno["probs"][:, 2:4], axis=1) - 1)).astype(float)
                naId = (np.amax(geno["probs"][:, :2], 1) + np.amax(geno["probs"][:, 2:4], 1)) < (
                            1 + minimumProbabilityStep)
                snp_df_dosage_t[naId] = float('NaN')
                snp_df_t[naId] = float('NaN')
            else:
                snp_df_dosage_t = ((geno["probs"][:, 0] * 2) + geno["probs"][:, 1]).astype(float)
                snp_df_t = (np.abs(np.argmax(geno["probs"][:, :3], axis=1) - 2)).astype(float)
                naId = np.amax(geno["probs"][:, :3], 1) < ((1 / 3) + minimumProbabilityStep)
                snp_df_dosage_t[naId] = float('NaN')
                snp_df_t[naId] = float('NaN')
            snp_df_dosage.loc[:, snp_names[rowNumber]] = snp_df_dosage_t
            snp_df.loc[:, snp_names[rowNumber]] = snp_df_t
        snp_df_dosage = snp_df_dosage.loc[individual_ids, :]

    # filter genotypes
    snp_df = snp_df[snp_df.index.isin(sample2individual_df["iid"].values)]
    snp_pheno_df = pd.merge(snp_df.reset_index(), sample2individual_df)
    snp_pheno_df = snp_pheno_df.merge(phenotype_corrected_ds,left_on="sample",right_index=True)

    # convert to variant
    variant_dic = dict([(0,str(bim[bim.index == row[0]]["a0"][0])+str(bim[bim.index == row[0]]["a0"][0])),
                    (1,str(bim[bim.index == row[0]]["a0"][0])+str(bim[bim.index == row[0]]["a1"][0])),
                    (2,str(bim[bim.index == row[0]]["a1"][0])+str(bim[bim.index == row[0]]["a1"][0]))])

    snp_pheno_df = snp_pheno_df.replace({snp_pheno_df.columns[1]:variant_dic})
    snp_pheno_df = snp_pheno_df.sort_values(by=snp_pheno_df.columns[1])
    if not os.path.exists(output_directory + "/plots/"):
        os.makedirs(output_directory + "/plots/")
    # plot phenotype corrected
    plt.figure()
    sns.set(style="whitegrid")
    ax = sns.violinplot(x=snp_pheno_df.columns[1], y="exp",
                   data=snp_pheno_df, alpha =.1, inner=None)
    plt.setp(ax.collections,alpha=.3)
    ax = sns.swarmplot(x= snp_pheno_df.columns[1], y="exp",
                       data = snp_pheno_df)
    ax.set(ylabel = "Covariates adjusted expression")
    figure = ax.get_figure()
    figure.savefig(output_directory + "/plots/"+snp_pheno_df.columns[1]+"_violinplot_corrected.pdf", dpi=400)

    # filter genotypes
    snp_df = snp_df[snp_df.index.isin(sample2individual_df["iid"].values)]
    snp_pheno_df = pd.merge(snp_df.reset_index(), sample2individual_df)
    snp_pheno_df = snp_pheno_df.merge(phenotype_ds,left_on="sample",right_index=True)

    # convert to variant
    variant_dic = dict([(0,str(bim[bim.index == row[0]]["a0"][0])+str(bim[bim.index == row[0]]["a0"][0])),
                    (1,str(bim[bim.index == row[0]]["a0"][0])+str(bim[bim.index == row[0]]["a1"][0])),
                    (2,str(bim[bim.index == row[0]]["a1"][0])+str(bim[bim.index == row[0]]["a1"][0]))])

    snp_pheno_df = snp_pheno_df.replace({snp_pheno_df.columns[1]:variant_dic})
    snp_pheno_df = snp_pheno_df.sort_values(by=snp_pheno_df.columns[1])
    # plot phenotype standard
    plt.figure()
    sns.set(style="whitegrid")
    ax = sns.violinplot(x=snp_pheno_df.columns[1], y="exp",
                   data=snp_pheno_df, alpha =.1, inner=None)
    plt.setp(ax.collections,alpha=.3)
    ax = sns.swarmplot(x= snp_pheno_df.columns[1], y="exp",
                       data = snp_pheno_df)
    ax.set(ylabel = "Covariates adjusted expression")
    figure = ax.get_figure()
    figure.savefig(output_directory + "/plots/"+snp_pheno_df.columns[1]+"_violinplot.pdf", dpi=400)


def run_plots(plinkGenotype, geno_prefix, annotation_filename, phenotype_filename,
              covariate_filename, top_qtl_results_filename, output_directory,
              sample_mapping_filename,randomeff_filename):

    # # Loading Files
    # phenotype_filename = "/Users/chaaya/dhonveli_dkfz/limix_qtl/limix_qtl-master/Limix_QTL/test_data/Expression/Geuvadis_CEU_YRI_Expr.txt"
    # annotation_filename = "/Users/chaaya/dhonveli_dkfz/limix_qtl/limix_qtl-master/Limix_QTL/test_data/Expression/Geuvadis_CEU_Annot.txt"
    # covariate_filename = "/Users/chaaya/dhonveli_dkfz/limix_qtl/limix_qtl-master/Limix_QTL/test_data/Expression/Geuvadis_CEU_YRI_covariates.txt"
    # sample_mapping_filename = "/Users/chaaya/dhonveli_dkfz/limix_qtl/limix_qtl-master/Limix_QTL/test_data/Geuvadis_CEU_gte.txt"
    # geno_prefix = "/Users/chaaya/dhonveli_dkfz/limix_qtl/limix_qtl-master/Limix_QTL/test_data/Genotypes/Geuvadis"
    # output_directory = "/Users/chaaya/dhonveli_dkfz/limix_qtl/limix_qtl-master/Limix_QTL/test_data/Output"
    # top_qtl_results_filename = "/Users/chaaya/dhonveli_dkfz/limix_qtl/limix_qtl-master/Limix_QTL/test_data/Output/top_qtl_results_all_FDR0.05.txt"
    # plinkGenotype = True


    [phenotype_df, kinship_df, readdepth_df, covariate_df, sample2individual_df,complete_annotation_df, annotation_df, snp_filter_df,
     snp_feature_filter_df, geneticaly_unique_individuals, minimum_test_samples, feature_list, bim, fam, bed, bgen,
     chromosome, selectionStart, selectionEnd, feature_variant_covariate_df]=\
    utils.run_QTL_analysis_load_intersect_phenotype_covariates_kinship_sample_mapping(pheno_filename=phenotype_filename,
                    anno_filename=annotation_filename, geno_prefix=geno_prefix, plinkGenotype=plinkGenotype, cis_mode=True,
                    skipAutosomeFiltering = False, minimum_test_samples= 10,
                    relatedness_score=None, snps_filename=None, feature_filename=None,
                    snp_feature_filename=None, selection='all',covariates_filename=covariate_filename,
                    randomeff_filename=randomeff_filename, sample_mapping_filename=sample_mapping_filename,
                    extended_anno_filename=None, feature_variant_covariate_filename=None)

    # results
    top_qtl_results_df = qtl_loader_utils.get_top_qtl_results(top_qtl_results_filename)

    for row in top_qtl_results_df.iterrows():

        # feature specific parameters for QS mixing
        rho1 = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
        Sigma = {}
        Sigma_qs = {}
        best = {}
        snpQcInfo = None
        currentFeatureNumber+= 1

        phenotype_ds = phenotype_df.loc[row[1]["feature_id"]]
        individual_ids = sample2individual_df.loc[phenotype_ds.index, 'iid'].values
        sample2individual_feature = sample2individual_df.loc[phenotype_ds.index]

        if kinship_df is not None and readdepth_df is None:
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

        randomeff_mix = False
        # combining the two matrices
        if kinship_df is not None and readdepth_df is not None:
            randomeff_mix = True
            kinship_mat = kinship_df.loc[individual_ids,individual_ids].values
            kinship_mat = kinship_mat.astype(float)
            ##GOWER normalization of Kinship matrix.
            kinship_mat *= (kinship_mat.shape[0] - 1) / (kinship_mat.trace() - kinship_mat.mean(0).sum())
            for rho in rho1:
                Sigma[rho] = rho * kinship_df + (1 - rho) * readdepth_df
                Sigma_qs[rho] = economic_qs(Sigma[rho])

        # creating a fake QS if none random effect is present or use the read depth one
        if kinship_df is None:
            if readdepth_df is None:
                K = np.eye(len(phenotype_ds.index))
                if(QS is None and not contains_missing_samples):
                    QS = economic_qs(K)
                elif (contains_missing_samples):
                    QS_tmp = QS
                    QS = economic_qs(K)
            else:
                if(QS is None and not contains_missing_samples):
                    QS = economic_qs(readdepth_df)
                elif (contains_missing_samples):
                    QS_tmp = QS
                    QS = economic_qs(readdepth_df)

        # covariance matrix
        cov_matrix = covariate_df.loc[sample2individual_feature['sample'], :].values if covariate_df is not None else None
        if covariate_df is None:
            cov_matrix = np.ones((len(individual_ids), 1))
        cov_matrix = cov_matrix.astype(float)

        phenotype = utils.force_normal_distribution(phenotype_ds.values,method="gaussnorm")

        # Prepare LMM
        phenotype = phenotype.astype(float)

        if randomeff_mix:
            # initialize best to minus infinite
            best["lml"] = - math.inf
            best["lmm"] = - math.inf
            best["rho1"] = - math.inf
            for rho, QS in Sigma_qs.items():
                lmm = LMM(phenotype, cov_matrix, QS)
                if not mixed:
                    lmm.delta = 1
                    lmm.fix('delta')
                lmm.fit(verbose=False)
                lml = lmm.lml()
                if lml > best["lml"]:
                    best["lml"] = lml
                    best["lmm"] = lmm
                    best["rho1"] = rho
            lmm = best["lmm"]
            print(best["rho1"])
            if best["rho1"] != 0:
                rho_log[(feature_id)] = best["rho1"]
                print("Read depth has actually an effect!")

        else:
            lmm = LMM(phenotype, cov_matrix, QS)
            if not mixed:
                lmm.delta = 1
                lmm.fix('delta')
            lmm.fit(verbose=False)

        # create phenotype_corrected_df for plotting
        phenotype_corrected = phenotype - cov_matrix[:, 1:].dot(lmm.beta[1:])
        phenotype_corrected_ds = pd.Series(data = phenotype_corrected, index = phenotype_ds.index, name="exp")

        qtl_plots(row, phenotype_ds, phenotype_corrected_ds, plinkGenotype, bim, fam, bed, annotation_df, sample2individual_df)


def parse_args():
    parser = argparse.ArgumentParser(description='Run standard QTL plots after minimal postprocess was run')
    parser.add_argument('--bgen','-bg',required=False)
    parser.add_argument('--plink','-pg',required=False)
    parser.add_argument('--annotation_file','-af', required=True)
    parser.add_argument('--phenotype_file','-pf', required=True)
    parser.add_argument('--randomeff_file', '-rf', required=False)
    parser.add_argument("--top_qtl_results","-tqr",required=True)
    parser.add_argument("--output_directory","-od",required=True)
    parser.add_argument('--covariates_file','-cf', required=False)
    parser.add_argument("--gaussianize_method","-gm",required=False)
    parser.add_argument("--sample_mapping_file","-smf",required=False)
    args = parser.parse_args()
    return args

if __name__=='__main__':
    args = parse_args()
    bgen  = args.bgen
    plink = args.plink
    annotation_filename = args.annotation_file
    phenotype_filename = args.phenotype_file
    top_qtl_results_filename = args.top_qtl_results
    output_directory  = args.output_directory
    covariates_filename = args.covariates_file
    sample_mapping_filename = args.sample_mapping_file
    randomeff_filename = args.randomeff_file

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

    run_plots(plinkGenotype, geno_prefix, annotation_filename,
              phenotype_filename, covariates_filename, top_qtl_results_filename,
              output_directory, sample_mapping_filename, randomeff_filename)
