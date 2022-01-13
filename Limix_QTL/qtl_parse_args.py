import argparse

def get_args():
    parser = argparse.ArgumentParser(description='Run QTL analysis given genotype, phenotype, and annotation.')
    parser.add_argument('--bgen','-bg',required=False)
    parser.add_argument('--plink','-pg',required=False)
    parser.add_argument('--annotation_file','-af', required=True)
    parser.add_argument('--phenotype_file','-pf', required=True)
    parser.add_argument('--output_directory','-od', required=True)
    parser.add_argument('--window','-w', required=False,
                        help=
                        'The size of the cis window to take SNPs from.'
                        'The window will extend between:                     '
                        '    (feature_start - (window))             '
                        ' and:                                               '
                        '    (feature_end + (window))               ',default=250000)
    parser.add_argument('--genomic_range','-gr',required=False,
                        help=
                        'A genomic range to do selecte features to be considered in the analysis.'
                        'Available options: all (default), a chromsome or chromosome:start-end.',default='all')
    parser.add_argument('--covariates_file','-cf',required=False,default=None)
    parser.add_argument('--randomeff_files','-rf',required=False,default=None,
                        help = """The random effect files divided by a comma. The first one should be the kinship file and the second one should contain the read depth. You can also
                               just input only one of them but the order is important""")
    parser.add_argument('--sample_mapping_file','-smf',required=False,default=None)
    parser.add_argument('--minor_allel_frequency','-maf',required=False,default=0.05)
    parser.add_argument('--hardy_weinberg_equilibrium','-hwe',required=False,default=0.0001)
    parser.add_argument('--call_rate','-cr',required=False,default=0.95)
    parser.add_argument('--block_size','-bs',required=False,default=1500)
    parser.add_argument('--number_of_permutations','-np',required=False,default=10)
    parser.add_argument('--variant_filter','-vf',required=False,default=None)
    parser.add_argument('--feature_variant_covariate','-fvc',required=False,default=None)
    parser.add_argument('--feature_variant_filter','-fvf',required=False,default=None)
    parser.add_argument('--feature_filter','-ff',required=False,default=None)
    parser.add_argument('--seed','-s',required=False)
    parser.add_argument('--extended_annotation_file','-eaf',
                        help=
                        'Secondary annotation file, to add a multiple locations to one feature.'
                        'This can be used to either link multiple test regions to one feature or exclude multiple regions while testing a feature.', required=False)
    parser.add_argument('--relatedness_score','-rs',required=False,default=None)
    parser.add_argument('--write_permutations','-wp',action="store_true",required=False,default=False)
    parser.add_argument('--write_zscores','-wz',action="store_true",required=False,default=False)
    parser.add_argument('--write_feature_top_permutations','-wftp',action='store_true',required=False,default=False)
    parser.add_argument('--minimum_test_samples','-mts',
                    help="The minimal number of samples with non-NA values to consider a feature for a QTL test, if covariates are used the number of covariates is added to this value.",required=False,default=10)
    parser.add_argument("--gaussianize_method","-gm",
                    help="Force normal distribution on phenotypes.", default=None)
    parser.add_argument("--cis","-c",
                        action="store_true",
                        help="Run cis analysis.", default=False)
    parser.add_argument("--trans","-t",
                        action="store_true",
                        help="Run trans analysis.", default=False)
    parser.add_argument("--no_chromosome_filter","-ncf",
                        action="store_true",
                        help="Don't filter on autosomes. By default only autosomes are selected, this is where the defaults are designed for."
                        "When running on X/Y/MT please be aware that these defaults might not be appropriate.", default=False)
    parser.add_argument("--regress_covariates","-rc",
                        action="store_true",
                        help="Regress-out covariates, using a LMM, before running the QTL mapping.", default=False)
    parser.add_argument("--low_rank_random_effect","-lrc",
                        action="store_true",
                        help="Use the low rank random effect optimizations.", default=False)
    parser.add_argument("--debugger","-d",
                        action="store_true",
                        help="Print the time needed for each step", default=False)
    args = parser.parse_args()
    return args

def get_interaction_args():
    parser = argparse.ArgumentParser(description='Run QTL analysis given genotype, phenotype, and annotation.')
    parser.add_argument('--bgen','-bg',required=False)
    parser.add_argument('--plink','-pg',required=False)
    parser.add_argument('--annotation_file','-af', required=True)
    parser.add_argument('--phenotype_file','-pf', required=True)
    parser.add_argument('--output_directory','-od', required=True)
    parser.add_argument('--interaction_term','-it',
                        help=
                        'Term to use for the interaction analysis, values are extracted from the covariate matrix.'
                        'The interaction term is also taken along in the covariate matrix.',required=True,default=None)
    parser.add_argument('--window','-w', required=False,
                        help=
                        'The size of the cis window to take SNPs from.'
                        'The window will extend between:                     '
                        '    (feature_start - (window))             '
                        ' and:                                               '
                        '    (feature_end + (window))               ',default=250000)
    parser.add_argument('--genomic_range','-gr',required=False,
                        help=
                        'A genomic range to do selecte features to be considered in the analysis.'
                        'Available options: all (default), a chromsome or chromosome:start-end.',default='all')
    parser.add_argument('--covariates_file','-cf',required=False,default=None)
    parser.add_argument('--randomeff_files','-rf',required=False,default=None,
                        help = """The random effect files divided by a comma. The first one should be the kinship file and the second one should contain the read depth. You can also
                               just input only one of them but the order is important""")
    parser.add_argument('--sample_mapping_file','-smf',required=False,default=None)
    parser.add_argument('--minor_allel_frequency','-maf',required=False,default=0.05)
    parser.add_argument('--hardy_weinberg_equilibrium','-hwe',required=False,default=0.0001)
    parser.add_argument('--call_rate','-cr',required=False,default=0.95)
    parser.add_argument('--block_size','-bs',required=False,default=1500)
    parser.add_argument('--number_of_permutations','-np',required=False,default=10)
    parser.add_argument('--variant_filter','-vf',required=False,default=None)
    parser.add_argument('--feature_variant_covariate','-fvc',required=False,default=None)
    parser.add_argument('--feature_variant_filter','-fvf',required=False,default=None)
    parser.add_argument('--feature_filter','-ff',required=False,default=None)
    parser.add_argument('--seed','-s',required=False)
    parser.add_argument('--extended_annotation_file','-eaf',
                        help=
                        'Secondary annotation file, to add a multiple locations to one feature.'
                        'This can be used to either link multiple test regions to one feature or exclude multiple regions while testing a feature.', required=False)
    parser.add_argument('--regress_snp_interaction','-rsi', action="store_true",required=False,default=False)
    parser.add_argument("--regress_covariates","-rc",
                        action="store_true",
                        help="Regress-out covariates, using a LMM, before running the QTL mapping.", default=False)
    parser.add_argument('--relatedness_score','-rs',required=False,default=None)
    parser.add_argument('--write_permutations','-wp',action="store_true",required=False,default=False)
    parser.add_argument('--write_zscores','-wz',action="store_true",required=False,default=False)
    parser.add_argument('--write_feature_top_permutations','-wftp',action='store_true',required=False,default=False)
    parser.add_argument('--minimum_test_samples','-mts',
                    help="The minimal number of samples with non-NA values to consider a feature for a QTL test, if covariates are used the number of covariates is added to this value.",required=False,default=10)
    parser.add_argument("--gaussianize_method","-gm",
                    help="Force normal distribution on phenotypes.", default=None)
    parser.add_argument("--cis","-c",
                        action="store_true",
                        help="Run cis analysis.", default=False)
    parser.add_argument("--trans","-t",
                        action="store_true",
                        help="Run trans analysis.", default=False)
    parser.add_argument("--no_chromosome_filter","-ncf",
                        action="store_true",
                        help="Don't filter on autosomes. By default only autosomes are selected, this is where the defaults are designed for."
                        "When running on X/Y/MT please be aware that these defaults might not be appropriate.", default=False)
    parser.add_argument("--low_rank_random_effect","-lrc",
                        action="store_true",
                        help="Use the low rank random effect optimizations.", default=False)
    parser.add_argument("--debugger","-d",
                        action="store_true",
                        help="Print the time needed for each step", default=False)
    args = parser.parse_args()

    return args

def get_struct_args():
    parser = argparse.ArgumentParser(description='Run structLMM QTL analysis given genotype, phenotype, environments and annotation.')
    parser.add_argument('--bgen','-bg', required=False)
    parser.add_argument('--plink','-pg', required=False)
    parser.add_argument('--annotation_file','-af', required=True)
    parser.add_argument('--phenotype_file','-pf', required=True)
    parser.add_argument('--output_directory','-od', required=True)
    parser.add_argument('--environment_file','-ef', required=True)
    parser.add_argument('--window','-w', required=False,
                        help=
                        'The size of the cis window to take SNPs from.'
                        'The window will extend between:                     '
                        '    (feature_start - (window))             '
                        ' and:                                               '
                        '    (feature_end + (window))               ', default=250000)
    parser.add_argument('--genomic_range','-gr', required=False,
                        help=
                        'A genomic range to do selecte features to be considered in the analysis.'
                        'Available options: all (default), a chromsome or chromosome:start-end.', default='all')
    parser.add_argument('--covariates_file','-cf', required=False, default=None)
    parser.add_argument('--kinship_file','-kf', required=False, default=None)
    parser.add_argument('--sample_mapping_file','-smf', required=False, default=None)
    parser.add_argument('--minor_allel_frequency','-maf', required=False, default=0.05)
    parser.add_argument('--hardy_weinberg_equilibrium','-hwe', required=False, default=0.0001)
    parser.add_argument('--call_rate','-cr', required=False, default=0.95)
    parser.add_argument('--block_size','-bs', required=False, default=1500)
    parser.add_argument('--number_of_permutations','-np', required=False, default=10)
    parser.add_argument('--variant_filter','-vf', required=False, default=None)
    parser.add_argument('--feature_variant_covariate','-fvc', required=False, default=None)
    parser.add_argument('--feature_variant_filter','-fvf', required=False, default=None)
    parser.add_argument('--feature_filter','-ff', required=False, default=None)
    parser.add_argument('--seed','-s', required=False)
    parser.add_argument('--extended_annotation_file','-eaf',
                        help=
                        'Secondary annotation file, to add a multiple locations to one feature.'
                        'This can be used to either link multiple test regions to one feature or exclude multiple regions while testing a feature.', required=False)
    parser.add_argument('--relatedness_score','-rs', required=False, default=None)
    parser.add_argument('--write_permutations','-wp', action="store_true", required=False, default=False)
    parser.add_argument('--minimum_test_samples','-mts',
                        help="The minimal number of samples with non-NA values to consider a feature for a QTL test, if covariates are used the number of covariates is added to this value.", required=False, default=10)
    parser.add_argument("--gaussianize_method","-gm",
                        help="Force normal distribution on phenotypes.", default=None)
    parser.add_argument("--cis","-c",
                        action="store_true",
                        help="Run cis analysis.", default=False)
    parser.add_argument("--trans","-t",
                        action="store_true",
                        help="Run trans analysis.", default=False)
    parser.add_argument("--interaction_mode","-im",
                        action="store_true",
                        help="Run structLMM in interaction mode.", default=False)
    parser.add_argument("--association_mode","-am",
                        action="store_true",
                        help="Run structLMM in association mode.", default=False)
    parser.add_argument("--no_chromosome_filter","-ncf",
                        action="store_true",
                        help="Don't filter on autosomes. By default only autosomes are selected, this is where the defaults are designed for."
                        "When running on X/Y/MT please be aware that these defaults might not be appropriate.", default=False)
    args = parser.parse_args()
    return args

def get_grsQtl_args():
    parser = argparse.ArgumentParser(description='Run GRS QT analysis given genotype, phenotype, and annotation.')
    parser.add_argument('--genetic_risk_scores','-grs',required=True)
    parser.add_argument('--annotation_file','-af', required=True)
    parser.add_argument('--phenotype_file','-pf', required=True)
    parser.add_argument('--output_directory','-od', required=True)
    parser.add_argument('--genomic_range','-gr',required=False,
                        help=
                        'A genomic range to do selecte features to be considered in the analysis.'
                        'Available options: all (default), a chromsome or chromosome:start-end.',default='all')
    parser.add_argument('--covariates_file','-cf',required=False,default=None)
    parser.add_argument('--randomeff_files','-rf',required=False,default=None)
    parser.add_argument('--sample_mapping_file','-smf',required=False,default=None)
    parser.add_argument('--call_rate','-cr',required=False,default=0.95)
    parser.add_argument('--block_size','-bs',required=False,default=1500)
    parser.add_argument('--number_of_permutations','-np',required=False,default=10)
    parser.add_argument('--variant_filter','-vf',required=False,default=None)
    parser.add_argument('--feature_variant_covariate','-fvc',required=False,default=None)
    parser.add_argument('--feature_variant_filter','-fvf',required=False,default=None)
    parser.add_argument('--feature_filter','-ff',required=False,default=None)
    parser.add_argument('--seed','-s',required=False)
    parser.add_argument('--relatedness_score','-rs',required=False,default=None)
    parser.add_argument('--write_permutations','-wp',action="store_true",required=False,default=False)
    parser.add_argument('--write_zscores','-wz',action="store_true",required=False,default=False)
    parser.add_argument('--write_feature_top_permutations','-wftp',action='store_true',required=False,default=False)
    parser.add_argument('--minimum_test_samples','-mts',
                    help="The minimal number of samples with non-NA values to consider a feature for a QTL test, if covariates are used the number of covariates is added to this value.",required=False,default=10)
    parser.add_argument("--gaussianize_method","-gm",
                    help="Force normal distribution on phenotypes.", default=None)
    parser.add_argument("--no_chromosome_filter","-ncf",
                        action="store_true",
                        help="Don't filter on autosomes. By default only autosomes are selected, this is where the defaults are designed for."
                        "When running on X/Y/MT please be aware that these defaults might not be appropriate.", default=False)
    parser.add_argument("--regress_covariates","-rc",
                        action="store_true",
                        help="Regress-out covariates, using a LMM, before running the QTL mapping.", default=False)
    parser.add_argument("--low_rank_random_effect","-lrc",
                        action="store_true",
                        help="Use the low rank random effect optimizations.", default=False)
    parser.add_argument("--debugger","-d",
                        action="store_true",
                        help="Print the time needed for each step", default=False)
    args = parser.parse_args()
    return args
