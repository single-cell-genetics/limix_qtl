import sys
import numpy as np
import pandas
import argparse

#sys.path.append('../../hipsci_pipeline/post-processing_QTL/')

from scripts.postprocess_functions_genome_results import *
from scripts.postprocess_functions_plots import *


def get_args():
    parser = argparse.ArgumentParser(description='Run QTL analysis given genotype, phenotype, and annotation.')
    parser.add_argument('-path_data','--path_data',required=True)
    parser.add_argument('-folder_destination','--folder_destination',required=False,default=None)
    parser.add_argument('-traits','--traits',required=True)
    parser.add_argument('-trait_labels','--trait_labels',required=False,default=None)
    parser.add_argument('-chromosomes','--chromosomes',required=False,default=','.join(np.arange(1,23).astype('U')))
    parser.add_argument('-run_type',required=False,default='summary_plots_power_replication_manhattan')
    parser.add_argument('-plot_name',required=False,default='summary_plots_power_replication_manhattan')
    parser.add_argument('-p_value_field',required=False,default='p_value')
    parser.add_argument('-local_adjustment_method',required=False,default=None)
    args = parser.parse_args()
    
    return args


if __name__=='__main__':
    print('reading arguments')
    args = get_args()
    print(args)

    path_data =args.path_data
    folder_destination = args.folder_destination if args.folder_destination is not None else path_data
    traits =args.traits.split(',')
    
    run_type= args.run_type  if args.run_type is not None else ''
    chromosomes = args.chromosomes.split(',')
    trait_labels =args.trait_labels.split(',') if args.trait_labels is not None else traits
    p_value_field =args.p_value_field
    local_adjustment_method=args.local_adjustment_method
    print ('in main')
    print (p_value_field)
#path_data='/Users/mirauta/Results/hipsci/QTL1/'
#traits=['param_protein_scaled_covar_gaussnorm_test','mrna']
#run_type='plots_power_replication'
#chromosomes=['21']
#folder_destination =  path_data
#trait_labels =  traits
#print(run_type)

if 'summary' in run_type:
    print('create summary')
    for trait in traits:
        summary_gene_feature(output_file='qtl_results_genome', feature_report='ensembl_gene_id', chr_list=chromosomes,\
                             p_value_field=p_value_field,p_value_raw_field='p_value', path_data=path_data,trait=trait,local_adjustment_method=local_adjustment_method)

elif 'feature_snp'in run_type:
    print('create summary')
    for trait in traits:
        summary_gene_feature_snp(output_file='qtl_results_genome_feature_snp', feature_report='ensembl_gene_id',chr_list=chromosomes,path_data=path_data,trait=trait,thr=0.5)
    
elif 'plots' not in run_type:
    print ('run only summary')
    sys.exit()

from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec
from matplotlib_venn import venn2, venn2_circles
from matplotlib_venn import venn3, venn3_circles
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import matplotlib as mpl

if 'power' in run_type:
    a=plot_summary(folder_name=path_data, folder_destination=folder_destination+'Images_pipeline/',plot_name='qtl_summary',\
                 traits=traits, trait_labels=trait_labels,
                 file_name_results_genome='ensembl_gene_id_qtl_results_genome',   qtl_results_file='qtl_results_', \
                 colors=np.array(['orange','darkblue','green','red','k']),cis=2.5*10**5,figsize=(8,8),gene_list=None,\
                 plot_calibration_flag=False)

if 'replication_beta' in run_type:
    print("replication_beta")
    rez_pro_pep=plot_replication_beta(rez=None,path_data =path_data,\
                                 traits=traits,trait_labels= trait_labels, qtl_results_file='qtl_results_',    snp_metadata_file='snp_metadata_',    feature_metadata_file='feature_metadata_',\
                                 results_genome_file='qtl_results_genome',    feature_report='ensembl_gene_id',\
                                 folder_destination=folder_destination+'Images_pipeline/', figsize=6,red_dots_features=None, \
                                 p_value_field=p_value_field,thr=0.1)

if 'replication_pv' in run_type:
    rez_pro_pep=plot_replication_pv(rez=None,path_data =path_data,\
                traits=traits,trait_labels= trait_labels, qtl_results_file='qtl_results_',\
                snp_metadata_file='snp_metadata_',    feature_metadata_file='feature_metadata_',\
                results_genome_file='qtl_results_genome',    feature_report='ensembl_gene_id',\
                folder_destination=folder_destination+'Images_pipeline/', figsize=6,red_dots_features=None, \
                p_value_field=p_value_field,thr=0.01)
    
    names=['empirical_feature_p_value', 'replicated_p_value', 'replicated_self_p_value','beta', 'replicated_beta', 'feature_id', 'gene_name', 'snp_id', 'chromosome', 'strand', 'position', 'ensembl_gene_id']
    df= pandas.DataFrame(data=np.array([rez_pro_pep[key] for key in names ]).T, index=rez_pro_pep['ensembl_gene_id'],columns=names)
    df.to_csv(path_or_buf=path_data+traits[0]+'_'+traits[1]+'_qtl_results.txt',mode='w', sep='\t', columns=None, header=True, index=True)
    pandas.DataFrame(df['snp_id'][df['empirical_feature_p_value'].values.astype(float)<0.01]).to_csv(path_or_buf=path_data+traits[0]+'_qtl_results_significant_snps.txt',mode='w', sep='\t', columns=None, header=False, index=False)
     
#    genes1=df.index[(df['p_value'].values.astype(float)<10**-4)&(df['replicated_p_value'].values.astype(float)<10**-4)]
#print (genes1.shape)

if 'manhattan' in run_type:
    df=pandas.read_table(path_data+traits[0]+'_'+traits[1]+'_qtl_results.txt', sep='\t',index_col=0)
#    genes1=df.index[(df['p_value'].values.astype(float)<10**-4)&(df['replicated_p_value'].values.astype(float)<10**-4)]
    genes1=df.index[(df[p_value_field].values.astype(float)<10**-7)]
    for i, g in enumerate(genes1):
        plot_manhatan_alone( gene_ensembl_id= genes1[i],folder_name=path_data,\
                 folder_destination=folder_destination+'Images_pipeline/'+'/Manhattan/',\
                 plot_name='manhattan_'+traits[0]+traits[1]+'_',traits=traits,trait_labels=trait_labels,file_name_results_genome='ensembl_gene_id_qtl_results_genome',\
                 qtl_results_file='qtl_results_',colors=np.array(['black','green','orange','blue']),p_value_field=p_value_field, figsize=4)

#        plt.show()



    