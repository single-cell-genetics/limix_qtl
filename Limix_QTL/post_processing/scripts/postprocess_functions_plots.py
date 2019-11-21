from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec
import sys

from matplotlib_venn import venn2, venn2_circles
from matplotlib_venn import venn3, venn3_circles
import matplotlib.pyplot as plt
import seaborn as sb
from matplotlib import colors as mcolors
import matplotlib as mpl
import h5py
import numpy as np
#
#import limix.stats.fdr as FDR
import os
import scipy.stats as scst
import statsmodels.sandbox.stats.multicomp as scst2
sys.path.append('../')
sys.path.append('/Users/mirauta/Git/hipsci_pipeline/post-processing_QTL/')
from scripts.postprocess_functions_genome_results import *
import copy

#
#
#
#folder_name='/Users/mirauta/Results/hipsci/QTL1/'
#file_name_qtl='qtl_results_'
#file_name_perm='perm_results_'
#traits=['param_protein_scaled_peer_gaussnorm_test'];
#chromosome='22'
def plot_summary_onechr_perm(plot_name='qtl_summary_onechr',folder_name=None,folder_destination=None,\
                 traits=[''],chromosome='21',
                 file_name_qtl='qtl_results_',file_name_perm='perm_results_', \
                 colors=np.array(['orange','darkblue','green','m']),cis=2.5*10**5, figsize=(8,5)):
    
    if folder_destination is None:
        folder_destination =folder_name
    if not os.path.exists(folder_destination):
        os.makedirs(folder_destination)
 
    featureh5=[h5py.File(folder_name+'/'+trait+'/'+file_name_qtl+chromosome+'.h5','r') for trait in traits]
    featurepermh5=[h5py.File(folder_name+'/'+trait+'/'+file_name_perm+chromosome+'.h5','r') for trait in traits]
    gene_list_common=np.array([np.array(list(fh5.keys()))  for indf,fh5 in enumerate(featureh5)])
#        return gene_list_common
    temp=np.unique(np.hstack(gene_list_common),return_counts=1)
    gene_list_common=temp[0][temp[1]==gene_list_common.shape[0]]
#    local_adjusted_common=np.array([np.array([fh5[gene]['summary_data/min_p_value_local_adjusted'][:][0] for gene in np.intersect1d(gene_list_common,np.array(list(fh5.keys())))])for indf,fh5 in enumerate(featureh5)])
    fh5=featurepermh5[0]
    perm_fields=np.array(fh5[list(fh5.keys())[0]].dtype.names)[['permutation' in field for field in np.array(fh5[list(fh5.keys())[0]].dtype.names)]]
    
    
    rez={}
    
    rez['p_value']=np.array([np.hstack([fh5[gene]['p_value']for gene in gene_list_common])  for indf,fh5 in enumerate(featureh5)])
    for field in perm_fields:
        rez[field]=np.array([np.hstack([fh5[gene][field] for gene in gene_list_common])  for indf,fh5 in enumerate(featurepermh5)])
    rez['p_value_min_bonferroni']=np.array([np.hstack([np.nanmin(fh5[gene]['p_value'])*fh5[gene]['p_value'].shape[0] for gene in gene_list_common])  for indf,fh5 in enumerate(featureh5)])[0]
    for field in perm_fields:
        rez[field+'_min_bonferroni']=np.array([np.hstack([np.nanmin(fh5[gene][field])*fh5[gene][field].shape[0] for gene in gene_list_common])  for indf,fh5 in enumerate(featurepermh5)])[0]
#    else:
#        local_adjusted=np.array([np.array([fh5[gene]['summary_data/min_p_value_local_adjusted'][:][0] for gene in np.intersect1d(gene_list,np.array(list(fh5.keys())))])for indf,fh5 in enumerate(featureh5)])

 
    fig=plt.figure(figsize=figsize)
    mpl.rcParams.update(mpl.rcParamsDefault)
    fig.patch.set_facecolor('white')
    axes = fig.add_subplot(1, 2, 1, axisbg='white')
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.yaxis.set_ticks_position('left')
    axes.xaxis.set_ticks_position('bottom')
    fdr=np.array([5,4,3,2,1])[::-1]
    trait_labels=['qval']+list(perm_fields)
    local_adjusted=[rez['p_value_min_bonferroni']]+[rez[field+'_min_bonferroni'] for field in perm_fields]
    for iy,yy in enumerate(local_adjusted):
#        print(yy)
 
#        plt.plot(fdr[::-1] ,[(-np.log10(FDR.qvalues1(yy))>thr).sum() for thr in fdr],color=colors[iy],markersize=3,label=trait_labels[iy])
        plt.plot(fdr[::-1] ,[(-np.log10(scst2.fdrcorrection0(yy)[1])>thr).sum() for thr in fdr],color=colors[iy],label=trait_labels[iy],lw=3)
    axes.plot((fdr[-1],fdr[-1]),(0,np.max([(a<10**-2).sum() for a in local_adjusted])) ,'--',color='k',lw=0.5)
    axes.plot((fdr[0],fdr[-1]),(np.max([(a<10**-2).sum() for a in local_adjusted]),np.max([(a<10**-2).sum() for a in local_adjusted])) ,'--',color='k',lw=0.5)
    plt.xlabel('-log10 FDR',fontsize=13);
    plt.ylabel('#pGenes',fontsize=13)
    plt.xticks(fdr,fdr[::-1])
    if trait_labels is not None: plt.legend(loc=2,fontsize=9)
    
#==============================================================================
#     power plot commmon
#==============================================================================
    plt.subplot(1,2,2)
    
    plt.title('Calibration',fontsize=10)
    
    for iy,yy in enumerate([rez['p_value']]+[rez[field] for field in perm_fields]):
        yy=np.sort(yy[(yy==yy)&(yy!=1)])
        plt.plot(-np.log10((0.5+np.arange(len(yy)))/len(yy))[::-1],-np.log10(yy)[::-1],'o',color=colors[iy],markersize=4,label=trait_labels[iy])
                
    plt.plot(-np.log10((0.5+np.arange(len(yy)))/len(yy))[::-1],-np.log10((0.5+np.arange(len(yy)))/len(yy))[::-1],'k',lw=1)
    plt.xlabel('-log10 PV',fontsize=13);
    plt.ylabel('-log10 random',fontsize=13)
    plt.legend(loc=2,fontsize=8)
    for f in featureh5: f.close()
    plt.tight_layout()
    plt.savefig(folder_destination+plot_name+'_'+'_'.join(trait_labels)+'.pdf',dpi=600)
    return [gene_list_common,rez]


def plot_summary_vs2_2(plot_name='qtl_summary',folder_name=None,folder_destination=None,\
                 traits=[''],trait_labels=None,
                 file_name_results_genome='Feature_results_genome', qtl_results_file='QTL_results_',\
                 colors=np.array(['orange','darkblue','green','m','k']),cis=2.5*10**5,gene_list=None, \
                 fig=None,axes=None):
    
    if folder_destination is None:
        folder_destination =folder_name
    if not os.path.exists(folder_destination):
        os.makedirs(folder_destination)

 
    featureh5=[h5py.File(folder_name+'/'+trait+'_'+file_name_results_genome+'.h5','r') for trait in traits]
 
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.yaxis.set_ticks_position('left')
    axes.xaxis.set_ticks_position('bottom')
    def get_distance_tss(gene,fh5):
        if ( fh5[gene]['metadata/feature_strand'][:].astype('U')[0]=='+'):
            return fh5[gene]['summary_data/min_p_value_position'][:][0]-fh5[gene]['metadata/start'][:][0]
        else:
            return -fh5[gene]['summary_data/min_p_value_position'][:][0]+fh5[gene]['metadata/end'][:][0]

    if gene_list is None:
        distance_tss=  [np.array([get_distance_tss(gene,fh5) for gene in  np.array(list(fh5.keys()))[ local_adjusted[indf]<10**-3]])\
                        for indf,fh5 in enumerate(featureh5)]
    else:
        distance_tss=  [np.array([get_distance_tss(gene,fh5) for gene in \
                                  np.intersect1d(gene_list[indf],np.array(list(fh5.keys())))])\
                        for indf,fh5 in enumerate(featureh5)]
        
    axes.hist([np.clip(d[np.isfinite(d)],-cis,cis) for d in distance_tss], bins=7,width=20000,\
               label=trait_labels,color=colors[:len(traits)],normed=1)
    plt.ylabel('density',fontsize=13);
    plt.xlabel('distance from TSS',fontsize=13);
    plt.legend(loc=2,fontsize=9)
    plt.xticks(np.linspace(-cis,cis,5),np.array([str(int(l))+'k' for l in np.linspace(-cis,cis,5)/1e3]))
#    plt.yticks('off')
                    
    for f in featureh5: f.close()
    
    return 1

def plot_summary_vs2_1(plot_name='qtl_summary',folder_name=None,folder_destination=None,traits=[''],trait_labels=None,
                 file_name_results_genome='Feature_results_genome', qtl_results_file='QTL_results_',\
                 colors=np.array(['orange','darkblue','green','m','k']),cis=2.5*10**5, gene_list=None,\
                 fig=None,axes=None):
    
    if folder_destination is None:
        folder_destination =folder_name
    if not os.path.exists(folder_destination):
        os.makedirs(folder_destination)

    featureh5=[h5py.File(folder_name+'/'+trait+'_'+file_name_results_genome+'.h5','r') for trait in traits]

    featureh5=[h5py.File(folder_name+'/'+trait+'_'+file_name_results_genome+'.h5','r') for trait in traits]

    if gene_list is None:
        local_adjusted=np.array([np.array([fh5[gene]['summary_data/min_p_value_local_adjusted'][:][0] for gene in fh5.keys()]) \
                                 for indf,fh5 in enumerate(featureh5)])
 
    else:
        local_adjusted=np.array([np.array([fh5[gene]['summary_data/min_p_value_local_adjusted'][:][0] for gene in np.intersect1d(gene_list,np.array(list(fh5.keys())))])for indf,fh5 in enumerate(featureh5)])

    fdr=np.array([5,4,3,2,1])[::-1]
 
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.yaxis.set_ticks_position('left')
    axes.xaxis.set_ticks_position('bottom')
    for iy,yy in enumerate(local_adjusted):

        axes.plot(fdr[::-1] ,[(-np.log10(scst2.fdrcorrection0(yy)[1])>thr).sum() for thr in fdr],color=colors[iy],\
                  label=trait_labels[iy],lw=3)
    axes.plot((fdr[-1],fdr[-1]),(0,np.max([(a<10**-2).sum() for a in local_adjusted])) ,'--',color='k',lw=0.5)
    axes.plot((fdr[0],fdr[-1]),(np.max([(a<10**-2).sum() for a in local_adjusted]),np.max([(a<10**-2).sum() for a in local_adjusted])) ,'--',color='k',lw=0.5)
    plt.xlabel('-log10 FDR',fontsize=13);
    plt.ylabel('#pGenes',fontsize=13)
    plt.xticks(fdr,fdr[::-1])
    if trait_labels is not None: plt.legend(loc=2,fontsize=9)

    for f in featureh5: f.close()
    
    return local_adjusted


def plot_summary(plot_name='qtl_summary',folder_name=None,folder_destination=None, traits=[''],trait_labels=None,
                 file_name_results_genome='Feature_results_genome', qtl_results_file='QTL_results_',\
                 colors=np.array(['orange','darkblue','green','m','k']),cis=2.5*10**5, figsize=(12,12),\
                                gene_list=None,\
                 plot_tss_distance_flag=False,plot_calibration_flag=False):
    
    if folder_destination is None:
        folder_destination =folder_name
    if not os.path.exists(folder_destination):
        os.makedirs(folder_destination)
    
    if trait_labels is None: trait_labels=traits
#    
    fig=plt.figure(figsize=figsize)
    mpl.rcParams.update(mpl.rcParamsDefault)
    fig.patch.set_facecolor('white')
 
    featureh5=[h5py.File(folder_name+'/'+trait+'_'+file_name_results_genome+'.h5','r') for trait in traits]
 
    if gene_list is None:
        local_adjusted=np.array([np.array([fh5[gene]['summary_data/min_p_value_local_adjusted'][:][0] for gene in fh5.keys()]) \
                                 for indf,fh5 in enumerate(featureh5)])
        gene_list_common=np.array([np.array(list(fh5.keys()))  for indf,fh5 in enumerate(featureh5)])
        temp=np.unique(np.hstack(gene_list_common),return_counts=1)
        gene_list_common=temp[0][temp[1]==gene_list_common.shape[0]]
        local_adjusted_common=np.array([np.array([fh5[gene]['summary_data/min_p_value_local_adjusted'][:][0] for gene in np.intersect1d(gene_list_common,np.array(list(fh5.keys())))])for indf,fh5 in enumerate(featureh5)])
        
 
    else:
        local_adjusted=np.array([np.array([fh5[gene]['summary_data/min_p_value_local_adjusted'][:][0] for gene in np.intersect1d(gene_list,np.array(list(fh5.keys())))])for indf,fh5 in enumerate(featureh5)])
        gene_list_common=gene_list
        temp=np.unique(np.hstack(gene_list_common),return_counts=1)
        gene_list_common=temp[0][temp[1]==gene_list_common.shape[0]]
        local_adjusted_common=np.array([np.array([fh5[gene]['summary_data/min_p_value_local_adjusted'][:][0] for gene in \
                                                  np.intersect1d(gene_list_common,np.array(list(fh5.keys())))])for indf,fh5 in enumerate(featureh5)])
        
 
    axes = fig.add_subplot(2, 2, 1, axisbg='white')
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.yaxis.set_ticks_position('left')
    axes.xaxis.set_ticks_position('bottom')
    fdr=np.array([5,4,3,2,1])[::-1]
 
    for iy,yy in enumerate(local_adjusted):
#        print(yy)
 
#        plt.plot(fdr[::-1] ,[(-np.log10(FDR.qvalues1(yy))>thr).sum() for thr in fdr],color=colors[iy],markersize=3,label=trait_labels[iy])
        plt.plot(fdr[::-1] ,[(-np.log10(scst2.fdrcorrection0(yy)[1])>thr).sum() for thr in fdr],color=colors[iy],label=trait_labels[iy],lw=3)
    axes.plot((fdr[-1],fdr[-1]),(0,np.max([(a<10**-2).sum() for a in local_adjusted])) ,'--',color='k',lw=0.5)
    axes.plot((fdr[0],fdr[-1]),(np.max([(a<10**-2).sum() for a in local_adjusted]),np.max([(a<10**-2).sum() for a in local_adjusted])) ,'--',color='k',lw=0.5)
    plt.xlabel('-log10 FDR',fontsize=13);
    plt.ylabel('#pGenes',fontsize=13)
    plt.xticks(fdr,fdr[::-1])
    if trait_labels is not None: plt.legend(loc=2,fontsize=9)
    
#==============================================================================
#     power plot commmon
#==============================================================================
    axes = fig.add_subplot(2, 2, 2, axisbg='white')
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.yaxis.set_ticks_position('left')
    axes.xaxis.set_ticks_position('bottom')
    fdr=np.array([5,4,3,2,1])[::-1]
 
    for iy,yy in enumerate(local_adjusted_common):
#        print(yy)
 
#        plt.plot(fdr[::-1] ,[(-np.log10(FDR.qvalues1(yy))>thr).sum() for thr in fdr],color=colors[iy],markersize=3,label=trait_labels[iy])
        plt.plot(fdr[::-1] ,[(-np.log10(scst2.fdrcorrection0(yy)[1])>thr).sum() for thr in fdr],color=colors[iy],label=trait_labels[iy],lw=3)
    axes.plot((fdr[-1],fdr[-1]),(0,np.max([(a<10**-2).sum() for a in local_adjusted_common])) ,'--',color='k',lw=0.5)
    axes.plot((fdr[0],fdr[-1]),(np.max([(a<10**-2).sum() for a in local_adjusted_common]),np.max([(a<10**-2).sum() for a in local_adjusted_common])) ,'--',color='k',lw=0.5)
    plt.xlabel('-log10 FDR',fontsize=13);
    plt.ylabel('#pGenes',fontsize=13)
    plt.xticks(fdr,fdr[::-1])
    plt.title('Common genes')
    if trait_labels is not None: plt.legend(loc=2,fontsize=9)


    if plot_tss_distance_flag:
        axes = fig.add_subplot(2, 2, 3, axisbg='white')
        axes.spines['top'].set_visible(False)
        axes.spines['right'].set_visible(False)
        axes.yaxis.set_ticks_position('left')
        axes.xaxis.set_ticks_position('bottom')
        def get_distance_tss(gene,fh5):
            if ( fh5[gene]['metadata/feature_strand'][:].astype('U')[0]=='+'):
                return fh5[gene]['summary_data/min_p_value_position'][:][0]-fh5[gene]['metadata/start'][:][0]
            else:
                return -fh5[gene]['summary_data/min_p_value_position'][:][0]+fh5[gene]['metadata/end'][:][0]
    
        if gene_list is None:
            distance_tss=  [np.array([get_distance_tss(gene,fh5) for gene in  np.array(list(fh5.keys()))[ local_adjusted[indf]<10**-3]]) for indf,fh5 in enumerate(featureh5)]
        else:
            distance_tss=  [np.array([get_distance_tss(gene,fh5) for gene in np.intersect1d(gene_list,np.array(list(fh5.keys())))[ local_adjusted[indf]<10**-3]]) for indf,fh5 in enumerate(featureh5)]
        axes.hist([np.clip(d[np.isfinite(d)],-cis,cis) for d in distance_tss], bins=7,width=15000,label=trait_labels,color=colors[:len(traits)])
        plt.ylabel('#pGenes \n(PV<0.001)',fontsize=13);plt.xlabel('distance from TSS',fontsize=13);
        plt.legend(loc=2,fontsize=9)
        plt.xticks(np.linspace(-cis,cis,5),np.array([str(int(l))+'k' for l in np.linspace(-cis,cis,5)/1e3]))
            

    if plot_calibration_flag:
        plt.subplot(2,2,4)
    
        featureh5=[h5py.File(folder_name+'/'+trait+'_'+file_name_results_genome+'.h5','r') for trait in traits]
    
        y=[np.hstack([np.hstack([fh5[feature_id]['data/p_value_raw'][f][:] for f in fh5[feature_id]['data/p_value_raw'].keys()]).flatten()  for feature_id in list(fh5.keys())]) for indf,fh5 in enumerate(featureh5)]

        plt.title('Calibration',fontsize=10)
        
        for iy,yy in enumerate(y):
            yy=np.sort(yy[(yy==yy)&(yy!=1)])
            try:
                plt.plot(-np.log10((0.5+np.arange(len(yy)))/len(yy))[::-1],-np.log10(yy)[::-1],'o',color=colors[iy],markersize=4,label=feature_label[iy])
            except:
                plt.plot(-np.log10((0.5+np.arange(len(yy)))/len(yy))[::-1],-np.log10(yy)[::-1],'o',color=colors[iy],markersize=4)
                    
        plt.plot(-np.log10((0.5+np.arange(len(yy)))/len(yy))[::-1],-np.log10((0.5+np.arange(len(yy)))/len(yy))[::-1],'k',lw=1)
        plt.legend(loc=2,fontsize=8)

    for f in featureh5: f.close()
    plt.tight_layout()
    plt.savefig(folder_destination+plot_name+'_'+'_'.join(trait_labels)+'.png',dpi=600)
    return local_adjusted


def plot_manhatan_alone(folder_name='/Users/mirauta/Data/MS/hipsci/TMT/',folder_destination='/Users/mirauta/Data/MS/hipsci/TMT/Images',plot_name='manhattan',\
                        traits=None,trait_labels=None,file_name_results_genome='ensembl_gene_id_qtl_results_genome',   qtl_results_file='qtl_results_',colors=np.array(['k','b','g','m']), figsize=4, gene_ensembl_id= 'ENSG00000182154',\
                        p_value_field='p_value',log_flag=True,ylim=None,savefig=True,fplot=None,ax=None,pdf=True,\
                        ann_snp=None):
    if folder_destination is None:
        folder_destination =folder_name+'/manhattan/'
    if not os.path.exists(folder_destination):
        os.makedirs(folder_destination)


#    print(traits)
    featureh5=[h5py.File(folder_name+'/'+trait+'_'+file_name_results_genome+'.h5','r')[gene_ensembl_id] for trait in traits]
    nfeatures=len(featureh5)
#    print ('featureh5')   
    rez={}
    temppos=[np.hstack([fh5['data']['position'][f][:] for f in fh5['data']['position'].keys()]) for fh5 in featureh5]
    temp=np.unique(np.hstack(temppos).flatten(),return_counts=1); 
    commonpos=temp[0][temp[1]==np.max(temp[1])]
    
    for dat in [p_value_field,'position']:
        rez[dat]=[np.array([fh5['data'][dat][f][:][np.in1d(fh5['data']['position'][f][:],commonpos)] for f in fh5['data'][dat].keys()]) for fh5 in featureh5]
 
    #==============================================================================
    #  modifty positions for plot
    #==============================================================================
    xpos=(commonpos-commonpos.min())/10.0**6;
    axxpos=np.array([i for i in np.linspace(min(xpos),max(xpos),5)])
    showxpos=np.array([int(i) for i in np.linspace(min(commonpos),max(commonpos),5)])
    startstop=np.array([featureh5[0]['metadata']['start'][:][0],featureh5[0]['metadata']['end'][:][0]]);
    startstop=(startstop-commonpos.min())/10.0**6

    if fplot is None:
        fig=plt.figure(figsize=(figsize*3,figsize*nfeatures))
        fig.set_facecolor('white')
        fplot = gridspec.GridSpec(nfeatures*3,8)
#        fplot.update(hspace=200, wspace=10)    
        ax = [plt.subplot(fplot[(i*3):(i*3+3),:10]) for i in np.arange(nfeatures)]
    
    gene_name=str(featureh5[0]['metadata']['gene_name'][:].astype('U')[0])
    
    feature_name=str(featureh5[0]['metadata']['feature_id'][:].astype('U')[0])
    colors2=colors[1:]
    firsttrait=rez[p_value_field][0][np.argsort([tt.min() for tt in rez[p_value_field][0]])[0]]
    
    for indf, a in enumerate(ax):
        if ylim is not None: a.set_ylim(ylim)
        a.spines['top'].set_visible(False);    a.spines['right'].set_visible(False);a.yaxis.set_ticks_position('left'); a.xaxis.set_ticks_position('bottom');a.set_axis_bgcolor('white');
        a.add_patch(Rectangle((startstop[0], 0), startstop[1]-startstop[0], 1, facecolor="grey",    alpha=0.35))
        a.set_xticks( axxpos)
        a.set_xticklabels([str(np.around(i,1))+' Mb' for i in showxpos/10.0**6])
        for indt,t in  enumerate(rez[p_value_field][indf][np.argsort([tt.min() for tt in rez[p_value_field][indf]])[:14]]):
            if log_flag: a.plot(xpos,-np.log10(t),'o',color=colors[indt],markersize=2.25)
            else: a.plot(xpos, t,'o',color=colors[indt],markersize=2.25)
        a.set_ylabel(trait_labels[indf]+'\n'+gene_name+"\n -log10PV",labelpad=30,fontsize=10,rotation=0,horizontalalignment= 'center' ,verticalalignment= 'center')
    
        for indt,t in enumerate(rez[p_value_field][indf][np.argsort([tt.min() for tt in rez[p_value_field][indf]])[:14]]): 
            a.plot(xpos[np.argsort(t)[:5]],-np.log10(t[np.argsort(t)[:5]]),'*',color=colors[indt],markersize=5)
            
#            for indf2, a2 in enumerate(ax):
#                if indf!=indf2:
        
        
        for indt,t in enumerate(rez[p_value_field][indf][np.argsort([tt.min() for tt in rez[p_value_field][indf]])[:14]]): 
            a.plot(xpos[np.argsort(firsttrait)[:10]],-np.log10(t[np.argsort(firsttrait)[:10]],),'ro',markersize=3)
#        plt.ylim(0,5)
 
#    print (str(featureh5[0]['metadata']['gene_name'][:].astype('U')[0]))
        if ann_snp is not None:
            if (ann_snp[indf].loc['gwas_snp']==ann_snp[indf].loc['gwas_snp']):
                x=(np.array([s.split('_')[1] for s in ann_snp[indf].loc['gwas_snp'].split(';')]).astype(int)[0]-commonpos.min())/10.0**6
                a.plot(x,-np.log10(np.nanmin(rez[p_value_field][indf])),'ro')
                a.annotate(ann_snp[indf].loc['gwas_trait'].split(';')[0], xy=(x, -np.log10(np.nanmin(rez[p_value_field][indf]))),fontsize=11)

    plt.tight_layout()
    if savefig:
        if pdf:
            plt.savefig(folder_destination+plot_name+str(featureh5[0]['metadata']['gene_name'][:].astype('U')[0])+'_manhattan.pdf',bbox_inches='tight')
        plt.savefig(folder_destination+plot_name+str(featureh5[0]['metadata']['gene_name'][:].astype('U')[0])+'_manhattan.png',dpi=600,bbox_inches='tight')



def plot_manhatan_genes(folder_name='/Users/mirauta/Data/MS/hipsci/TMT/',\
                     folder_destination='/Users/mirauta/Data/MS/hipsci/TMT/Images',plot_name='manhattan',\
                      trait=None,trait_labels=None,file_name_results_genome='ensembl_gene_id_qtl_results_genome',\
                      qtl_results_file='qtl_results_',colors=np.array(['b','k','g','m']), figsize=4, genes_ensembl_id= [],\
                      p_value_field='p_value_raw',log_flag=True,ylim=None,savefig=False,fplot=None,ax=None):
    
    if folder_destination is None:
        folder_destination =folder_name+'/manhattan/'
    if not os.path.exists(folder_destination):
        os.makedirs(folder_destination)

    featureh5=[h5py.File(folder_name+'/'+trait+'_'+file_name_results_genome+'.h5','r')[gene] for gene in genes_ensembl_id]

    genes_names=[fh5['metadata/gene_name'][:][0].astype('U') for fh5 in featureh5]  

    rez={}
    for fh5 in featureh5:
        rez[fh5]={}
        for subf in fh5['data']['position'].keys():
            rez[fh5][subf] =pd.DataFrame(data=np.array([fh5['data'][dat][subf][:].astype('U') for dat in [p_value_field,'position','snp_id']]).T,\
               columns=[p_value_field,'position','snp_id'])
        
   
    
    for fh5 in featureh5:
        for subf in fh5['data']['position'].keys():
#            rez[fh5][subf]['position']=rez[fh5][subf]['position'].astype(int)
            rez[fh5][subf]=rez[fh5][subf].set_index('snp_id',drop=0)
            
    
    
    snps=np.unique(np.hstack([np.hstack([rez[fh5][subf].index   for subf in fh5['data']['position'].keys()]) for fh5 in featureh5]))
    for fh5 in featureh5:
        for subf in fh5['data']['position'].keys():            
            rez[fh5][subf]=rez[fh5][subf].loc[snps]
            rez[fh5][subf]['position']=np.array([s.split('_')[1] for s in  rez[fh5][subf].index]).astype(int)
            rez[fh5][subf][p_value_field][rez[fh5][subf][p_value_field]!=rez[fh5][subf][p_value_field]]=1
            rez[fh5][subf]['chromosome']=np.array([s.split('_')[0] for s in rez[fh5][subf].index]).astype(int)
            rez[fh5][subf]['position2']=copy.deepcopy(rez[fh5][subf]['position'])
            lengths={};curr=0;lengths[1]=0
            for chr in np.unique(rez[fh5][subf]['chromosome']): 
                
                lengths[chr+1]= curr+np.max(rez[fh5][subf]['position'] [rez[fh5][subf]['chromosome'] ==chr])
                curr=int(lengths[chr+1])
  

            for ichr, chr in enumerate(np.unique(rez[fh5][subf]['chromosome'])):
                rez[fh5][subf]['position2'][rez[fh5][subf]['chromosome']==chr]=rez[fh5][subf]['position'][rez[fh5][subf]['chromosome']==chr].values.astype(int)+lengths[chr]
        rez[fh5]['gene']= fh5['metadata']['start'][:][0]+lengths[fh5['metadata']['chromosome'][:][0]]

    if fplot is None:
        fig=plt.figure(figsize=(figsize*2,figsize*len(featureh5)))
        fig.set_facecolor('white')
        fplot = gridspec.GridSpec(len(featureh5)*6,8)
        fplot.update(hspace=10, wspace=10)    
        ax = [plt.subplot(fplot[(i*3):(i*3+3),:10]) for i in np.arange(len(featureh5))]
    
    
    colors2=colors[1:]
    
    for indf, a in enumerate(ax):
        fh5=featureh5[indf]
        subf=list(fh5['data']['position'].keys())[0] 
        if ylim is not None: a.set_ylim(ylim)
        a.spines['top'].set_visible(False);    a.spines['right'].set_visible(False);a.yaxis.set_ticks_position('left'); a.xaxis.set_ticks_position('bottom');a.set_axis_bgcolor('white');
        a.add_patch(Rectangle((rez[fh5]['gene']-10**7, 0), 2*10**7, 4, facecolor="grey",    alpha=0.99935))
        a.set_xticks(np.hstack([(np.array(list(lengths.values()))[:-1]+np.array(list(lengths.values()))[1:])/2,lengths[22]]))
#        order=np.argsort(rez['position2'][indf][0])
        a.set_xticklabels(np.array(list(lengths.keys())))
        for subf in fh5['data']['position'].keys():    
            y=-np.log10(rez[fh5][subf][p_value_field].values.astype(float))
            for thr in [0,1,2,3,4,6]:
                a.plot(rez[fh5][subf]['position2'][y>=thr],y[y>=thr],'o',color=colors[indf],markersize=min(thr/1.5+1,4))
#            print(rez[fh5][subf]['snp_id'][np.argmax(y)])
#            print(np.nanmax(y))
#            a.annotate(rez[fh5][subf]['snp_id'][np.argmax(y)], xy=(rez[fh5][subf]['position2'][np.argmax(y)], np.nanmax(y)),fontsize=9)
            
        a.set_ylabel(genes_names[indf]+"\n -log10PV", labelpad=40,fontsize=10,rotation=0,horizontalalignment= 'center' ,verticalalignment= 'center')
       
        

    if savefig:
        plt.savefig(folder_destination+plot_name+str(featureh5[0]['metadata']['gene_name'][:].astype('U')[0])+'.pdf',dpi=600)
#colors=np.array(['b','k','g','m'])
#path_data='/Users/mirauta/Results/hipsci/QTL_new/';folder_name=path_data
#ensembl_genes_id=np.array(['ENSG00000010292','ENSG00000109805', 'ENSG00000113810', 'ENSG00000136824'])
#trait='param_protein_scaled_peer_log_trans001'
#genes_ensembl_id=np.intersect1d(ensembl_genes_id ,\
#                                list(h5py.File(path_data+'/'+trait+'_ensembl_gene_id_qtl_results_genome.h5','r').keys()))
#trans_ensembl_id=genes_ensembl_id
#
#fig=plt.figure(figsize=(10,10))
#fig.set_facecolor('white')
#fplot = gridspec.GridSpec(3*6,10)
#fplot.update(hspace=10, wspace=10)
#cis_ensembl_gene_id=trans_ensembl_id[0]
#ax = [plt.subplot(fplot[:5,:3]) for i in np.arange(len([trait,trait]))]
#plot_manhatan_alone(folder_name=path_data,folder_destination=path_data+'Images_pipeline/'+'/Manhattan/',plot_name='complexes_manhattan',\
#                    traits=[trait,trait],trait_labels=[trait,trait], gene_ensembl_id= cis_ensembl_gene_id,\
#                    p_value_field='p_value_raw', savefig=0,fplot=fplot,ax=ax)
#ax = [plt.subplot(fplot[(i*3):(i*3+3),3:10]) for i in np.arange(len(trans_ensembl_id))]
#plot_manhatan_genes( genes_ensembl_id= trans_ensembl_id,folder_name=path_data,\
#                    folder_destination=path_data+'Images_pipeline/'+'/Manhattan/',\
#                    plot_name='complex_manhattan_'+trait,trait=trait,\
#                    p_value_field='p_value_raw', figsize=4,fplot=fplot,ax=ax)
#plt.savefig(path_data+'Images_pipeline/'+'/Manhattan/summary_test.png')
#
#plt.show()
#==============================================================================
#==============================================================================
# # 
#==============================================================================
#==============================================================================
 

def plot_replication_beta(rez=None,path_data ='/Users/mirauta/Data/MS/hipsci/TMT/', path_data2 =None,   traits=['protein_test','peptide_test'],trait_labels=['protein_test','peptide_test'],\
    qtl_results_file='qtl_results_',    snp_metadata_file='snp_metadata_',    feature_metadata_file='feature_metadata_',\
    results_genome_file='qtl_results_genome',    feature_report='ensembl_gene_id',folder_destination='/Users/mirauta/Results/hipsci/Images_pipeline',\
    figsize=5, red_dots_features=None,red_dot='ro',plot_name='',p_value_field='p_value',thr=0.01):
    if rez is None:
        rez=replication_two_features(path_data =path_data, path_data2 =path_data2,    traits=np.array(traits), qtl_results_file=qtl_results_file,    snp_metadata_file=snp_metadata_file,   \
                                     feature_metadata_file=feature_metadata_file, results_genome_file=results_genome_file,    feature_report= feature_report,p_value_field=p_value_field,thr=thr)
   
    fig=plt.figure(figsize=(figsize,figsize))
    mpl.rcParams.update(mpl.rcParamsDefault)
    axes = fig.add_subplot(1, 1, 1, axisbg='white')
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.spines['top'].set_visible(False)
    axes.spines['left'].set_visible(True)
    axes.yaxis.set_ticks_position('left')
    axes.xaxis.set_ticks_position('bottom')
    fig.patch.set_facecolor('white')
    plt.plot((rez['beta']),(rez['replicated_beta']),'o',color='grey',markersize=6,label='Spearman: '+\
             str(np.around(scst.spearmanr(rez['beta'][np.isfinite(rez['replicated_beta']+rez['beta'])],rez['replicated_beta'][np.isfinite(rez['replicated_beta']+rez['beta'])])[0],2)))
    plt.plot((rez['beta']),(rez['replicated_beta']),'o',color='grey',markersize=6,label='Pearson: '+\
             str(np.around(np.corrcoef(rez['beta'][np.isfinite(rez['replicated_beta']+rez['beta'])],rez['replicated_beta'][np.isfinite(rez['replicated_beta']+rez['beta'])])[0,1],2)))
#    thrs=0.001;
#    plt.plot((rez['beta'][rez['replicated_self_beta']<thrs]),(rez['replicated_beta'][rez['replicated_self_beta']<thrs]),'.',markersize=8,color='darkorange')
    
    if red_dots_features is not None:        
        plt.plot((rez['beta'][np.in1d(rez['feature_id'],red_dots_features)]), (rez['replicated_beta'][np.in1d(rez['feature_id'],red_dots_features)]),red_dot,markersize=6)#, mfc='none')
        
    plt.plot((np.nanmin(rez['beta']),(np.nanmax(rez['beta']))),(np.nanmin(rez['beta']),(np.nanmax(rez['beta']))),'k--',lw=0.25)

    plt.xlabel(trait_labels[0]+'\n beta ');plt.ylabel(trait_labels[1]+'\n beta',rotation=90 )
    plt.legend(loc=1)
    plt.savefig(folder_destination+'replication_beta_'+traits[0]+'_'+traits[1]+plot_name+'.png')
    return rez

def plot_replication_pv(rez=None,path_data ='/Users/mirauta/Data/MS/hipsci/TMT/', path_data2 =None,   traits=['protein_test','peptide_test'],trait_labels=['protein_test','peptide_test'],\
    qtl_results_file='qtl_results_',    snp_metadata_file='snp_metadata_',    feature_metadata_file='feature_metadata_',\
    results_genome_file='qtl_results_genome',    feature_report='ensembl_gene_id',folder_destination='/Users/mirauta/Results/hipsci/Images_pipeline',\
    figsize=5, red_dots_features=None,red_dot='ro',plot_name='',p_value_field='p_value',thr=0.001):
    if rez is None:
        rez=replication_two_features(path_data =path_data, path_data2 =path_data2,    traits=np.array(traits), qtl_results_file=qtl_results_file, snp_metadata_file=snp_metadata_file,\
                                     feature_metadata_file=feature_metadata_file, results_genome_file=results_genome_file,  feature_report= feature_report,p_value_field=p_value_field,thr=thr)
   
    fig=plt.figure(figsize=(figsize,figsize))
    mpl.rcParams.update(mpl.rcParamsDefault)
    axes = fig.add_subplot(1, 1, 1, axisbg='white')
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.spines['top'].set_visible(False)
    axes.spines['left'].set_visible(True)
    axes.yaxis.set_ticks_position('left')
    axes.xaxis.set_ticks_position('bottom')
    fig.patch.set_facecolor('white')
    plt.plot(-np.log10(rez[p_value_field]),-np.log10(rez['replicated_p_value']),'o',color='grey',markersize=6,label='Spearman: '+\
             str(np.around(scst.spearmanr(rez[p_value_field][np.isfinite(rez['replicated_p_value']+rez[p_value_field])],\
                                          rez['replicated_p_value'][np.isfinite(rez['replicated_p_value']+rez[p_value_field])])[0],2)))
    plt.plot(-np.log10(rez[p_value_field]),-np.log10(rez['replicated_p_value']),'o',color='grey',markersize=6,label='Pearson: '+\
             str(np.around(np.corrcoef(np.log(rez[p_value_field][np.isfinite(rez['replicated_p_value']+rez[p_value_field])]),np.log(rez['replicated_p_value'][np.isfinite(rez['replicated_p_value']+rez[p_value_field])]))[0,1],2)))

    plt.plot(-np.log10(rez[p_value_field]),-np.log10(rez['replicated_p_value']),'o',color='grey',markersize=6)
    thrs=0.001;
    plt.plot(-np.log10(rez[p_value_field][rez['replicated_self_p_value']<thrs]),-np.log10(rez['replicated_p_value'][rez['replicated_self_p_value']<thrs]),'.',markersize=8,color='darkorange')
    plt.ylim(0,np.nanmax(-np.log10(rez['replicated_p_value'])))
    if red_dots_features is not None:        
        plt.plot(-np.log10(rez[p_value_field][np.in1d(rez['feature_id'],red_dots_features)]), -np.log10(rez['replicated_p_value'][np.in1d(rez['feature_id'],red_dots_features)]),red_dot,markersize=6)#, mfc='none')
        
    plt.plot((2,-np.log10(np.nanmin(rez[p_value_field]))),(2,-np.log10(np.nanmin(rez[p_value_field]))),'k--',lw=0.25)

    plt.xlabel(trait_labels[0]+'\n - log10 PV');plt.ylabel(trait_labels[1]+'\n - log10 PV',rotation=90 )
    plt.legend(loc=1)
    plt.savefig(folder_destination+'replication_pv_'+traits[0]+'_'+traits[1]+plot_name+'.png')
    return rez



def render_mpl_table(data, col_width=3.0, row_height=0.625, font_size=14,
                     header_color='#40466e', row_colors=['#f1f1f2', 'w'], edge_color='w',
                     bbox=[0, 0, 1, 1], header_columns=0,
                     ax=None, **kwargs):
    import six
    if ax is None:
        size = (np.array(data.shape[::-1]) + np.array([0, 1])) * np.array([col_width, row_height])
        fig, ax = plt.subplots(figsize=size)
        ax.axis('off')

    mpl_table = ax.table(cellText=data.values, bbox=bbox, colLabels=data.columns, **kwargs)

    mpl_table.auto_set_font_size(False)
    mpl_table.set_fontsize(font_size)

    for k, cell in  six.iteritems(mpl_table._cells):
        cell.set_edgecolor(edge_color)
        if k[0] == 0 or k[1] < header_columns:
            cell.set_text_props(weight='bold', color='w')
            cell.set_facecolor(header_color)
        else:
            cell.set_facecolor(row_colors[k[0]%len(row_colors) ])
    return ax

    
def plot_intersect3_venn(folder_dest,qv,thrs,comment='',labels=['eQTL','pQTL_1','pQTL_2']):
    
    fig = plt.figure(figsize=(5, 5))
    fig.patch.set_facecolor('white')
    mpl.rcParams.update(mpl.rcParamsDefault)
    axes = fig.add_subplot(1, 1, 1, axisbg='white')
#    ids=np.intersect1d(e1['ensembl_ID'],np.intersect1d(e2['ensembl_ID'],e3['ensembl_ID']))
#    qv=np.vstack([np.array([e['eqv_min'][np.in1d(e['ensembl_ID'],ids)]for e in [e1]]),np.array([e['qv_min'][np.in1d(e['ensembl_ID'],ids)]for e in [e2,e3]])])
#    fracs=[np.sum((ee>np.min(ee[eq>thrs]))&(pp<np.min(pp[pq>thrs]))),np.sum((pp>np.min(pp[pq>thrs]))&(ee<np.min(ee[eq>thrs]))),np.sum((pp>np.min(pp[pq>thrs]))&(ee>np.min(ee[eq>thrs])))]
    
    fracs=[np.sum((qv[0]<thrs)&(qv[1]>thrs)&(qv[2]>thrs)),np.sum((qv[0]>thrs)&(qv[1]<thrs)&(qv[2]>thrs)),np.sum((qv[0]<thrs)&(qv[1]<thrs)&(qv[2]>thrs)),\
    np.sum((qv[0]>thrs)&(qv[1]>thrs)&(qv[2]<thrs)),np.sum((qv[0]<thrs)&(qv[1]>thrs)&(qv[2]<thrs)),np.sum((qv[0]>thrs)&(qv[1]<thrs)&(qv[2]<thrs)),\
    np.sum((qv[0]<thrs)&(qv[1]<thrs)&(qv[2]<thrs))]
    
    v=venn3(subsets=fracs, set_labels =labels)
    v.get_patch_by_id('100').set_alpha(1.0)
    v.get_patch_by_id('100').set_color('lightskyblue')
    #v.get_label_by_id('A').set_text('Set "A"')
    c = venn3_circles(subsets=fracs, lw=0.5)
    for text in v.set_labels:    text.set_fontsize(10)
#    for text in out.subset_labels:    text.set_fontsize(16)
  
    #plt.annotate('Unknown set', xy=v.get_label_by_id('100').get_position() - np.array([0, 0.05]), xytext=(-70,-70),
#             ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
#             arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'))

    plt.savefig(folder_dest+"venn3_eql_pqtl"+comment+".png",dpi=800)
      




def plot_manhatan_custom(path_data='/Users/mirauta/Results/hipsci/QTL_jan/',path_destination='/Users/mirauta/Results/hipsci/manuscript_images/Manhattan',\
                         plot_name='manhattan',\
                        traits=None,trait_labels=None,file_name_results_genome='ensembl_gene_id_qtl_results_genome',\
                        colors=np.array(['darkgreen','grey','darkblue','deepskyblue']), figsize=4, \
                        gene_ensembl_id= 'ENSG00000182154',\
                        p_value_field='p_value',log_flag=True,ylim=None,savefig=True,fplot=None,ax=None,pdf=True,transcripts=None,\
                        ann_snp=None):

    if not os.path.exists(path_destination):
        os.makedirs(path_destination)


    featureh5=[h5py.File(path_data+'/'+trait+'_'+file_name_results_genome+'.h5','r')[gene_ensembl_id] for trait in traits]
    
    rez={}
    temppos=[np.hstack([fh5['data']['position'][f][:] for f in fh5['data']['position'].keys()]) for fh5 in featureh5]
    temp=np.unique(np.hstack(temppos).flatten(),return_counts=1); 
    commonpos=temp[0][temp[1]==np.max(temp[1])]
    
    for dat in [p_value_field,'position']:
        rez[dat]=[np.array([fh5['data'][dat][f][:][np.in1d(fh5['data']['position'][f][:],commonpos)] for f in fh5['data'][dat].keys()]) for fh5 in featureh5]
    if transcripts is not None:
        fh5=featureh5[3]
        for dat in [p_value_field,'position']:
            rez[dat][3]=np.array([fh5['data'][dat][f][:][np.in1d(fh5['data']['position'][f][:],commonpos)] for f in transcripts])
#    return(rez)
    bestrez0=np.argmin(np.min(rez[p_value_field][0],1))
    #==============================================================================
    #  modifty positions for plot
    #==============================================================================
    xpos=(commonpos-commonpos.min())/10.0**6;
    axxpos=np.array([i for i in np.linspace(min(xpos),max(xpos),5)])
    showxpos=np.array([int(i) for i in np.linspace(min(commonpos),max(commonpos),5)])
    startstop=np.array([featureh5[0]['metadata']['start'][:][0],featureh5[0]['metadata']['end'][:][0]]);
    startstop=(startstop-commonpos.min())/10.0**6

    nfeatures=2
    fig=plt.figure(figsize=(figsize*3,figsize*nfeatures))
    fig.set_facecolor('white')
    fplot = gridspec.GridSpec(nfeatures*3,8)
    ax = [plt.subplot(fplot[(i*3):(i*3+3),:10]) for i in np.arange(nfeatures)]
    
    gene_name=str(featureh5[0]['metadata']['gene_name'][:].astype('U')[0])  
    feature_name=str(featureh5[0]['metadata']['feature_id'][:].astype('U')[0])
    colors2=colors[1:]
    firsttrait=rez[p_value_field][0][np.argsort([tt.min() for tt in rez[p_value_field][0]])[0]]
    
    for indf, a in enumerate(ax):
        if ylim is not None: a.set_ylim(ylim)
        a.spines['top'].set_visible(False);    a.spines['right'].set_visible(False);a.yaxis.set_ticks_position('left'); a.xaxis.set_ticks_position('bottom');a.set_axis_bgcolor('white');
        a.add_patch(Rectangle((startstop[0], 0), startstop[1]-startstop[0], 1, facecolor="grey",    alpha=0.35))
        a.set_xticks( axxpos)
        a.set_xticklabels([str(np.around(i,1))+' Mb' for i in showxpos/10.0**6])
    
    t=rez[p_value_field][0][bestrez0];
    ax[0].plot(xpos,-np.log10(t),'o',color='steelblue',markersize=4, mfc='none')
    ax[0].set_ylabel(gene_name+'\nProtein \n -log10PV',labelpad=30,fontsize=10,rotation=0,horizontalalignment= 'center' ,verticalalignment= 'center')
    
 
    

    markers=['P','*']*100
    for indt,t in  enumerate(rez[p_value_field][3]):
         ax[1].plot(xpos,-np.log10(t),markers[indt],color=colors[indt],markersize=4, mfc='none')
        
    if ann_snp is not None:
       if (ann_snp[indf].loc['gwas_snp']==ann_snp[indf].loc['gwas_snp']):
                x=(np.array([s.split('_')[1] for s in ann_snp[indf].loc['gwas_snp'].split(';')]).astype(int)[0]-commonpos.min())/10.0**6
                a.plot(x,-np.log10(np.nanmin(rez[p_value_field][indf])),'ro')
                a.annotate(ann_snp[indf].loc['gwas_trait'].split(';')[0], xy=(x, -np.log10(np.nanmin(rez[p_value_field][indf]))),fontsize=11)
    t=rez[p_value_field][2][0];
    ax[1].plot(xpos,-np.log10(t),'o',color='grey',markersize=4, mfc='none')
    ax[1].set_ylabel('mRNA\n(gene & transcript)\n -log10PV',labelpad=30,fontsize=10,rotation=0,horizontalalignment= 'center' ,verticalalignment= 'center')
    
    plt.tight_layout()
    if savefig:
#        if pdf:
#            plt.savefig(path_destination+plot_name+str(featureh5[0]['metadata']['gene_name'][:].astype('U')[0])+'_manhattan.pdf',bbox_inches='tight')
        plt.savefig(path_destination+plot_name+str(featureh5[0]['metadata']['gene_name'][:].astype('U')[0])+'_manhattan.png',dpi=600,bbox_inches='tight')


def plot_manhatan_custom2(path_data='/Users/mirauta/Results/hipsci/QTL_jan/',path_destination='/Users/mirauta/Results/hipsci/manuscript_images/Manhattan',\
                         plot_name='manhattan',\
                        traits=None,trait_labels=None,file_name_results_genome='ensembl_gene_id_qtl_results_genome',\
                        colors=np.array(['darkgreen','grey','darkblue','deepskyblue']), figsize=4, \
                        gene_ensembl_id= 'ENSG00000182154',\
                        p_value_field='p_value',log_flag=True,ylim=None,savefig=True,fplot=None,ax=None,pdf=True,\
                        ann_snp=None,transcripts=None):

    if not os.path.exists(path_destination):
        os.makedirs(path_destination)


    featureh5=[h5py.File(path_data+'/'+trait+'_'+file_name_results_genome+'.h5','r')[gene_ensembl_id] for trait in traits]
    
    rez={}
    temppos=[np.hstack([fh5['data']['position'][f][:] for f in fh5['data']['position'].keys()]) for fh5 in featureh5]
    temp=np.unique(np.hstack(temppos).flatten(),return_counts=1); 
    commonpos=temp[0][temp[1]==np.max(temp[1])]
    
    for dat in [p_value_field,'position']:
        rez[dat]=[np.array([fh5['data'][dat][f][:][np.in1d(fh5['data']['position'][f][:],commonpos)] for f in fh5['data'][dat].keys()]) for fh5 in featureh5]
    if transcripts is not None:
        fh5=featureh5[3]
        for dat in [p_value_field,'position']:
            rez[dat][3]=np.array([fh5['data'][dat][f][:][np.in1d(fh5['data']['position'][f][:],commonpos)] for f in transcripts])
            #    return(rez)
    bestrez0=np.argmin(np.min(rez[p_value_field][0],1))
    #==============================================================================
    #  modifty positions for plot
    #==============================================================================
    xpos=(commonpos-commonpos.min())/10.0**6;
    axxpos=np.array([i for i in np.linspace(min(xpos),max(xpos),5)])
    showxpos=np.array([int(i) for i in np.linspace(min(commonpos),max(commonpos),5)])
    startstop=np.array([featureh5[0]['metadata']['start'][:][0],featureh5[0]['metadata']['end'][:][0]]);
    startstop=(startstop-commonpos.min())/10.0**6

    nfeatures=3
    fig=plt.figure(figsize=(figsize*3,figsize*nfeatures))
    fig.set_facecolor('white')
    fplot = gridspec.GridSpec(nfeatures*3,8)
    ax = [plt.subplot(fplot[(i*2):(i*2+2),:10]) for i in np.arange(nfeatures)]
    
    gene_name=str(featureh5[0]['metadata']['gene_name'][:].astype('U')[0])  
    feature_name=str(featureh5[0]['metadata']['feature_id'][:].astype('U')[0])
    colors2=colors[1:]
    firsttrait=rez[p_value_field][0][np.argsort([tt.min() for tt in rez[p_value_field][0]])[0]]
    
    for indf, a in enumerate(ax):
        if ylim is not None: a.set_ylim(ylim)
        a.spines['top'].set_visible(False);    a.spines['right'].set_visible(False);a.yaxis.set_ticks_position('left'); a.xaxis.set_ticks_position('bottom');a.set_axis_bgcolor('white');
        a.add_patch(Rectangle((startstop[0], 0), startstop[1]-startstop[0], 1, facecolor="grey",    alpha=0.35))
        a.set_xticks( axxpos)
        a.set_xticklabels([str(np.around(i,1))+' Mb' for i in showxpos/10.0**6])
    
    t=rez[p_value_field][0][bestrez0];
    ax[0].plot(xpos,-np.log10(t),'o',color='steelblue',markersize=4, mfc='none')
    ax[0].set_ylabel(gene_name+'\nProtein \n -log10PV',labelpad=30,fontsize=10,rotation=0,horizontalalignment= 'center' ,verticalalignment= 'center')
    
 

    ax[1].set_ylabel('mRNA\n(gene)\n -log10PV',labelpad=30,fontsize=10,rotation=0,horizontalalignment= 'center' ,verticalalignment= 'center')
    ax[2].set_ylabel('mRNA isoform\n -log10PV',labelpad=30,fontsize=10,rotation=0,horizontalalignment= 'center' ,verticalalignment= 'center')
    
    t=rez[p_value_field][2][0];
    ax[1].plot(xpos,-np.log10(t),'o',color='grey',markersize=4, mfc='none')
    ax[1].set_ylim(0,np.nanmax(-np.log10(rez[p_value_field][3])))   
    markers=['P','*']*100
    for indt,t in  enumerate(rez[p_value_field][3]):
         ax[2].plot(xpos,-np.log10(t),markers[indt],color=colors[indt],markersize=4, mfc='none')
        
    if ann_snp is not None:
       if (ann_snp[indf].loc['gwas_snp']==ann_snp[indf].loc['gwas_snp']):
                x=(np.array([s.split('_')[1] for s in ann_snp[indf].loc['gwas_snp'].split(';')]).astype(int)[0]-commonpos.min())/10.0**6
                a.plot(x,-np.log10(np.nanmin(rez[p_value_field][indf])),'ro')
                a.annotate(ann_snp[indf].loc['gwas_trait'].split(';')[0], xy=(x, -np.log10(np.nanmin(rez[p_value_field][indf]))),fontsize=11)

    plt.tight_layout()
    if savefig:
#        if pdf:
#            plt.savefig(path_destination+plot_name+str(featureh5[0]['metadata']['gene_name'][:].astype('U')[0])+'_manhattan.pdf',bbox_inches='tight')
        plt.savefig(path_destination+plot_name+str(featureh5[0]['metadata']['gene_name'][:].astype('U')[0])+'_manhattan.png',dpi=600,bbox_inches='tight')

           