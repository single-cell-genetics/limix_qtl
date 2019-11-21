import sys
import h5py
import numpy as np
import pandas
import pandas as pd
import glob
import statsmodels.sandbox.stats.multicomp as scst2

def local_adjustment(pv, N=1,  method=''):

    if method is None:
        return np.hstack(pv)*pv.shape[0]
    if method=='Bonferroni': 
        N=np.hstack(pv).shape[0]
        return pv*N
    else:
        print('Valid multiple testing correction  methods are: None;Bonferroni')


        
def summary_gene_feature_divided(qtl_results_file='qtl_results_',snp_metadata_file='snp_metadata_', \
                                 feature_metadata_file='feature_metadata_',output_file='qtl_results_genome',\
                                 feature_report='ensembl_gene_id',feature_report_inh5='ensembl_gene_id',chr_list=[9],path_data=None,path_dataout=None,folder_qtl=None, \
                                 p_value_field='p_value',p_value_raw_field='p_value',local_adjustment_method='Bonferroni', exclude_snps=['']):
    if path_dataout is None: path_dataout=path_data
    _doc=" aggregates qtl results to feature_report level"
    iichr=0
    for ichr,chr in enumerate(chr_list):
        
        print ('chromosome: '+str(chr))
        files=glob.glob(path_data+folder_qtl+'/'+feature_metadata_file+str(chr)+'_*')
        files=np.array([f.replace(path_data+folder_qtl+'/'+feature_metadata_file,'').replace('.txt','') for f in files])
        print (files)
        for file in files:
            print (file)
            try:
                frez=h5py.File(path_data+'/'+folder_qtl+'/'+qtl_results_file+file+'.h5','r')
                frezkeys= np.array([k.replace('_i_','') for k in list(frez.keys())])
                ffea= pandas.read_table(path_data+'/'+folder_qtl+'/'+feature_metadata_file+ file+'.txt', sep='\t')
                fsnp= pandas.read_table(path_data+'/'+folder_qtl+'/'+snp_metadata_file+ file+'.txt', sep='\t').set_index('snp_id',drop=False).transpose()
                print (ffea.columns.values)
            except:
                print('chromosome'+file+' missing')
                continue
            
            if ~np.in1d(feature_report,ffea.columns.values):
                ffea[feature_report]=ffea[feature_report_inh5]
            
        
            indexnan=np.where((ffea[feature_report])!=(ffea[feature_report]))[0]
            for i in indexnan:
                ffea[feature_report][i]='gene_'+str(ffea['chromosome'][i])+'_'+str(ffea['start'][i])
            ffea_feature=ffea.set_index('feature_id', drop=False).transpose()
            ffea_report_feature=ffea.set_index(feature_report, drop=False).transpose()
            
                                 
            if iichr==0:
                fOut=h5py.File(path_dataout+'/'+folder_qtl+'_'+feature_report+'_'+output_file+'.h5','w')
                iichr=1
            else:
                fOut=h5py.File(path_dataout+'/'+folder_qtl+'_'+feature_report+'_'+output_file+'.h5','r+')
        
            # for each report_feature  create h5 groups
            count=0
            for ifea,report_feature in enumerate(np.unique(list(ffea_report_feature))):
#                print (report_feature)
           
                #select features for which qtl was computed
                features=np.intersect1d(np.array( ffea_report_feature[report_feature].transpose()['feature_id']), frezkeys)
                if len(features) >=1:
                    
                    try:
                        fg=fOut.create_group(report_feature)
                    
                        pv= np.array([frez[f][p_value_field] for f in  features ])
                        for i, ppv in  enumerate(pv): ppv[ppv!=ppv]=1
                        pv2= np.array([frez[f][p_value_raw_field] for f in  features ])
                        for i, ppv in  enumerate(pv2): ppv[ppv!=ppv]=1
                        beta= np.array([frez[f]['beta'] for f in  features ])
                        for i, b in  enumerate(beta): b [b!=b]=1
                        beta_se= np.array([frez[f]['beta_se']**2 for f in  features ])
                        for i, b in  enumerate(beta_se): b [b!=b]=1
                        snp_id= np.array([frez[f]['snp_id'] for f in  features ])
                        position= np.array([fsnp[snp_id[indf].astype('U')].transpose()['position'] for indf, f in  enumerate(features) ])
                        for i, p in  enumerate(position): p [p!=p]=1
            
                        fgm=fg.create_group('metadata')
                        for key in ffea_feature[features[0]].keys():
                 
                            if isinstance(ffea_feature[features[0]][key],int) :
                                fgm.create_dataset(key,data=np.array([ffea_feature[f][key] for f in  features ]))
                            else:
                                fgm.create_dataset(key,data=np.array([ffea_feature[f][key] for f in  features ]).astype('S'))
                        
                        fgd=fg.create_group('data')
                        fgd.create_dataset('features',data= features.astype('S'))
                        fgdp=fgd.create_group(p_value_field) 
                        for indf,f in enumerate(features): fgdp.create_dataset(f,data= pv[indf])
        
                        fgdp2=fgd.create_group('p_value_raw')
                        for indf,f in enumerate(features): fgdp2.create_dataset(f,data= pv2[indf])
                        
                        fgdb=fgd.create_group('beta') 
                        for indf,f in enumerate(features): fgdb.create_dataset(f,data=beta[indf])
                        fgdbse=fgd.create_group('beta_se') 
                        for indf,f in enumerate(features): fgdbse.create_dataset(f,data=beta_se[indf])
                        fgdpo=fgd.create_group('position') 
                        for indf,f in enumerate(features): fgdpo.create_dataset(f,data=position[indf].astype(int))
                        fgds=fgd.create_group('snp_id') 
                        for indf,f in enumerate(features): fgds.create_dataset(f,data=snp_id[indf])
                                                                                  
                                                         
                        fgs=fg.create_group('summary_data')
                        fgs.create_dataset('min_p_value',data= np.nanmin(np.hstack(pv)) [None])
                        fgs.create_dataset('min_p_value_beta',data=np.hstack(beta)[np.nanargmin(np.hstack(pv))][None])
                        fgs.create_dataset('min_p_value_beta_se',data=np.hstack(beta_se)[np.nanargmin(np.hstack(pv))][None])
                        
                        p_bonf=np.nanmin(local_adjustment((pv),method=local_adjustment_method));
                        if p_bonf>1:p_bonf=np.array(1)
                        fgs.create_dataset('min_p_value_local_adjusted',data=p_bonf[None])
                        
                        minfeatureindex=np.argmin([np.nanmin(ppv) *len(ppv) for ppv in pv])
                        fgs.create_dataset('min_p_value_feature_id',data= np.array(fgm ['feature_id'][minfeatureindex].astype('S'))[None])
                
                        min_snp_id=frez[features[minfeatureindex]]['snp_id'][np.nanargmin(pv[minfeatureindex])].astype('U')
                        fgs.create_dataset('min_p_value_snp_id',data= np.array(min_snp_id).astype('S')[None])
                        fgs.create_dataset('min_p_value_position',data= np.array(fsnp[min_snp_id]['position'])[None])
                
                        count+=1
                    except:1
    
            frez.close()
            fOut.close()

    
def summary_gene_feature_snp_all_files(qtl_results_file='qtl_results_',snp_metadata_file='snp_metadata_', feature_metadata_file='feature_metadata_', output_file='qtl_results_genome_feature_snp',\
                                       feature_report_inh5='ensembl_gene_id',feature_report='ensembl_gene_id',chr_list=None,path_data=None,path_dataout=None,folder_qtl=None,thr=0.1,df_snp_file=None,\
                                       df_filter_field='',data_filter_field='',append=False):

    _doc=" aggregates qtl results to feature_report level"
    if chr_list is None:
        files=glob.glob(path_data+folder_qtl+'/'+feature_metadata_file+'*')
        files=np.array([f.replace(path_data+folder_qtl+'/'+feature_metadata_file,'').replace('.txt','') for f in files])
    else:
        files=chr_list
    print ('files')
    print (files)
    for ichr,chr in enumerate(files):
        if ichr%10==0:print ('chromosome: '+str(ichr)+'from'+str(len(files)))
    
        try:
            frez=h5py.File(path_data+'/'+folder_qtl+'/'+qtl_results_file+str(chr)+'.h5','r')
            frezkeys= np.array([k.replace('_i_','') for k in list(frez.keys())])
            ffea= pandas.read_table(path_data+'/'+folder_qtl+'/'+feature_metadata_file+ str(chr)+'.txt', sep='\t')
            fsnp= pandas.read_table(path_data+'/'+folder_qtl+'/'+snp_metadata_file+ str(chr)+'.txt', sep='\t').set_index('snp_id',drop=False).transpose()
            if ichr%10: print (ffea.columns.values)
        except:
            print('chromosome'+str(chr)+' missing')
            continue

        if(ffea.shape[0]==0):
            print('No features in feature annotation, for block '+chr)
            continue
        if(fsnp.shape[0]==0):
            print('No features in snps annotation, for block '+chr)
            continue
        
        if ~np.in1d(feature_report,ffea.columns.values):
            ffea[feature_report]=ffea[feature_report_inh5]
        
            
        indexnan=np.where((ffea[feature_report])!=(ffea[feature_report]))[0]
        for i in indexnan:
            ffea[feature_report][i]='gene_'+str(ffea['chromosome'][i])+'_'+str(ffea['start'][i])
        ffea_feature=ffea.set_index('feature_id', drop=False).transpose()

        ffea_report_feature=ffea.set_index(feature_report, drop=False).transpose()
        
        count=0
        data={}
        for key in ['empirical_feature_p_value','p_value','beta','beta_se','snp_id','feature_id',feature_report, 'n_samples']:
            data[key]=np.zeros(len(list(ffea_report_feature)),dtype='object')+np.nan
#        if feature_report!='ensembl_gene_id':
#            data['ensembl_gene_id']=np.zeros(len(list(ffea_report_feature)),dtype='object')+np.nan
                
        for ifea,report_feature in enumerate(np.unique(list(ffea_report_feature))):

            features=np.intersect1d(np.array( ffea_report_feature[report_feature].transpose()['feature_id']), frezkeys)
            if len(features) >=1:
                for key in ['empirical_feature_p_value','p_value','beta','beta_se','snp_id', 'n_samples']:
                    temp = np.array([frez[f][key] for f in  features ])
                    data[key][ifea]=np.hstack(temp).astype('U')
            data['feature_id'][ifea]=np.hstack([np.repeat(f,len(frez[f][key])) for f in  features ])
            if feature_report!='feature_id':
                data[feature_report][ifea]=np.hstack([np.repeat(report_feature,len(frez[f][key])) for f in  features ])

        for key in data.keys():
            data[key]=np.hstack(data[key]) 
#        if feature_report!='ensembl_gene_id':
#            if 'ENSG00' in data['feature_id'][ifea]:
#                data['ensembl_gene_id']=np.array([g.split('.')[0].split('_')[0] for g in data['feature_id']])
#            data['q_value_empirical_feature_p_value']=scst2.fdrcorrection0(
        if df_snp_file is None:
            temp=pd.DataFrame(data).iloc[data['p_value'].astype(float)<thr]
            temp=temp.set_index('snp_id')
            temp.to_csv(path_or_buf=path_dataout+folder_qtl+'_qtl_results_feature_snp_'+str(thr)+'.txt', mode='w'if ((ichr==0)&(~append)) else 'a', sep='\t', columns=None, header=ichr==0, index=True)
        else:
            df_snp=pandas.read_table(path_data+df_snp_file+'.txt', sep='\t',index_col=0)
            temp=pd.DataFrame(data).set_index(data_filter_field,drop=0)
            temp=temp.loc[np.intersect1d(temp.index,df_snp[df_filter_field][df_snp[df_filter_field]==df_snp[df_filter_field]])]
            temp=temp.set_index('snp_id')
            temp.to_csv(path_or_buf=path_dataout+folder_qtl+'_qtl_results_feature_snp_'+df_snp_file+'.txt',\
                    mode='w'if ((ichr==0)&(~append)) else 'a', sep='\t', columns=None, header=ichr==0, index=True)

    
def summary_gene_feature_snp(qtl_results_file='qtl_results_',snp_metadata_file='snp_metadata_', feature_metadata_file='feature_metadata_',output_file='qtl_results_genome_feature_snp',
                         feature_report='ensembl_gene_id',chr_list=[9],path_data=None,trait=None,thr=0.1):

    _doc=" aggregates qtl results to feature_report level"
    
    for ichr,chr in enumerate(chr_list):
        print ('chromosome: '+str(chr))
    
        try:
            frez=h5py.File(path_data+'/'+trait+'/'+qtl_results_file+str(chr)+'.h5','r')
            frezkeys= np.array([k.replace('_i_','') for k in list(frez.keys())])
            ffea= pandas.read_table(path_data+'/'+trait+'/'+feature_metadata_file+ str(chr)+'.txt', sep='\t')
            fsnp= pandas.read_table(path_data+'/'+trait+'/'+snp_metadata_file+ str(chr)+'.txt', sep='\t').set_index('snp_id',drop=False).transpose()
            print (ffea.columns.values)
        except:
            print('chromosome'+str(chr)+' missing')
            continue
        indexnan=np.where((ffea[feature_report])!=(ffea[feature_report]))[0]
        for i in indexnan:
            ffea[feature_report][i]='gene_'+str(ffea['chromosome'][i])+'_'+str(ffea['start'][i])
        ffea_feature=ffea.set_index('feature_id', drop=False).transpose()
        ffea_report_feature=ffea.set_index(feature_report, drop=False).transpose()
        
        count=0
        data={}
        for key in ['empirical_feature_p_value','p_value','beta','snp_id','feature_id',feature_report]:
            data[key]=np.zeros(len(list(ffea_report_feature)),dtype='object')+np.nan
            
        for ifea,report_feature in enumerate(np.unique(list(ffea_report_feature))):

            features=np.intersect1d(np.array( ffea_report_feature[report_feature].transpose()['feature_id']), frezkeys)
            if len(features) >=1:
                for key in ['empirical_feature_p_value','p_value','beta','snp_id']:
                    temp = np.array([frez[f][key] for f in  features ])
                    data[key][ifea]=np.hstack(temp).astype('U')
            data['feature_id'][ifea]=np.hstack([np.repeat(f,len(frez[f][key])) for f in  features ])
            data[feature_report][ifea]=np.hstack([np.repeat(report_feature,len(frez[f][key])) for f in  features ])
            
        for key in data.keys():
            data[key]=np.hstack(data[key]) 
            
        pd.DataFrame(data).iloc[data['p_value'].astype(float)<thr].to_csv(path_or_buf=path_data+trait+'_qtl_results_feature_snp_'+str(thr)+'.txt',\
                    mode='w'if ichr==0 else 'a', sep='\t', columns=None, header=ichr==0, index=True)

        
def summary_gene_feature(qtl_results_file='qtl_results_',snp_metadata_file='snp_metadata_', feature_metadata_file='feature_metadata_',output_file='qtl_results_genome',\
                         feature_report='ensembl_gene_id',chr_list=[9],path_data=None,folder_qtl=None, \
                         p_value_field='p_value',p_value_raw_field='p_value',local_adjustment_method='Bonferroni', exclude_snps=['']):

    _doc=" aggregates qtl results to feature_report level"
    iichr=0
    for ichr,chr in enumerate(chr_list):
        print ('chromosome: '+str(chr))
    
        try:
            frez=h5py.File(path_data+'/'+folder_qtl+'/'+qtl_results_file+str(chr)+'.h5','r')
            frezkeys= np.array([k.replace('_i_','') for k in list(frez.keys())])
            ffea= pandas.read_table(path_data+'/'+folder_qtl+'/'+feature_metadata_file+ str(chr)+'.txt', sep='\t')
            fsnp= pandas.read_table(path_data+'/'+folder_qtl+'/'+snp_metadata_file+ str(chr)+'.txt', sep='\t').set_index('snp_id',drop=False).transpose()
            print (ffea.columns.values)
        except:
            print('chromosome'+str(chr)+' missing')
            continue
        
        if ~np.in1d(feature_report,ffea.columns.values):
            ffea[feature_report]=ffea['feature_id']
            
        
        indexnan=np.where((ffea[feature_report])!=(ffea[feature_report]))[0]
        for i in indexnan:
            ffea[feature_report][i]='gene_'+str(ffea['chromosome'][i])+'_'+str(ffea['start'][i])
        ffea_feature=ffea.set_index('feature_id', drop=False).transpose()
        ffea_report_feature=ffea.set_index(feature_report, drop=False).transpose()
        
                             
        if iichr==0:
            fOut=h5py.File(path_data+'/'+folder_qtl+'_'+feature_report+'_'+output_file+'.h5','w')
            iichr=1
        else:
            fOut=h5py.File(path_data+'/'+folder_qtl+'_'+feature_report+'_'+output_file+'.h5','r+')
    
        # for each report_feature  create h5 groups
        count=0
        for ifea,report_feature in enumerate(np.unique(list(ffea_report_feature))):
#            print (report_feature)
       
            #select features for which qtl was computed
            features=np.intersect1d(np.array( ffea_report_feature[report_feature].transpose()['feature_id']), frezkeys)
            if len(features) >=1:
                
                fg=fOut.create_group(report_feature)
                
                pv= np.array([frez[f][p_value_field] for f in  features ])
                for i, ppv in  enumerate(pv): ppv[ppv!=ppv]=1
                pv2= np.array([frez[f][p_value_raw_field] for f in  features ])
                for i, ppv in  enumerate(pv2): ppv[ppv!=ppv]=1
                beta= np.array([frez[f]['beta'] for f in  features ])
                for i, b in  enumerate(beta): b [b!=b]=1
                beta_se= np.array([frez[f]['beta_se'] for f in  features ])
                for i, b in  enumerate(beta_se): b [b!=b]=1
                snp_id= np.array([frez[f]['snp_id'] for f in  features ])
                position= np.array([fsnp[snp_id[indf].astype('U')].transpose()['position'] for indf, f in  enumerate(features) ])
                for i, p in  enumerate(position): p [p!=p]=1
    
                fgm=fg.create_group('metadata')
                for key in ffea_feature[features[0]].keys():
         
                    if isinstance(ffea_feature[features[0]][key],int) :
                        fgm.create_dataset(key,data=np.array([ffea_feature[f][key] for f in  features ]))
                    else:
                        fgm.create_dataset(key,data=np.array([ffea_feature[f][key] for f in  features ]).astype('S'))
                
                fgd=fg.create_group('data')
                fgd.create_dataset('features',data= features.astype('S'))
                fgdp=fgd.create_group(p_value_field) 
                for indf,f in enumerate(features): fgdp.create_dataset(f,data= pv[indf])

                fgdp2=fgd.create_group('p_value_raw')
                for indf,f in enumerate(features): fgdp2.create_dataset(f,data= pv2[indf])
                
                fgdb=fgd.create_group('beta') 
                for indf,f in enumerate(features): fgdb.create_dataset(f,data=beta[indf])
                fgdbse=fgd.create_group('beta_se') 
                for indf,f in enumerate(features): fgdbse.create_dataset(f,data=beta_se[indf])
                fgdpo=fgd.create_group('position') 
                for indf,f in enumerate(features): fgdpo.create_dataset(f,data=position[indf].astype(int))
                fgds=fgd.create_group('snp_id') 
                for indf,f in enumerate(features): fgds.create_dataset(f,data=snp_id[indf])
                                                                          
                                                 
                fgs=fg.create_group('summary_data')
                fgs.create_dataset('min_p_value',data= np.nanmin(np.hstack(pv)) [None])
                fgs.create_dataset('min_p_value_beta',data=np.hstack(beta)[np.nanargmin(np.hstack(pv))][None])
                fgs.create_dataset('min_p_value_beta_se',data=np.hstack(beta_se)[np.nanargmin(np.hstack(pv))][None])
                        
                p_bonf=np.nanmin(local_adjustment((pv),method=local_adjustment_method));
                if p_bonf>1:p_bonf=np.array(1)
                fgs.create_dataset('min_p_value_local_adjusted',data=p_bonf[None])
                
                minfeatureindex=np.argmin([np.nanmin(ppv) *len(ppv) for ppv in pv])
                fgs.create_dataset('min_p_value_feature_id',data= np.array(fgm ['feature_id'][minfeatureindex].astype('S'))[None])
        
                min_snp_id=frez[features[minfeatureindex]]['snp_id'][np.nanargmin(pv[minfeatureindex])].astype('U')
                fgs.create_dataset('min_p_value_snp_id',data= np.array(min_snp_id).astype('S')[None])
                fgs.create_dataset('min_p_value_position',data= np.array(fsnp[min_snp_id]['position'])[None])
        
                count+=1

        frez.close()
        fOut.close()

def replication_two_features(path_data =None,  path_data2=None,  traits=None,
    qtl_results_file='qtl_results_',    snp_metadata_file='snp_metadata_',    feature_metadata_file='feature_metadata_',
    results_genome_file='qtl_results_genome',    feature_report='ensembl_gene_id',p_value_field='empirical_feature_p_value',thr=0.01):
    
    _doc=" aggregates qtl results from two traits at feature_report level; return replication of pvalues for trait1  signigicant snps in trait2 "

    if path_data2 is None:
        path_data2 =path_data 
    featureh5=[h5py.File(path_data+'/'+traits[0]+'_'+feature_report+'_'+results_genome_file+'.h5','r'),h5py.File(path_data2+'/'+traits[1]+'_'+feature_report+'_'+results_genome_file+'.h5','r')]
    
    feature_ids=np.intersect1d(list(featureh5[0].keys()),list(featureh5[1].keys()))
 
    rez={}
    rez[p_value_field]=np.zeros(len(feature_ids))+np.nan
    rez['beta']=np.zeros(len(feature_ids))+np.nan
    rez['number_features']=np.zeros(len(feature_ids))+np.nan
    rez['replicated_number_features']=np.zeros(len(feature_ids))+np.nan
    rez['replicated_p_value']=np.zeros(len(feature_ids))+np.nan
    rez['p_value_raw']=np.zeros(len(feature_ids))+np.nan
    rez['replicated_p_value_raw']=np.zeros(len(feature_ids))+np.nan
    rez['replicated_beta']=np.zeros(len(feature_ids))+np.nan
    rez['replicated_self_p_value']=np.zeros(len(feature_ids))+np.nan
    rez['replicated_self_beta']=np.zeros(len(feature_ids))+np.nan
    rez['feature_id']=np.zeros(len(feature_ids),dtype='|S32')
    rez['gene_name']=np.zeros(len(feature_ids),dtype='|S32')
    rez['snp_id']=np.zeros(len(feature_ids),dtype='|S32')
    rez['chromosome']=np.zeros(len(feature_ids),dtype='|S5')
    rez['strand']=np.zeros(len(feature_ids),dtype='|S5')
    rez['position']=np.zeros(len(feature_ids))+np.nan
 
    for indf, feature in enumerate(feature_ids):
        
    
        try:temp=featureh5[0][feature]['summary_data/min_p_value'][0]
        except: continue
        if temp<thr:
            
            rez[p_value_field][indf]=temp
            rez['number_features'][indf]=len(featureh5[0][feature]['data/beta'])
 
            rez['beta'][indf]=featureh5[0][feature]['summary_data/min_p_value_beta'][0]
            rez['feature_id'][indf]=featureh5[0][feature]['summary_data/min_p_value_feature_id'][:][0]
#            print (rez['feature_id'][indf])
            rez['chromosome'][indf]=featureh5[0][feature]['metadata/chromosome'][:][0]
            try: rez['strand'][indf]=featureh5[0][feature]['metadata/feature_strand'][:][0]
            except: 1
            try:rez['gene_name'][indf]=featureh5[0][feature]['metadata/gene_name'][:][0]
            except:1
            
            rez['snp_id'][indf]=featureh5[0][feature]['summary_data/min_p_value_snp_id'][:][0]
            rez['position'][indf]=featureh5[0][feature]['summary_data/min_p_value_position'][:][0]
            
            try:
                temppvraw=[featureh5[0][feature]['data']['p_value_raw'][f1][:][featureh5[0][feature]['data/snp_id'][f1][:]==featureh5[0][feature]['summary_data/min_p_value_snp_id'][0]] \
                    for f1 in featureh5[0][feature]['data/features']]
                rez['p_value_raw'][indf]=np.nanmin( np.hstack(temppvraw))     
            except: 1
            try:
                rez['replicated_p_value_raw'][indf]=np.nanmin( np.hstack([featureh5[1][feature]['data']['p_value_raw'][f1][:][featureh5[1][feature]['data/snp_id'][f1][:]==featureh5[0][feature]['summary_data/min_p_value_snp_id'][0]] for f1 in featureh5[1][feature]['data/features']])     )
            except: 1

            try:
                rez['replicated_p_value'][indf]=np.nanmin(np.hstack([featureh5[1][feature]['data'][p_value_field][f1][:]\
                    [featureh5[1][feature]['data/snp_id'][f1][:]==featureh5[0][feature]['summary_data/min_p_value_snp_id'][0]] for f1 in featureh5[1][feature]['data/features']])     )
                rez['replicated_number_features'][indf]=len(featureh5[1][feature]['data/beta'])
            except: 1
            
            try:
                temp=np.hstack([featureh5[1][feature]['data/beta'][f1][:][featureh5[1][feature]['data/snp_id'][f1][:]==featureh5[0][feature]['summary_data/min_p_value_snp_id'][0]]     for f1 in featureh5[1][feature]['data/features']])
                rez['replicated_beta'][indf]=np.nanmin(temp) if featureh5[0][feature]['summary_data/min_p_value_beta'][0]<0 else np.nanmax(temp)    
            except: 1
            
            try:
                rez['replicated_self_p_value'][indf]=np.min([featureh5[0][feature]['data'][p_value_field][f1][:][featureh5[0][feature]['data/snp_id'][f1][:]==featureh5[0][feature]['summary_data/min_p_value_snp_id'][0]]\
                    for f1 in np.setdiff1d(featureh5[0][feature]['data/features'],featureh5[0][feature]['summary_data/min_p_value_feature_id'][0])]) 
            except: 1
            try:
                temp=[featureh5[0][feature]['data/beta'][f1][:][featureh5[0][feature]['data/snp_id'][f1][:]==featureh5[0][feature]['summary_data/min_p_value_snp_id'][0]]\
                    for f1 in np.setdiff1d(featureh5[0][feature]['data/features'],featureh5[0][feature]['summary_data/min_p_value_feature_id'][0])] 
                rez['replicated_self_beta'][indf]=np.nanmin(temp) if featureh5[0][feature]['summary_data/min_p_value_beta'][0]<0 else np.nanmax(temp)    
            except: 1
            
            if (indf%500)==0:print(feature +'_'+str(indf) +'_'+ str(rez['number_features'][indf])+'_'+ str(rez['replicated_number_features'][indf]))
            
    for f in featureh5: f.close()
    rez['ensembl_gene_id']=feature_ids.astype('U')
    for key in ['feature_id','chromosome','strand','snp_id','gene_name']:
        rez[key]=rez[key].astype('U')
    rez['traits']=traits
    return rez


def replication_two_features_allsnps(path_data =None,  path_data2=None,  traits=None,trait_labels=None,
    qtl_results_file='qtl_results_',    snp_metadata_file='snp_metadata_',    feature_metadata_file='feature_metadata_',
    results_genome_file='qtl_results_genome',    feature_report='ensembl_gene_id',p_value_field='empirical_feature_p_value',thr0=0.001,thr1=0.05):
    
    _doc=" aggregates qtl results from two traits at feature_report level; return replication of pvalues for trait1  signigicant snps in trait2 "

    if path_data2 is None:
        path_data2 =path_data 
    featureh5=[h5py.File(path_data+'/'+traits[0]+'_'+feature_report+'_'+results_genome_file+'.h5','r'),h5py.File(path_data2+'/'+traits[1]+'_'+feature_report+'_'+results_genome_file+'.h5','r')]
    
    feature_ids=np.intersect1d(list(featureh5[0].keys()),list(featureh5[1].keys()))
    rez={}
    for i in range(2):
        rez[i]={}
    
        for key in ['snp_id','beta',p_value_field]:
            rez[i][key]=[featureh5[i][feature]['data'][key][list(featureh5[i][feature]['data']['snp_id'].keys())[0]][:] for indf, feature in enumerate(feature_ids)]
            
    
    key='snp_id'
    snps=[np.intersect1d(rez[0][key][indf],rez[1][key][indf]) for indf, feature in enumerate(feature_ids)]
    features=[np.repeat(feature,len(snps[indf])) for indf, feature in enumerate(feature_ids)]
    gene=[np.repeat(featureh5[0][feature]['metadata/gene_name'][:].astype('U')[0],len(snps[indf])) for indf, feature in enumerate(feature_ids)]
    
    
    df=pd.DataFrame(data=np.vstack([np.hstack(features).astype('U'),np.hstack(snps).astype('U'),np.hstack( gene).astype('U')]).T,columns=['feature_id','snps','gene_name'])
    
    for i in range(2):
        for key in ['snp_id']:
            df[key+'_'+trait_labels[i]]=np.hstack([rez[i][key][indf][np.in1d(rez[i]['snp_id'][indf],snps[indf])] for indf, feature in enumerate(feature_ids)]).astype('U')
       
    for i in range(2):
        for key in ['beta',p_value_field]:
            df[key+'_'+trait_labels[i]]=np.hstack([rez[i][key][indf][np.in1d(rez[i]['snp_id'][indf],snps[indf])] for indf, feature in enumerate(feature_ids)]).astype(float)
        
 
    only0=df.iloc[np.where((df[p_value_field+'_'+trait_labels[0]].values<thr0)&(df[p_value_field+'_'+trait_labels[1]].values>thr1))[0]].set_index('feature_id',drop=0)
    only0=pd.concat([only0,only0])
    
    df.to_csv(path_or_buf=path_data+traits[0]+'_'+traits[1]+'_qtl_results_allsnps.txt',mode='w', sep='\t', columns=None, header=True, index=True)
    np.unique(only0['feature_id']).shape
             
    
    df1= df.iloc[np.where((df[p_value_field+'_'+trait_labels[0]].values<thr0)|(df[p_value_field+'_'+trait_labels[1]].values<thr1))[0]]
    df1.to_csv(path_or_buf=path_data+traits[0]+'_'+traits[1]+'_qtl_results_relevantsnps.txt',mode='w', sep='\t', columns=None, header=True, index=True)
    np.unique(only0['feature_id']).shape
    
    df2=[only0.loc[f].iloc[np.argmin(only0.loc[f][p_value_field+'_'+trait_labels[0]].values)] for f in np.unique(only0['feature_id'])]
    df2=pd.concat(df2,1).transpose()
    df2['number_snps']=np.array([only0.loc[f].shape[0] for f in np.unique(only0['feature_id'])])
     
    df2.to_csv(path_or_buf=path_data+traits[0]+'_'+traits[1]+'_qtl_results_ONLYPROTEIN.txt',mode='w', sep='\t', columns=None, header=True, index=True)
     
          
    return df

 
