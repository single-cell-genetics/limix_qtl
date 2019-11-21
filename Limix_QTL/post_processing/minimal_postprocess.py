import h5py
import glob
import os.path
import numpy as np
import pandas as pd
import argparse
import glob
import pdb

def minimal_qtl_processing(QTL_Dir, OutputDir, writeToOneFile=True, compressed = False, overWrite=True, minimalPValue = 1, minimalFeaturePValue = 1, topMode = False, debugMode = False):
    qtl_results_file='qtl_results_'
    snp_metadata_file='snp_metadata_'
    feature_metadata_file='feature_metadata_'
    if topMode:
        output_file='top_qtl_results_'
    else:
        output_file='qtl_results_'
    
    
    h5FilesToProcess = (glob.glob(QTL_Dir+"/qtl_*.h5"))
    
    #print(h5FilesToProcess)
    #print(os.path.dirname(h5FilesToProcess[1]))
    for file in h5FilesToProcess :
        partTmp = os.path.basename(file).replace(qtl_results_file,"").replace(".h5","")
        if(debugMode):
            print(partTmp)
        if(writeToOneFile): 
            outputFile = OutputDir+output_file+"all.txt"
        else:
            outputFile = OutputDir+output_file+partTmp+".txt"
        
        #print(outputFile)
        if(((os.path.isfile(outputFile) or os.path.isfile(outputFile+".gz")) and not overWrite) and not writeToOneFile):
            #print("Skipping: "+partTmp)
            continue
        #else :
            #print('Processing: '+partTmp)
        #print(partTmp)
        if not os.path.isfile(QTL_Dir+"/"+snp_metadata_file+partTmp+".txt"):
            print("Skipping: " +partTmp + " not all necessary files are present.")
            continue
        if not os.path.isfile(QTL_Dir+"/"+feature_metadata_file+partTmp+".txt"):
            print("Skipping: " +partTmp + " not all necessary files are present.")
            continue
        try :
            #print(QTL_Dir+"/"+feature_metadata_file+partTmp+".txt")
            #print(QTL_Dir+"/"+snp_metadata_file+partTmp+".txt")
            ffea= pd.read_table(QTL_Dir+"/"+feature_metadata_file+partTmp+".txt", sep='\t')
            fsnp= pd.read_table(QTL_Dir+"/"+snp_metadata_file+partTmp+".txt", sep='\t')
        except:
            print("Issue in features or snp annotation.\n Skipping: "+partTmp)
            continue
        
        ffea = ffea.rename(index=str, columns={"chromosome": "feature_chromosome", "start": "feature_start", "end": "feature_end"})
        fsnp = fsnp.rename(index=str, columns={"chromosome": "snp_chromosome", "position": "snp_position"})
        #pdb.set_trace()
        frez=h5py.File(file,'r')
        frezkeys= np.array([k.replace('_i_','') for k in list(frez.keys())])
        
        data={}
        for key in ['feature_id','snp_id','p_value','beta','beta_se','empirical_feature_p_value']:
            data[key]=np.zeros(len(np.unique(list(frezkeys))),dtype='object')+np.nan
        
        for ifea,report_feature in enumerate(np.unique(list(frezkeys))):
            for key in ['snp_id','p_value','beta','beta_se','empirical_feature_p_value']:
                temp = np.array(frez[report_feature][key])
                data[key][ifea]=np.hstack(temp).astype('U')
            data['feature_id'][ifea]=np.hstack(np.repeat(report_feature,len(frez[report_feature][key])))
        #pdb.set_trace()
        for key in data.keys():
            data[key]=np.hstack(data[key])
        
        temp=pd.DataFrame(data)
        #print(temp.head())
        
        temp = pd.merge(temp, ffea, on='feature_id', how='left')
        #print(temp.head())
        
        #pdb.set_trace()
        if(len(glob.glob(QTL_Dir+'snp_qc_metrics_naContaining_feature_*.txt'))>0):
            ##Here we need to check, we can based on the output of glob do this quicker.
            temp2 = pd.DataFrame(columns=temp.columns)
            for key in frezkeys:
                if os.path.isfile(QTL_Dir+"/snp_qc_metrics_naContaining_feature_"+key+".txt"):
                    fsnp_rel = pd.read_table(QTL_Dir+"/snp_qc_metrics_naContaining_feature_"+key+".txt", sep='\t')
                    temp_t = temp.loc[temp["feature_id"]==key]
                    fsnp_t = fsnp.loc[:,["snp_id","snp_chromosome","snp_position","assessed_allele"]]
                    fsnp_t = pd.merge(fsnp_t, fsnp_rel, on='snp_id', how='right')
                    temp_t = pd.merge(temp_t, fsnp_t, on='snp_id', how='left')
                    temp2 = temp2.append(temp_t,sort=False)
                else:
                    temp_t = temp.loc[temp["feature_id"]==key]
                    temp_t = pd.merge(temp_t, fsnp, on='snp_id', how='left')
                    temp2 = temp2.append(temp_t,sort=False)
                data[key]=np.zeros(len(np.unique(list(frezkeys))),dtype='object')+np.nan
            temp = temp2
            temp2 = None
        else : 
            temp = pd.merge(temp, fsnp, on='snp_id', how='left')
        
        #pdb.set_trace()
        #print(temp.head())
        temp['empirical_feature_p_value'] = temp['empirical_feature_p_value'].astype(float)
        temp['p_value'] = temp['p_value'].astype(float)
        if(minimalPValue<1):
            temp=pd.DataFrame(data).iloc[data['p_value'].astype(float)<minimalPValue]
        
        if(minimalFeaturePValue<1):
            temp=pd.DataFrame(data).iloc[data['empirical_feature_p_value'].astype(float)<minimalFeaturePValue]
        
        temp = temp.sort_values(by =['empirical_feature_p_value',"p_value"], ascending=[True,True])
        
        if topMode:
            temp = temp.groupby(temp['feature_id']).first()
            temp['feature_id'] = temp.index
        #print(outputFile)#
        #temp.to_csv(path_or_buf=outputFile, mode='w', sep='\t', columns=None,index=None)
        #print('w'if not os.path.isfile(outputFile) else 'a')
        if(not compressed):
            temp.to_csv(path_or_buf=outputFile, mode='w'if not os.path.isfile(outputFile) else 'a', sep='\t', columns=None,index=None, header= True if not os.path.isfile(outputFile) else False)
        else:
            temp.to_csv(path_or_buf=outputFile+".gz", mode='w'if not os.path.isfile(outputFile) else 'a', sep='\t', columns=None,index=None,compression='gzip', header= True if not os.path.isfile(outputFile) else False )

def parse_args():
    parser = argparse.ArgumentParser(description='Run QTL analysis given genotype, phenotype, and annotation.')
    parser.add_argument('--input_dir','-id',required=True)
    parser.add_argument('--ouput_dir','-od',required=True)
    parser.add_argument('--single_file_output','-sfo', action="store_true", required=False,default=False)
    parser.add_argument('--write_compressed','-wc', action="store_true", required=False,default=False)
    parser.add_argument('--overwrite_old','-oo', action="store_true", required=False)
    parser.add_argument('--minimimal_reporting_p','-mrp', required=False, default=1.0)
    parser.add_argument('--minimal_reporting_featureP', '-mrf', required=False, default=1.0)
    parser.add_argument('--top_feature_based', '-tfb', action="store_true", required=False, default=False)
    parser.add_argument('--debug', '-d', action="store_true", required=False, default=False)
    args = parser.parse_args()
    return args

if __name__=='__main__':
    args = parse_args()
    inputDir  = args.input_dir
    outputDir = args.ouput_dir
    writeToOneFile = args.single_file_output
    compressed = args.write_compressed
    overWrite = args.overwrite_old
    minimalPValue = args.minimimal_reporting_p
    minimalFeaturePValue = args.minimal_reporting_featureP
    topMode = args.top_feature_based
    debugMode = args.debug

    minimal_qtl_processing(inputDir, outputDir, writeToOneFile, compressed, overWrite, float(minimalPValue), float(minimalFeaturePValue),topMode, debugMode)

