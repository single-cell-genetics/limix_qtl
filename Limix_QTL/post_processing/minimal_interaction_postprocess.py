import h5py
import glob
import os.path
import numpy as np
import pandas as pd
import argparse
import glob
import pdb
import pickle
from pathlib import Path
import re

def minimal_iqtl_processing(QTL_Dir, OutputDir, writeToOneFile=True, compressed = False, overWrite=True, minimalPValue = 1, minimalFeaturePValue = 1, topMode = False, debugMode = False):
    qtl_results_file='iqtl_results_'
    snp_metadata_file='snp_metadata_'
    feature_metadata_file='feature_metadata_'
    if topMode:
        output_file='top_iqtl_results_'
    else:
        output_file='iqtl_results_'
    
    if not OutputDir.endswith(os.path.sep):
        OutputDir += os.path.sep
    
    allH5Files = (glob.glob(QTL_Dir+"/"+qtl_results_file+"*.h5"))

    pattern = re.compile(r"^(?!.*_snp_)(?!.*_cov_).*$")
    h5FilesToProcess = [i for i in allH5Files if pattern.match(i)]

    if(len(h5FilesToProcess)*3!=len(allH5Files)):
        print("There is a mismatch between the number of h5 files and the expected number of h5 files.\n    Please check your input folder. ")
        exit

    #iterate over h5files
    #print(h5FilesToProcess)
    #print(os.path.dirname(h5FilesToProcess[1]))

    wroteData = False
    for file in h5FilesToProcess :
        #print(file)
        partTmp = os.path.basename(file).replace(qtl_results_file,"").replace(".h5","")
        if(debugMode):
            print(partTmp)
        if(writeToOneFile):
            outputFile = OutputDir+output_file+"all.txt"
        else:
            outputFile = OutputDir+output_file+partTmp+".txt"

        if compressed:
            outputFile = outputFile+".gz"
        
        #print(outputFile)
        if((os.path.isfile(outputFile) and not overWrite) and not writeToOneFile):
            #print("Skipping: "+partTmp)
            continue
        #else :
            #print('Processing: '+partTmp)
        #print(partTmp)
        if not os.path.isfile(QTL_Dir+"/"+snp_metadata_file+partTmp+".txt") and not os.path.isfile(QTL_Dir+"/"+snp_metadata_file+partTmp+".txt.gz"):
            print("Skipping: " +partTmp + " not all necessary files are present.")
            continue

        try :
            #print(QTL_Dir+"/"+feature_metadata_file+partTmp+".txt")
            ffea= pd.read_table(QTL_Dir+"/"+feature_metadata_file+partTmp+".txt", sep='\t')
        except :
            try :
                ffea= pd.read_table(QTL_Dir+"/"+feature_metadata_file+partTmp+".txt.gz", sep='\t')
            except :
                print("Issue in features annotation.\n Skipping: "+partTmp)
                continue
                
        if not os.path.isfile(QTL_Dir+"/"+feature_metadata_file+partTmp+".txt") and not os.path.isfile(QTL_Dir+"/"+feature_metadata_file+partTmp+".txt.gz"):
            print("Skipping: " +partTmp + " not all necessary files are present.")
            continue
        try :
            #print(QTL_Dir+"/"+snp_metadata_file+partTmp+".txt")
            fsnp= pd.read_table(QTL_Dir+"/"+snp_metadata_file+partTmp+".txt", sep='\t')
        except:
            try:
                fsnp= pd.read_table(QTL_Dir+"/"+snp_metadata_file+partTmp+".txt.gz", sep='\t')
            except :
                print("Issue in snp annotation.\n Skipping: "+partTmp)
                continue

        ffea = ffea.rename(index=str, columns={"chromosome": "feature_chromosome", "start": "feature_start", "end": "feature_end"})
        fsnp = fsnp.rename(index=str, columns={"chromosome": "snp_chromosome", "position": "snp_position"})

        ##Interactions
        frez=h5py.File(file,'r')
        frezkeys= np.array([k.replace('_i_','') for k in list(frez.keys())])

        data={}
        for key in ['feature_id','snp_id','p_value','beta','beta_se','empirical_feature_p_value']:
            data[key]=np.zeros(len(np.unique(list(frezkeys))),dtype='object')+np.nan

        for ifea,report_feature in enumerate(np.unique(list(frezkeys))):
            #if(debugMode):
                #print(report_feature)
            for key in ['snp_id','p_value','beta','beta_se','empirical_feature_p_value']:
                temp = np.array(frez[report_feature][key])
                data[key][ifea]=np.hstack(temp).astype('U')
            data['feature_id'][ifea]=np.hstack(np.repeat(report_feature,len(frez[report_feature][key])))
        for key in data.keys():
            data[key]=np.hstack(data[key])

        pattern = re.compile(r"iqtl_results_")

        ##Covariate 
        frez=h5py.File(pattern.sub("iqtl_results_cov_", file),'r')
        frezkeys= np.array([k.replace('_i_','') for k in list(frez.keys())])

        dataCov={}
        for key in ['feature_id','snp_id','beta','beta_se']:
            dataCov[key]=np.zeros(len(np.unique(list(frezkeys))),dtype='object')+np.nan

        for ifea,report_feature in enumerate(np.unique(list(frezkeys))):
            #if(debugMode):
                #print(report_feature)
            for key in ['snp_id','beta','beta_se']:
                temp = np.array(frez[report_feature][key])
                dataCov[key][ifea]=np.hstack(temp).astype('U')
            dataCov['feature_id'][ifea]=np.hstack(np.repeat(report_feature,len(frez[report_feature][key])))
        for key in dataCov.keys():
            dataCov[key]=np.hstack(dataCov[key])
            
        ##SNP 
        frez=h5py.File(pattern.sub("iqtl_results_snp_", file),'r')
        frezkeys= np.array([k.replace('_i_','') for k in list(frez.keys())])

        dataSnp={}
        for key in ['feature_id','snp_id','beta','beta_se']:
            dataSnp[key]=np.zeros(len(np.unique(list(frezkeys))),dtype='object')+np.nan

        for ifea,report_feature in enumerate(np.unique(list(frezkeys))):
            #if(debugMode):
                #print(report_feature)
            for key in ['snp_id','beta','beta_se']:
                temp = np.array(frez[report_feature][key])
                dataSnp[key][ifea]=np.hstack(temp).astype('U')
            dataSnp['feature_id'][ifea]=np.hstack(np.repeat(report_feature,len(frez[report_feature][key])))
        for key in dataSnp.keys():
            dataSnp[key]=np.hstack(dataSnp[key])

        temp=pd.DataFrame(data)

        temp = pd.merge(temp, ffea, on='feature_id', how='left')
        del data
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
                    temp2 = pd.concat([temp2,temp_t]) #temp2.append(temp_t,sort=False)
                else:
                    temp_t = temp.loc[temp["feature_id"]==key]
                    temp_t = pd.merge(temp_t, fsnp, on='snp_id', how='left')
                    temp2 = pd.concat([temp2,temp_t]) #temp2.append(temp_t,sort=False)
                #data[key]=np.zeros(len(np.unique(list(frezkeys))),dtype='object')+np.nan ##
            temp = temp2
            del temp2
        elif(os.path.exists(QTL_Dir+'na_snp_qc_metrics_features_'+partTmp+'.pkl')):
            naQcInfo = None
            with open(QTL_Dir+'na_snp_qc_metrics_features_'+partTmp+'.pkl', 'rb') as f:
                naQcInfo = pickle.load(f)
            
            temp2 = pd.DataFrame(columns=temp.columns)
            for key in frezkeys:
                if key in naQcInfo.keys():
                    temp_t = temp.loc[temp["feature_id"]==key]
                    fsnp_t = fsnp.loc[:,["snp_id","snp_chromosome","snp_position","assessed_allele"]]
                    fsnp_t = pd.merge(fsnp_t, naQcInfo[key], on='snp_id', how='right')
                    temp_t = pd.merge(temp_t, fsnp_t, on='snp_id', how='left')
                    temp2 = pd.concat([temp2,temp_t])
                else:
                    temp_t = temp.loc[temp["feature_id"]==key]
                    temp_t = pd.merge(temp_t, fsnp, on='snp_id', how='left')
                    temp2 = pd.concat([temp2,temp_t])
                data[key]=np.zeros(len(np.unique(list(frezkeys))),dtype='object')+np.nan
            del temp2
        else :
            temp = pd.merge(temp, fsnp, on='snp_id', how='left')
        
        temp = pd.merge(temp, pd.DataFrame(dataCov), on=['feature_id', 'snp_id'], how='outer',suffixes=('','_COV'))
        temp = pd.merge(temp, pd.DataFrame(dataSnp), on=['feature_id', 'snp_id'], how='outer',suffixes=('','_SNP'))
        del dataCov, dataSnp
        #print(temp.head())
        temp['empirical_feature_p_value'] = temp['empirical_feature_p_value'].astype(float)
        temp['p_value'] = temp['p_value'].astype(float)
        if(minimalPValue<1):
            temp=temp.loc[temp['p_value']<minimalPValue]

        if(minimalFeaturePValue<1):
            temp=temp.iloc[temp.loc['empirical_feature_p_value']<minimalFeaturePValue]

        temp = temp.sort_values(by =['empirical_feature_p_value',"p_value"], ascending=[True,True])

        if topMode:
            temp = temp.groupby(temp['feature_id']).first()
            temp['feature_id'] = temp.index
        #print(outputFile)#
        #temp.to_csv(path_or_buf=outputFile, mode='w', sep='\t', columns=None,index=None)
        #print('w'if not os.path.isfile(outputFile) else 'a')
        outputFileExists = os.path.isfile(outputFile)
        writeMode = 'w' if not outputFileExists else 'a'
        writeHeader = True if not outputFileExists else False
        if(compressed):
            temp.to_csv(path_or_buf=outputFile, mode=writeMode, sep='\t', columns=None,index=None, header=writeHeader,compression='gzip')
        else:
            temp.to_csv(path_or_buf=outputFile, mode=writeMode, sep='\t', columns=None,index=None, header=writeHeader)
        wroteData = True
    
    if writeToOneFile and not wroteData:
        if(compressed):
            Path(OutputDir+output_file+"all.txt.gz").touch()
        else:
            Path(OutputDir+output_file+"all.txt").touch()

def parse_args():
    parser = argparse.ArgumentParser(description='Process the h5 and annoation files to flat files after ieQTL mapping.')
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

    minimal_iqtl_processing(inputDir, outputDir, writeToOneFile, compressed, overWrite, float(minimalPValue), float(minimalFeaturePValue),topMode, debugMode)
