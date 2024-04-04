import h5py
import glob
import os.path
import numpy as np
import pandas as pd
import argparse
import glob
import pdb
import qtl_fdr_utilities

#writeToOneFile=True; compressed = False; overWrite=True; minimalPValue = 1; minimalFeaturePValue = 1; topMode = False; debugMode = False
#QTL_Dir = "./"; OutputDir = "./"; featureGroupFile="./featureGrouping.txt"

def minimal_qtl_processing(QTL_Dir, OutputDir, featureGroupFile, cis_mode, writeToOneFile=True, compressed = False, overWrite=True, minimalPValue = 1, minimalFeaturePValue = 1, topMode = False, topGroupMode = False , debugMode = False):
    qtl_results_file='qtl_results_'
    snp_metadata_file='snp_metadata_'
    feature_metadata_file='feature_metadata_'
    if topMode:
        output_file='top_qtl_results_'
    elif topGroupMode : 
        output_file='top_group_qtl_results_'
    else:
        output_file='qtl_results_'

    if os.path.isfile(featureGroupFile):
        feature_grouping_df = pd.read_csv(featureGroupFile,sep='\t')
    else :
        print("Error: feature grouping file doesn't excist.")
        return(-1)
    
    h5FilesToProcess = (glob.glob(QTL_Dir+"/qtl_*.h5"))

    #iterate over h5files
    #print(h5FilesToProcess)
    #print(os.path.dirname(h5FilesToProcess[1]))
    for file in h5FilesToProcess :
        print(file)
        
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

        data=pd.DataFrame(data)
        
        data = pd.merge(data, feature_grouping_df, on='feature_id', how='left')
        
        #print(data.head())
        data = pd.merge(data, ffea, on='feature_id', how='left')
        
        #print(data.head())
        
        if(len(glob.glob(QTL_Dir+'snp_qc_metrics_naContaining_feature_*.txt'))>0):
            ##Here we need to check, we can based on the output of glob do this quicker.
            temp2 = pd.DataFrame(columns=data.columns)
            for key in frezkeys:
                if os.path.isfile(QTL_Dir+"/snp_qc_metrics_naContaining_feature_"+key+".txt"):
                    fsnp_rel = pd.read_table(QTL_Dir+"/snp_qc_metrics_naContaining_feature_"+key+".txt", sep='\t')
                    temp_t = data.loc[data["feature_id"]==key]
                    fsnp_t = fsnp.loc[:,["snp_id","snp_chromosome","snp_position","assessed_allele"]]
                    fsnp_t = pd.merge(fsnp_t, fsnp_rel, on='snp_id', how='right')
                    temp_t = pd.merge(temp_t, fsnp_t, on='snp_id', how='left')
                    temp2 = temp2.append(temp_t,sort=False)
                else:
                    temp_t = data.loc[data["feature_id"]==key]
                    temp_t = pd.merge(temp_t, fsnp, on='snp_id', how='left')
                    temp2 = temp2.append(temp_t,sort=False)
                data[key]=np.zeros(len(np.unique(list(frezkeys))),dtype='object')+np.nan
            data = temp2
            temp2 = None
        else :
            data = pd.merge(data, fsnp, on='snp_id', how='left')
        
        ##make sure the pvalues are seen as floats.
        data['empirical_feature_p_value'] = data['empirical_feature_p_value'].astype(float)
        data['p_value'] = data['p_value'].astype(float)
        
        for featGroup in np.unique(data["feature_group"]) :
            #pdb.set_trace()
            if(debugMode):
                print(featGroup)
            selection = np.unique(data.loc[data["feature_group"] == featGroup, "feature_id"])
            permutationValues = None
            for relEntry in selection : 
                if os.path.isfile(QTL_Dir+"/Permutation.pValues."+relEntry+".txt"):
                    pValTmp = pd.read_csv(QTL_Dir+"/Permutation.pValues."+relEntry+".txt",sep='\t',header=None)
                    if permutationValues is not None:
                        permutationValues[0][np.where(permutationValues[0] > pValTmp[0])[0]] = pValTmp[0][np.where(permutationValues[0] > pValTmp[0])[0]]
                    else :
                        permutationValues = pValTmp
                else :
                    print("Warning, relevant permutation pValue file not observed!")
            ##Here we need to refit the beta distribution.
            correction_function, alpha_para, beta_para = qtl_fdr_utilities.define_correction_function(permutationValues[0],cis_mode)
            for relEntry in selection : 
                #Replace empirical_feature_p_value, alpha and beta param for the relevant feature.
                data.loc[np.where(data["feature_id"]==relEntry)[0],"empirical_feature_p_value"] = correction_function(data.loc[np.where(data["feature_id"]==relEntry)[0],"p_value"])
                data.loc[np.where(data["feature_id"]==relEntry)[0],"alpha_param"] = alpha_para
                data.loc[np.where(data["feature_id"]==relEntry)[0],"beta_param"] = beta_para
            
        
        #print(data.head())
        
        if(minimalPValue<1):
            data=pd.DataFrame(data).iloc[data['p_value'].astype(float)<minimalPValue]

        if(minimalFeaturePValue<1):
            data=pd.DataFrame(data).iloc[data['empirical_feature_p_value'].astype(float)<minimalFeaturePValue]

        data = data.sort_values(by =['empirical_feature_p_value',"p_value"], ascending=[True,True])

        if topMode:
            data = data.groupby(data['feature_id']).first()
            data['feature_id'] = data.index
        elif topGroupMode:
            data = data.groupby(data['feature_group']).first()
            data['feature_group'] = data.index
        #print(outputFile)#
        #data.to_csv(path_or_buf=outputFile, mode='w', sep='\t', columns=None,index=None)
        #print('w'if not os.path.isfile(outputFile) else 'a')
        outputFileExists = os.path.isfile(outputFile)
        writeMode = 'w' if not outputFileExists else 'a'
        writeHeader = True if not outputFileExists else False
        if(compressed):
            temp.to_csv(path_or_buf=outputFile, mode=writeMode, sep='\t', columns=None,index=None, header=writeHeader,compression='gzip')
        else:
            temp.to_csv(path_or_buf=outputFile, mode=writeMode, sep='\t', columns=None,index=None, header=writeHeader)
        wroteData = True

def parse_args():
    parser = argparse.ArgumentParser(description='Run QTL analysis given genotype, phenotype, and annotation.')
    parser.add_argument('--input_dir','-id',required=True)
    parser.add_argument('--ouput_dir','-od',required=True)
    parser.add_argument('--feature_grouping_file','-fgf',required=True)
    parser.add_argument('--cis','-c', action="store_true")
    parser.add_argument('--single_file_output','-sfo', action="store_true", required=False,default=False)
    parser.add_argument('--write_compressed','-wc', action="store_true", required=False,default=False)
    parser.add_argument('--overwrite_old','-oo', action="store_true", required=False)
    parser.add_argument('--minimimal_reporting_p','-mrp', required=False, default=1.0)
    parser.add_argument('--minimal_reporting_featureP', '-mrf', required=False, default=1.0)
    parser.add_argument('--top_feature_based', '-tfb', action="store_true", required=False, default=False)
    parser.add_argument('--top_feature_group_based', '-tfgb', action="store_true", required=False, default=False)
    parser.add_argument('--debug', '-d', action="store_true", required=False, default=False)
    args = parser.parse_args()
    return args

if __name__=='__main__':
    args = parse_args()
    inputDir  = args.input_dir
    outputDir = args.ouput_dir
    featureGroupingFile = args.feature_grouping_file
    cis = args.cis
    writeToOneFile = args.single_file_output
    compressed = args.write_compressed
    overWrite = args.overwrite_old
    minimalPValue = args.minimimal_reporting_p
    minimalFeaturePValue = args.minimal_reporting_featureP
    topMode = args.top_feature_based
    topGroupMode = args.top_feature_group_based
    debugMode = args.debug
    
    if topMode & topGroupMode:
        print("Error, cant do both group and feature top mode.")
        exit()
    
    minimal_qtl_processing(inputDir, outputDir, featureGroupingFile, cis, writeToOneFile, compressed, overWrite, float(minimalPValue), float(minimalFeaturePValue),topMode, topGroupMode, debugMode)
