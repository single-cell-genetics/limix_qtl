import sys
import os
import tables
import numpy as np
import qtl_fdr_utilities

#V0.1.1

class hdf5_writer:

    def __init__(self,output_filename):
        self.h5file = tables.open_file(output_filename,'w')

    def close(self):
        self.h5file.close()

    def add_result_df(self,qtl_results_df):
        assert(len(set(qtl_results_df['feature_id'].values))==1)
        feature_id = qtl_results_df['feature_id'].values[0]
        column_names = ['snp_id','p_value','beta','beta_se','empirical_feature_p_value']
        try:
            #get the existing table for this feature
            table = self.h5file.get_node('/'+feature_id)
        except tables.exceptions.NoSuchNodeError:
            #this table doesn't exist yet - create it
            table = self.h5file.create_table(self.h5file.root,feature_id,QTL_result_hdf5,"QTL analysis results")
            pass
        qtl_result = table.row
        for idx,df_row in qtl_results_df.iterrows():
            for col_name in column_names:
                qtl_result[col_name] = df_row[col_name]
            qtl_result.append()
        table.flush()

    def apply_pval_correction(self,feature_id,top_pvalues_perm,cis_mode):
        '''Function to correct p values based on nominal p values and the top
        hits from permutation runs for the given feature.'''
        table = self.h5file.get_node('/'+feature_id)
        if(np.mean(top_pvalues_perm)>=0.999999999 and np.var(top_pvalues_perm)==0):
            for row in table:
                row['empirical_feature_p_value'] = row['p_value']
                row.update()
            alpha_para=-9
            beta_para=-9
        else:
            correction_function, alpha_para, beta_para = qtl_fdr_utilities.define_correction_function(top_pvalues_perm,cis_mode)
            for row in table:
                row['empirical_feature_p_value'] = correction_function(row['p_value'])
                row.update()
        table.flush()
        return [alpha_para, beta_para]


class text_writer:

    def __init__(self,output_filename):
        self.column_names = ['feature_id','snp_id','p_value','beta','beta_se','empirical_feature_p_value']
        with open(output_filename,'w') as f:
            header = '\t'.join(self.column_names)
            f.write(header+'\n')
        self.outfile = open(output_filename,'a')

    def close(self):
        self.outfile.close()

    def add_result_df(self,qtl_results_df):
        qtl_results_df.loc[:,self.column_names].to_csv(self.outfile,header=None,mode='a',index=False,sep='\t')


class hdf5_permutations_writer:

    def __init__(self,output_filename,n_permutations):
        self.h5file = tables.open_file(output_filename,'w',max_columns=(n_permutations+2))
        self.column_names = ['snp_id'] + ['permutation_'+str(x) for x in range(n_permutations)]
        #define the permutation result object on-the-fly, depending on the number of permutations that will be performed
        self.permutation_result_definition = dict([(x,tables.Float64Col()) for x in self.column_names if x.split('_')[0]=='permutation'])
        self.permutation_result_definition['snp_id'] = tables.StringCol(100)

    def close(self):
        self.h5file.close()

    def add_permutation_results_df(self,permutation_results_df,feature_id):
        '''Takes as input permutation_results_df and feature_id.
        permutation_results_df must contain a "snp_id" column, and
        columns labelled ""permutation_1","permutation_2",...,"permutation_n",
        where n=the number of permutations specified when initialising the
        hdf5_permutations_writer.'''
        try:
            #get the existing table for this feature
            table = self.h5file.get_node('/'+feature_id)
        except tables.exceptions.NoSuchNodeError:
            #this table doesn't exist yet - create it
            table = self.h5file.create_table(self.h5file.root,
                                             feature_id,
                                             self.permutation_result_definition,
                                             "Permutation analysis results")
            pass
        permutation_result = table.row
        for idx,df_row in permutation_results_df.iterrows():
            for col_name in self.column_names:
                permutation_result[col_name] = df_row[col_name]
            permutation_result.append()
        table.flush()


class QTL_result_hdf5(tables.IsDescription):
    snp_id  = tables.StringCol(100)   # 100-character String
    p_value = tables.Float64Col()    # double (double-precision)
    beta = tables.Float64Col()    # double (double-precision)
    beta_se = tables.Float64Col()    # double (double-precision)
    empirical_feature_p_value = tables.Float64Col()
