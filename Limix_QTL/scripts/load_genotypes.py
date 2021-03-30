import limix.io
import numpy as np

def load_genotypes_plink(plink_prefix):
    '''Return the genotype from a set of plink files
    
    >>> plink_prefix = '../data/geuvadis_CEU_YRI_test_data/Geuvadis_chr1'
    >>> bim,fam,bed,genotypes = load_genotypes_plink(plink_prefix)
    >>> genotypes[0,4]
    2.0
    >>> genotypes[1,-1]
    2.0
    >>> genotypes[-1,-4]
    1.0

    '''
    (bim, fam, bed) = limix.io.plink.read(plink_prefix,verbose=False)

    #make missing values 0 for now (i.e. homozygous for second allele)
    #needs to be replaced with a snp-specific value
    missing_value = 2
    
    def bed_to_genotype_value(input_value):
    	value_mapping = [0,1,2,missing_value]
    	return value_mapping[input_value]
    
    #(nN,nS)
    genotypes = np.array(np.vectorize(bed_to_genotype_value)(bed.compute()),dtype=float).transpose()
    
    return bim,fam,bed,genotypes

if __name__ == "__main__":
    import doctest
    doctest.testmod()
