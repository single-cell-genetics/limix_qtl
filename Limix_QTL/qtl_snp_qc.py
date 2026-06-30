import pandas as pd
import numpy as np
import math

#V0.1.2

def hwe_exact_p(n_hom_ref, n_het, n_hom_alt):
    """
    Pure-Python exact Hardy-Weinberg P value
    (Wigginton et al. 2005).
    """
    obs_hets = int(n_het)
    obs_hom_ref = int(n_hom_ref)
    obs_hom_alt = int(n_hom_alt)

    genotypes = obs_hom_ref + obs_hets + obs_hom_alt
    if genotypes == 0:
        return np.nan

    rare = 2 * min(obs_hom_ref, obs_hom_alt) + obs_hets

    mid = (rare * (2 * genotypes - rare)) // (2 * genotypes)
    if (mid % 2) != (rare % 2):
        mid += 1

    het_probs = np.zeros(rare + 1, dtype=np.float64)
    het_probs[mid] = 1.0
    sum_p = 1.0

    # Left tail
    curr_hets = mid
    curr_hom_r = (rare - mid) // 2
    curr_hom_c = genotypes - mid - curr_hom_r
    while curr_hets >= 2:
        het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (
            4.0 * (curr_hom_r + 1.0) * (curr_hom_c + 1.0)
        )
        sum_p += het_probs[curr_hets - 2]
        curr_hom_r += 1
        curr_hom_c += 1
        curr_hets -= 2

    # Right tail
    curr_hets = mid
    curr_hom_r = (rare - mid) // 2
    curr_hom_c = genotypes - mid - curr_hom_r
    while curr_hets <= rare - 2:
        het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_hom_r * curr_hom_c / (
            (curr_hets + 2.0) * (curr_hets + 1.0)
        )
        sum_p += het_probs[curr_hets + 2]
        curr_hom_r -= 1
        curr_hom_c -= 1
        curr_hets += 2

    het_probs /= sum_p
    threshold = het_probs[obs_hets]
    p_val = het_probs[het_probs <= threshold].sum()
    return p_val

def do_snp_qc(snp_df, min_call_rate, min_maf, min_hwe_P):
    arr = snp_df.to_numpy(dtype=np.float64, copy=False)
    n_samples, n_snps = arr.shape
    names = snp_df.columns

    # -------------------------
    # CALL RATE
    # -------------------------
    valid = np.sum(~np.isnan(arr), axis=0)
    call_rate = valid / float(n_samples)
    fail_cr = call_rate < min_call_rate

    failed_snp_names = list(names[fail_cr])
    keep_cr = ~fail_cr

    call_rate_df = pd.DataFrame(call_rate, index=names, columns=['call_rate'])

    if not np.any(keep_cr):
        return names[keep_cr], failed_snp_names, call_rate_df, None, None

    arr_cr = arr[:, keep_cr]
    names_cr = names[keep_cr]

    # -------------------------
    # GENOTYPE COUNTS (vectorized)
    # -------------------------
    G = np.rint(arr_cr)
    n0 = np.sum(G == 0, axis=0)
    n1 = np.sum(G == 1, axis=0)
    n2 = np.sum(G == 2, axis=0)

    # -------------------------
    # MAF (vectorized)
    # -------------------------
    gc = n0 + n1 + n2
    minor_copies = np.minimum(2*n2 + n1, 2*n0 + n1)

    with np.errstate(divide='ignore', invalid='ignore'):
        maf = minor_copies / (2.0 * gc)

    maf_df = pd.DataFrame(maf, index=names_cr, columns=['maf'])
    fail_maf = maf < min_maf

    failed_snp_names.extend(list(names_cr[fail_maf]))
    keep_maf = ~fail_maf

    if not np.any(keep_maf):
        return names_cr[keep_maf], failed_snp_names, call_rate_df, maf_df, None

    names_maf = names_cr[keep_maf]

    # -------------------------
    # HWE exact test (much faster loop)
    # -------------------------
    n0_h = n0[keep_maf]
    n1_h = n1[keep_maf]
    n2_h = n2[keep_maf]

    hwe_p = np.array([hwe_exact_p(n0_h[i], n1_h[i], n2_h[i])
                      for i in range(len(n0_h))])

    hwe_df = pd.DataFrame(hwe_p, index=names_maf, columns=['hwe_p'])
    fail_hwe = hwe_p < min_hwe_P
    failed_snp_names.extend(list(names_maf[fail_hwe]))

    kept_final = names_maf[~fail_hwe]

    return kept_final, failed_snp_names, call_rate_df, maf_df, hwe_df

def do_snp_qc_org(snp_df, min_call_rate, min_maf, min_hwe_P):
    
    #Determine call rate.
    call_rate = 1-snp_df.isnull().sum()/len(snp_df.index)
    #Return call_rate
    selection = call_rate < min_call_rate
    #print("Pass CR: ");print(snp_df.isnull().sum())
    #print(call_rate)
    #print(call_rate < min_call_rate)
    call_rate = pd.DataFrame(data=call_rate.values, index=snp_df.columns, columns=['call_rate'])
    failed_snp_names  = list(snp_df.columns[selection])
    snp_df_c = snp_df.loc[:,list(snp_df.columns[~selection])]
    #print("Pass CR: ");print(snp_df_c.shape)
    if(len(snp_df_c.columns)==0):
        return snp_df_c.columns, failed_snp_names, call_rate, None, None
    #Determine MAF.
    genotypeCounter = np.zeros((len(snp_df_c.columns),3), dtype=int)
    for allele_index in [0,1,2]:
        genotypeCounter[:,allele_index] = np.nansum(np.around(snp_df_c.values)==allele_index,0)

    #Here we make sure that the major allele is temporarly 'coded' as 0 & directly calculate the MAF (based on allele counts and non NA samples)
    
    mac = np.zeros((len(snp_df_c.columns)), dtype=int)
    gc = np.zeros((len(snp_df_c.columns)), dtype=int)
    maf = np.zeros((len(snp_df_c.columns)), dtype=float)

    for snp in range(0, len(snp_df_c.columns)):
        if genotypeCounter[snp,0]<genotypeCounter[snp,2]:
            tmp = genotypeCounter[snp,0]
            genotypeCounter[snp,0] = genotypeCounter[snp,2]
            genotypeCounter[snp,2] = tmp
        mac[snp] = int((genotypeCounter[snp,2]*2)+genotypeCounter[snp,1])
        gc[snp] = int(genotypeCounter[snp,2]+genotypeCounter[snp,1]+genotypeCounter[snp,0])
        maf[snp] = mac[snp] / (float(2*gc[snp]))
    #Return maf
    selection = maf < min_maf
    maf = pd.DataFrame(data=maf, index=snp_df_c.columns, columns=['maf'])
    failed_snp_names.extend(list(snp_df_c.columns[selection]))
    snp_df_c = snp_df_c.iloc[:,~selection]
    #print("Pass MAF: ");print(snp_df_c.shape)
    if(len(snp_df_c.columns)==0):
        return snp_df_c.columns, failed_snp_names, call_rate, maf, None
    
    mac=mac[~selection]
    gc=gc[~selection]
    genotypeCounter = genotypeCounter[~selection,]
    
    #Determine HWE.
    hweP = np.zeros((len(snp_df_c.columns)), dtype=float)
    #print(len(snp_df_c.columns))
    #print(len(mac))
    #This can also be multi-threaded if we place it in an separate function. (And multi-threading is as easy as it seems)
    for snp in range(0, len(snp_df_c.columns)):
        rare_copies = mac[snp]
        genotypes = gc[snp]
        
        
        het_probs = np.zeros((rare_copies+1))
        # start at midpoint #

        mid = int(math.floor(rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes)))

        # check to ensure that midpoint and rare alleles have same parity #
        if int(mid % 2) != int(rare_copies % 2) :
            mid+=1
        
        curr_homr = int(math.floor((rare_copies - mid) / 2))
        curr_homc = genotypes - mid - curr_homr

        #print(het_probs.shape)
        het_probs[mid] = 1.0
        sum_values = het_probs[mid]
        for curr_hets in range(mid, 0, -2) :
            het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0))
            sum_values += het_probs[curr_hets - 2]
            ## 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote ##
            curr_homr+=1
            curr_homc+=1

        curr_homr = int(math.floor((rare_copies - mid) / 2))
        curr_homc = genotypes - mid - curr_homr

        for curr_hets in range(mid, (rare_copies), 2) :
            het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
            sum_values += het_probs[curr_hets + 2]
            curr_homr-=1
            curr_homc-=1
        
        p_hwe = 0.0
        het_probs[int(genotypeCounter[snp,1])]/= sum_values
        for index in range(0,rare_copies+1) :
            if index != int(genotypeCounter[snp,1]) :
                het_probs[index] /= sum_values
            if het_probs[index] <= het_probs[int(genotypeCounter[snp,1])] :
                p_hwe += het_probs[index];

        hweP[snp] = 1 if p_hwe > 1.0 else p_hwe;
    #Return hweP
    selection = hweP < min_hwe_P
    hweP = pd.DataFrame(data=hweP, index=snp_df_c.columns, columns=['hwe_p'])
    failed_snp_names.extend(list(snp_df_c.columns[selection]))
    snp_df_c = snp_df_c.loc[:,list(snp_df_c.columns[~selection])]
    #print("Pass HWE: ");print(snp_df_c.shape)
    return snp_df_c.columns, failed_snp_names, call_rate, maf, hweP

def do_snp_qc_stringent(snp_df, min_call_rate, min_maf, min_hwe_P, min_hmachR2):
    arr = snp_df.to_numpy(dtype=np.float64, copy=False)
    n_samples, n_snps = arr.shape
    names = snp_df.columns

    # ------------------------------------
    # CALL RATE
    # ------------------------------------
    valid = np.sum(~np.isnan(arr), axis=0)
    call_rate = valid / float(n_samples)
    fail_cr = call_rate < min_call_rate

    failed_snp_names = list(names[fail_cr])
    keep_cr = ~fail_cr
    call_rate_df = pd.DataFrame(call_rate, index=names, columns=['call_rate'])

    if not np.any(keep_cr):
        return names[keep_cr], failed_snp_names, call_rate_df, None, None, None

    arr_cr = arr[:, keep_cr]
    names_cr = names[keep_cr]

    # ------------------------------------
    # GENOTYPE COUNTS (vectorized)
    # ------------------------------------
    G = np.rint(arr_cr)
    n0 = np.sum(G == 0, axis=0)
    n1 = np.sum(G == 1, axis=0)
    n2 = np.sum(G == 2, axis=0)

    # ------------------------------------
    # MAF (vectorized)
    # ------------------------------------
    gc = n0 + n1 + n2
    minor_copies = np.minimum(2*n2 + n1, 2*n0 + n1)

    with np.errstate(divide='ignore', invalid='ignore'):
        maf = minor_copies / (2.0 * gc)

    maf_df = pd.DataFrame(maf, index=names_cr, columns=['maf'])
    fail_maf = maf < min_maf
    failed_snp_names.extend(list(names_cr[fail_maf]))
    keep_maf = ~fail_maf

    if not np.any(keep_maf):
        return names_cr[keep_maf], failed_snp_names, call_rate_df, maf_df, None, None

    names_maf = names_cr[keep_maf]

    # ------------------------------------
    # HWE P-values (optimized loop)
    # ------------------------------------
    n0_h = n0[keep_maf]
    n1_h = n1[keep_maf]
    n2_h = n2[keep_maf]

    hwe_p = np.array([
        hwe_exact_p(n0_h[i], n1_h[i], n2_h[i])
        for i in range(len(n0_h))
    ])

    hwe_df = pd.DataFrame(hwe_p, index=names_maf, columns=['hwe_p'])
    fail_hwe = hwe_p < min_hwe_P
    failed_snp_names.extend(list(names_maf[fail_hwe]))
    keep_hwe = ~fail_hwe

    if not np.any(keep_hwe):
        return names_maf[keep_hwe], failed_snp_names, call_rate_df, maf_df, hwe_df, None

    names_hwe = names_maf[keep_hwe]

    # ------------------------------------
    # MaCH R² (vectorized)
    # ------------------------------------
    n0_r2 = n0_h[keep_hwe]
    n1_r2 = n1_h[keep_hwe]
    n2_r2 = n2_h[keep_hwe]
    gc_r2 = n0_r2 + n1_r2 + n2_r2

    # Dosage mean and variance
    dosage_sum = n1_r2 * 1 + n2_r2 * 2
    dosage_sqr = n1_r2 * 1 + n2_r2 * 4

    with np.errstate(divide='ignore', invalid='ignore'):
        p = dosage_sum / (2 * gc_r2)  # allele frequency
        var_obs = (dosage_sqr / gc_r2) - (dosage_sum / gc_r2)**2
        var_exp = 2 * p * (1 - p)
        mach_r2 = var_obs / var_exp

    # Clean up invalid or >1
    mach_r2[(p <= 0) | (p >= 1)] = 1.0
    mach_r2 = np.clip(mach_r2, 0, 1)

    machr2_df = pd.DataFrame(mach_r2, index=names_hwe, columns=['mach_r2'])

    fail_r2 = mach_r2 < min_hmachR2
    failed_snp_names.extend(list(names_hwe[fail_r2]))

    kept_final = names_hwe[~fail_r2]

    return kept_final, failed_snp_names, call_rate_df, maf_df, hwe_df, machr2_df

def do_snp_qc_stringent_org(snp_df, min_call_rate, min_maf, min_hwe_P, min_hmachR2):
   
    #Determine call rate.
    call_rate = 1-snp_df.isnull().sum()/len(snp_df.index)
    selection = call_rate < min_call_rate
    #print(call_rate)
    #print(call_rate < min_call_rate)
    failed_snp_names  = list(snp_df.columns[selection])
    snp_df_c = snp_df.loc[:,list(snp_df.columns[~selection])]
    if(len(snp_df_c.columns)==0):
        return snp_df_c.columns, failed_snp_names
    #Determine MAF.
    genotypeCounter = np.zeros((len(snp_df_c.columns),3), dtype=int)
    for allele_index in [0,1,2]:
        genotypeCounter[:,allele_index] = np.nansum(np.around(snp_df_c.values)==allele_index,0)

    #Here we make sure that the major allele is temporarly 'coded' as 0 & directly calculate the MAF (based on allele counts and non NA samples)
    mac = np.zeros((len(snp_df_c.columns)), dtype=int)
    gc = np.zeros((len(snp_df_c.columns)), dtype=int)
    maf = np.zeros((len(snp_df_c.columns)), dtype=float)
    for snp in range(0, len(snp_df_c.columns)):
        if genotypeCounter[snp,0]<genotypeCounter[snp,2]:
            tmp = genotypeCounter[snp,0]
            genotypeCounter[snp,0] = genotypeCounter[snp,2]
            genotypeCounter[snp,2] = tmp
        mac[snp] = int((genotypeCounter[snp,2]*2)+genotypeCounter[snp,1])
        gc[snp] = int(genotypeCounter[snp,2]+genotypeCounter[snp,1]+genotypeCounter[snp,0])
        maf[snp] = mac[snp] / (float(2*gc[snp]))
    
    selection = maf < min_maf
    
    failed_snp_names.extend(list(snp_df_c.columns[selection]))
    snp_df_c = snp_df_c.loc[:,list(snp_df_c.columns[~selection])]
    
    if(len(snp_df_c.columns)==0):
        return snp_df_c.columns, failed_snp_names
    
    mac=mac[~selection]
    gc=gc[~selection]
    genotypeCounter = genotypeCounter[~selection,]
    
    #Determine HWE.
    hweP = np.zeros((len(snp_df_c.columns)), dtype=float)
    
    for snp in range(0, len(snp_df_c.columns)):
        rare_copies = mac[snp]
        genotypes = gc[snp]
        
        
        het_probs = np.zeros((rare_copies+1))
        # start at midpoint #

        mid = int(math.floor(rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes)))

        # check to ensure that midpoint and rare alleles have same parity #
        if int(mid % 2) != int(rare_copies % 2) :
            mid+=1
        
        curr_homr = int(math.floor((rare_copies - mid) / 2))
        curr_homc = genotypes - mid - curr_homr

        #print(het_probs.shape)
        het_probs[mid] = 1.0
        sum_values = het_probs[mid]
        for curr_hets in range(mid, 0, -2) :
            het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0))
            sum_values += het_probs[curr_hets - 2]
            ## 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote ##
            curr_homr+=1
            curr_homc+=1

        curr_homr = int(math.floor((rare_copies - mid) / 2))
        curr_homc = genotypes - mid - curr_homr

        for curr_hets in range(mid, (rare_copies), 2) :
            het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
            sum_values += het_probs[curr_hets + 2]
            curr_homr-=1
            curr_homc-=1
        
        p_hwe = 0.0
        het_probs[int(genotypeCounter[snp,1])]/= sum_values
        for index in range(0,rare_copies+1) :
            if index != int(genotypeCounter[snp,1]) :
                het_probs[index] /= sum_values
            if het_probs[index] <= het_probs[int(genotypeCounter[snp,1])] :
                p_hwe += het_probs[index];

        hweP[snp] = 1 if p_hwe > 1.0 else p_hwe;
    selection = hweP < min_hwe_P
    failed_snp_names.extend(list(snp_df_c.columns[selection]))
    snp_df_c = snp_df_c.loc[:,list(snp_df_c.columns[~selection])]
    
    if(len(snp_df_c.columns)==0):
        return snp_df_c.columns, failed_snp_names
    
    gc = gc[~selection]
    genotypeCounter = genotypeCounter[~selection,]
    
    machR2 = np.zeros((len(snp_df_c.columns)), dtype=float)
    for snp in range(0, len(snp_df_c.columns)):
        dosageSum = 0.0
        dosageSqrSum = 0.0
        dosageSum += genotypeCounter[snp,1]
        dosageSum += genotypeCounter[snp,2]*2
        dosageSqrSum += math.pow(1, 2)*genotypeCounter[snp,1]
        dosageSqrSum += math.pow(2, 2)*genotypeCounter[snp,2]
        nonMissingCount = gc[snp]
        estimatedAlleleFrequency = dosageSum / (2 * nonMissingCount)
        if (estimatedAlleleFrequency <= 0 or estimatedAlleleFrequency >= 1) :
            machR2[snp] = 1
        tmpR2 = ((dosageSqrSum / nonMissingCount) - math.pow((dosageSum / nonMissingCount), 2)) / (2 * estimatedAlleleFrequency * (1 - estimatedAlleleFrequency))
        machR2[snp] = 1 if tmpR2 > 1.0 else tmpR2
    selection = machR2 < min_hmachR2
    
    failed_snp_names.extend(list(snp_df_c.columns[selection]))
    snp_df_c = snp_df_c.loc[:,list(snp_df_c.columns[~selection])]
    
    return snp_df_c.columns, failed_snp_names, call_rate, maf, hweP, machR2

