#!/usr/bin/env python3

#----------------------------------------------------------------------
# Estimates Ne 
#
# Inputs:
#   - Z.by.year.{prefix}.tsv - Z values per year
# Output:
#   - Z.by.year.{prefix}.tsv - Z values per year
#   - dz.max_interval.{prefix}.tsv - diff in Z values in consecutive year intervals 
#----------------------------------------------------------------------

from math import sqrt, asin
from utils import load_paired_samples
import numpy as np

work_dir = "../../results/06_SNPs_stats"

# Input files
# Z.by.year.{prefix}.tsv - Z values per year (loaded dynamically)

# Output files
# Ne results file

minMAF = 0.05
z_low = 2.0*asin(sqrt(minMAF))
z_high= 2.0*asin(sqrt(1.0-minMAF))

def get_samples_years(prefixes):    
    prefixes_ne = {}
    for pre in prefixes:
        with open(f"{work_dir}/z_by_year.{pre}.tsv") as pre_fh:
            header = pre_fh.readline().strip().split("\t")
            ncols = len(header)            
            if ncols > 5:                
                start_year = header[4].split("_")[-1]
                last_year = header[-1].strip().split("_")[-1]
                prefixes_ne[pre] = [start_year, last_year]    
    return(prefixes_ne)

def compute_dz_variance(prefix):    
    dzlist = []
    mean_evar =0
    mean_dz2 = 0    
    z_by_interval_in = f"{work_dir}/z_by_year.{prefix}.tsv"
    with open(z_by_interval_in, 'r') as z_by_interval_fh:
        next(z_by_interval_fh)
        n_row = 0
        for row in z_by_interval_fh:
            cols = row.strip().split("\t")
            n_row += 1
            # Get first and last year z values
            z_first =cols[4].split(",")[0]
            z_last = cols[len(cols)-1].split(",")[0]

            if z_first != "NA"  and z_last != "NA":
                z_first = float(z_first)
                z_last = float(z_last)
                z_mean = (z_first + z_last)/2.0
                if z_mean > z_low and z_mean < z_high:
                    dz = z_first - z_last                    
                    if n_row%2:
                        dz = -dz
                    # Get first and last year SE values and estimate variance
                    se_first = float(cols[4].split(",")[1])
                    se_last = float(cols[-1].split(",")[1])
                    evar = se_first**2 + se_last **2
                    mean_evar += evar
                    mean_dz2 += (dz**2)
                    dzlist.append([dz, evar])
    return(dzlist, mean_dz2, mean_evar)

def compute_quantiles(dzlist):
    dz_values = np.array([dz for dz, _ in dzlist])
    q25 = np.percentile(dz_values, 25)
    q75 = np.percentile(dz_values, 75)
    IQR = q75 - q25
    return(IQR)

def main():    
    paired_samples = load_paired_samples()
    prefixes = sorted(list({ "_".join(k.split("_")[:-1]) for k in paired_samples.keys() }))
    samples = get_samples_years(prefixes)    
    for key, value in samples.items():
        prefix = key 
        dzlist, mean_dz2, mean_evar = compute_dz_variance(prefix)

        nx = len(dzlist)
        if nx > 0:                        
            IQR = compute_quantiles(dzlist)
            Q25 = -IQR / 2
            start_year = value[0]
            last_year = value[1]
            t = int(value[1]) - int(value[0])
            print(prefix, start_year, last_year)
            print(f"IQR={IQR:.4f}, Q25={Q25:.4f}")            
            print(f"SNPs included: {nx}")
            mean_evar = mean_evar/float(nx) 
            mean_dz2 = mean_dz2/float(nx)
            print(f"cases: {nx}, mean dz2={mean_dz2:.4f}, Estimation var={mean_evar:.4f}")
            
            # Estimate Ne
            ne_iqr = t/(2* ((IQR/1.34896)**2.0 - mean_evar))
            ne_dz2 = t/(2* (mean_dz2 - mean_evar))
            print(f"IQR Ne={ne_iqr:.2f}")
            print(f"Ez2 Ne={ne_dz2:.2f}")            
        else:
            print(f"No SNPs passed filters for {prefix}")
        
        print("-"*30)

if __name__ == '__main__':
    main()