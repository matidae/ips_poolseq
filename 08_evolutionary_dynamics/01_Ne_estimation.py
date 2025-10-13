#!/usr/bin/env python3

#----------------------------------------------------------------------
# Estimates effective population size (Ne). Reads the z-transformed AF across years
# Calculates the changes in AF as dz between first and last year of each population
# It corrects for sampling error and estimates Ne using both IQR and mean-square 
#
# Inputs:
#   - z_year.{prefix}.tsv - z values per year per population
# Output:
#   - Ne_estimates.tsv - Ne estimates for each population
#----------------------------------------------------------------------

from math import sqrt, asin
from utils import load_paired_samples
import numpy as np

work_dir = "../../results/08_evolutionary_dynamics"

# Input files
# z_year.{prefix}.tsv - z values per year (loaded dynamically)

# Output files
ne_out = f"{work_dir}/Ne_estimates.tsv" # Ne results file

minMAF = 0.05
z_low = 2.0*asin(sqrt(minMAF))
z_high= 2.0*asin(sqrt(1.0-minMAF))

def get_samples_years(prefixes):    
    prefixes_ne = {}
    for pre in prefixes:
        with open(f"{work_dir}/z_year.{pre}.tsv") as pre_fh:
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
    z_interval_in = f"{work_dir}/z_year.{prefix}.tsv"
    dz_out = f"{work_dir}/dz_max.{prefix}.tsv"
    with open(z_interval_in, 'r') as z_interval_fh, open(dz_out, "w") as dz_fh:
        next(z_interval_fh)
        n_row = 0
        for row in z_interval_fh:
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
                    # Save dz values of each population
                    dz_fh.write(f"{dz:.5f}\n")
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
    iqr = q75 - q25
    return(iqr)

def main():    
    paired_samples = load_paired_samples()
    prefixes = sorted(list({ "_".join(k.split("_")[:-1]) for k in paired_samples.keys() }))
    samples = get_samples_years(prefixes)
    with open(ne_out, "w") as ne_fh:
         # Write header 
        ne_fh.write("\t".join(["prefix", "start", "end", "t", "SNPs",
            "mean_dz2", "mean_evar", "IQR", "Ne_IQR", "Ne_dz2"]) + "\n")
        for key, value in samples.items():
            prefix = key 
            dzlist, mean_dz2, mean_evar = compute_dz_variance(prefix)

            nx = len(dzlist)
            if nx > 0:                        
                iqr = compute_quantiles(dzlist)
                q25 = -iqr / 2
                start_year = value[0]
                last_year = value[1]
                t = int(value[1]) - int(value[0])                
                mean_evar = mean_evar/float(nx) 
                mean_dz2 = mean_dz2/float(nx)                
                # Estimate Ne
                ne_iqr = t/(2* ((iqr/1.34896)**2.0 - mean_evar))
                ne_dz2 = t/(2* (mean_dz2 - mean_evar))                
                # Save Ne to table                
                ne_fh.write("\t".join(map(str,[prefix, start_year, last_year, t, nx,
                f"{mean_dz2:.4f}", f"{mean_evar:.4f}", f"{iqr:.4f}", int(ne_iqr), int(ne_dz2)])) + "\n")            

if __name__ == '__main__':
    main()