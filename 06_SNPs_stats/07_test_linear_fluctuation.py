#!/usr/bin/env python3

#----------------------------------------------------------------------
# 
#
# Inputs:
#   - 
# Output:
#   - 
#----------------------------------------------------------------------

from pathlib import Path
import numpy as np
from math import sqrt, asin
#from scipy.stats import multivariate_normal, chi2
from utils import load_paired_samples 


work_dir = "../../results/06_SNPs_stats"

# Input files
ne_in = f"{work_dir}/Ne_estimates.tsv"
#z_by_year_in = f"{work_dir}/Z.by.year.{prefix}.tsv"   # Z values per year

# Output files


minMAF = 0.05
z_low = 2.0 * asin(sqrt(minMAF))
z_high = 2.0 * asin(sqrt(1.0 - minMAF))


def get_years_and_ne(prefix, ne_in):
    with open(ne_in) as ne_fh:
        next(ne_fh)
        for line in ne_fh:
            cols = line.strip().split("\t")
            if cols[0] == prefix:
                n_years = int(cols[3])
                ne      = float(cols[9])
                return (n_years, ne)
    return (None, None)


def tests(z_by_year_in, n_years, var_drift):
    
#    with open(ne_tests_out, "w") as ne_tests_fh:        
#       ne_tests_fh.write("\t".join([
#           "CHROM", "POS", "REF", "ALT", "tests","mu_null", "LL0", "mu98", "mu21", "LL1", 
#           "LRT_1df", "pval_1df", "LL_sat", "LRT_df", "pval_df" ]) + "\n")
        with open(z_by_year_in) as z_by_year_fh:
            cols = z_by_year_fh.readline().strip().split("\t")
            if len(cols) > 5:
                header = cols[4:]  # All years columns
                header_years = [int(h.split("_")[-1]) for h in header]
                
                for row in z_by_year_fh:
                    cols = row.rstrip().split("\t")
                    z = []
                    se = []
                    years = []
                    
                    for idx, value in enumerate(cols[4:]):
                        if not "NA" in value.split(",")[0]:
                            z.append(float(value.split(",")[0]))
                            se.append(float(value.split(",")[1]))
                            years.append(header_years[idx])                  

                    if len(z) >= n_years:
                        z_mean = sum(z)/len(z)
                        if z_mean > z_low and z_mean < z_high:
                            # Drift model
                            cov_matrix = []
                            for j in range(len(z)):
                                row_cmat = [0.0] * len(z)
                                row_cmat[j] = se[j] ** 2
                                cov_matrix.append(row_cmat)                            
                            
                            for x in range(1, len(z)):
                                yr_x = years[x] - years[0]
                                for y in range(x, len(z)):
                                    yr_y = years[y] - years[0]
                                    if y==x:
                                        cov_matrix[x][y]+=yr_y * var_drift
                                    else:
                                        cov_matrix[x][y]=float(yr_x) * var_drift
                                        cov_matrix[y][x]=cov_matrix[x][y]
                            
                            for row in cov_matrix:
                                print("\t".join(f"{v:.6e}" for v in row))
                            print("-" * 40)

def main():
    paired_samples = load_paired_samples()    
    prefixes = sorted(list({ "_".join(k.split("_")[:-1]) for k in paired_samples.keys()}))    

    # Process for each prefix
    for pre in prefixes:
        n_years, ne = get_years_and_ne(pre, ne_in)
        z_by_year_in = f"{work_dir}/z_by_year.{pre}.tsv"
       # ne_tests_out = f"{work_dir}/Ne_tests_{pre}.tsv"
        if ne is not None:
            var_drift = 1.0 / (2.0 * ne) if ne > 0 else 1.0 / (2.0 * 1.0)
            z_path=Path(z_by_year_in)
            if z_path.is_file():
                print(pre)       
                print(n_years)
                tests(z_by_year_in, n_years, var_drift)

if __name__ == "__main__":
    main()