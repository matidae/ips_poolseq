#!/usr/bin/env python3

#----------------------------------------------------------------------
# Calculates the values of z for each year and the difference (dz) between consecutive years
#
# Inputs:
#   - genic_m_and_z.filter.tsv: filtered allele frequencies and counts per SNP
# Output:
#   - Z.by.year.{prefix}.tsv — Z values per year
#   - DZ.by.interval.{prefix}.tsv — diff in Z values in consecutive year intervals 
#----------------------------------------------------------------------

from math import sqrt, asin
import os
import sys
import numpy as np
from scipy.stats import multivariate_normal, chi2
from utils import load_paired_samples 


work_dir = "../../results/06_SNPs_stats"

# Input files
ne_in = f"{work_dir}/Ne_estimates.tsv"
 
minMAF = 0.05
z_low = 2.0 * asin(sqrt(minMAF))
z_high = 2.0 * asin(sqrt(1.0 - minMAF))


Ne_assumed = 11790.0
VarDrift = 1.0 / (2.0 * Ne_assumed)

def get_years_and_ne(prefix, ne_in):
    with open(ne_in) as ne_fh:
        next(ne_fh)
        for line in ne_fh:
            cols = line.strip().split("\t")
            if cols[0] == prefix:
                n_years = int(cols[3])
                ne      = float(cols[9])
                return (n_years, ne)
    raise ValueError(f"Prefix '{prefix}' not found in {ne_in}")


def main():
    paired_samples = load_paired_samples()    
    prefixes = sorted(list({ "_".join(k.split("_")[:-1]) for k in paired_samples.keys()}))
    print(prefixes)

    # Process for each prefix
    for pre in prefixes:
        n_years, ne = get_years_and_ne(pre, ne_in)
        print(pre, n_years, ne)      

if __name__ == "__main__":
    main()