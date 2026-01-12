#!/usr/bin/env python3

#----------------------------------------------------------------------
# Calculates the mean value of z for each timepoint and filters based on depth and MAF
#
# Inputs:
#   - genic_m_and_z.filter.tsv: filtered allele frequencies and counts per SNP
# Output:
#   - z_mean.{prefix}.tsv — Z values per year_season
#----------------------------------------------------------------------

import sys
sys.path.append("./utils")
import os
from math import sqrt, asin
from utils import load_paired_samples, load_null_variance_recalc

work_dir = "../results/08_models"
input_dir = "../results/07_null_variance"

# Input files
m_and_z_in = f"{input_dir}/genic_m_and_z.filter.tsv"

# Output files
#z_mean_out = f"{work_dir}/z_mean.{prefix}.tsv"

#Filter rules
minMAF = 0.05
zlow = 2.0*asin(sqrt(minMAF))
zhigh= 2.0*asin(sqrt(1.0-minMAF))
min_depth = 100

def calculate_z_mean(prefixes, paired_samples, null_var, m_and_z_in):
    season_order = {"E": 0, "L": 1}

    for pre in prefixes:
        # Filter only the samples for this prefix
        filtered_samples = {k: v for k, v in paired_samples.items() if k.startswith(pre)}

        # Sort keys by year, then by season order
        ordered_samples = sorted(
            filtered_samples.keys(),
            key=lambda k: (int(k.split("_")[2]), season_order[k.split("_")[1]])
        )

        # Output files for this prefix
        z_mean_out = f"{work_dir}/z_mean.{pre}.tsv"
        
        with open(m_and_z_in, "r") as m_and_z_fh, \
            open(z_mean_out, "w") as z_mean_fh: 
            # Write file header 
            z_mean_fh.write("CHROM\tPOS\tREF\tALT")
            for sample in ordered_samples:
                z_mean_fh.write(f"\t{sample}")
            z_mean_fh.write("\n")

            # Start analyzing m and z file
            next(m_and_z_fh)
            for line in m_and_z_fh:
                cols = line.strip().split('\t')
                row_start = "\t".join(cols[:4])
                # Write in z_mean.{pre}.tsv: chromo, coor, ref, alt ...
                z_mean_fh.write(row_start)
                # Iterate over each group of replicates
                for sample in ordered_samples:
                    rep_a = cols[paired_samples[sample][0]].split(",")
                    m_a = int(rep_a[0])
                    rep_b = cols[paired_samples[sample][1]].split(",")
                    m_b = int(rep_b[0])
                    
                    if m_a >= min_depth and m_b >= min_depth:
                        z_a = float(rep_a[1])
                        z_b = float(rep_b[1])
                        z_mean=(z_a + z_b)/2.0
                        var = (null_var[sample] + 1.0/float(m_a) + 1.0/float(m_b))/4.0
                        std_error = sqrt(var)
                        # Write in z_mean.{pre}.tsv ... z_mean, SE, depth_rep_b, depth_rep_b
                        if z_mean > zlow and z_mean < zhigh:
                            z_mean_fh.write(f"\t{z_mean:.4f},{std_error:.4f},{m_a},{m_b}")
                        else:
                            z_mean_fh.write(f"\tNA,Z_FIL,{m_a},{m_b}")
                    else:
                        z_mean_fh.write(f"\tNA,D_FIL,{m_a},{m_b}")
                z_mean_fh.write("\n") 

def main():
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    null_var = load_null_variance_recalc()
    paired_samples = load_paired_samples()
    prefixes = sorted(list({ k.split("_")[0] for k in paired_samples.keys() }))
 
 
    calculate_z_mean(prefixes, paired_samples, null_var, m_and_z_in)

if __name__ == '__main__':
    main()
