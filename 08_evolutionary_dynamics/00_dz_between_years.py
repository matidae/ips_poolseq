#!/usr/bin/env python3

#----------------------------------------------------------------------
# Calculates the values of z for each year and the difference (dz) between consecutive years
#
# Inputs:
#   - genic_m_and_z.filter.tsv: filtered allele frequencies and counts per SNP
# Output:
#   - z_year.{prefix}.tsv — Z values per year
#   - dz_interval.{prefix}.tsv — diff in Z values in consecutive year intervals 
#----------------------------------------------------------------------

import sys
sys.path.append("../utils")
import os
from math import sqrt, asin
from utils import load_paired_samples, load_null_variance_recalc

work_dir = "../../results/08_evolutionary_dynamics"
input_dir = "../../results/07_null_variance"

# Input files
m_and_z_in = f"{input_dir}/genic_m_and_z.filter.tsv"     # Genic m and z data input

# Output files
#z_by_year_out = f"{work_dir}/z_year.{prefix}.tsv"       # z values per year
#dz_by_year_out = f"{work_dir}/z_interval.{prefix}.tsv"  # ziff in z values in consecutive intervals 

def calculate_z_dz(prefixes, paired_samples, null_var, m_and_z_in):
    minMAF = 0.05
    zlow = 2.0*asin(sqrt(minMAF))
    zhigh= 2.0*asin(sqrt(1.0-minMAF))
    min_depth = 80

    for pre in prefixes:
        paired_samples_set= {}
        paired_samples_year_list= []
        # Get all samples from same region and season
        for key, val in paired_samples.items():
            if key.startswith(pre):
                paired_samples_set[key] = val
                paired_samples_year_list.append(key.split("_")[2])
        paired_samples_year_list = sorted(list(set(paired_samples_year_list)))
        num_years = len(paired_samples_year_list)
        # Output files for this prefix
        z_by_year_out = f"{work_dir}/z_year.{pre}.tsv"
        dz_by_interval_out = f"{work_dir}/dz_interval.{pre}.tsv"
        
        with open(m_and_z_in, "r") as m_and_z_fh, \
            open(z_by_year_out, "w") as z_by_year_fh:            

            # Conditionally open dz file for writing
            dz_by_interval_fh = None
            if num_years > 1:
                dz_by_interval_fh = open(dz_by_interval_out, "w")

            # Write header (common fields for each sample)
            z_by_year_fh.write("CHROM\tPOS\tREF\tALT")            
            # Write header (sample-specific year field for Z values)
            for sample in paired_samples_set:
                z_by_year_fh.write(f"\t{sample}")
            z_by_year_fh.write("\n")
            
            # Write dz file header if applicable
            if dz_by_interval_fh:
                dz_by_interval_fh.write("CHROM\tPOS\tREF\tALT")
                # Write header (sample-specific year_interval field for dZ values)
                for i in range(len(paired_samples_year_list) - 1):                                    
                    pre_year = paired_samples_year_list[i]
                    pos_year = paired_samples_year_list[i+1]                
                    dz_by_interval_fh.write(f"\t{pre}_{pre_year}-{pos_year}")            
                dz_by_interval_fh.write("\n")
            # End writing headers
            # Start analyzing m and z file
            next(m_and_z_fh)
            for line in m_and_z_fh:                    
                cols = line.strip().split('\t')
                row_start = "\t".join(cols[:4])
                # Write in z_year.{pre}.tsv: chromo, coor, ref, alt ...
                z_by_year_fh.write(row_start)
                # Conditionally write if dz
                if dz_by_interval_fh:
                    dz_by_interval_fh.write(row_start)
                #Set last value of z to NA
                last_z = ["NA","NA"]                
                year_count = 0

                # Iterate over each group of replicates
                for sample in paired_samples_set:
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

                        # Write in z_year.{pre}.tsv ... z_mean, SE, depth_rep_b, depth_rep_b
                        if z_mean > zlow and z_mean < zhigh:
                            z_by_year_fh.write(f"\t{z_mean:.4f},{std_error:.4f},{m_a},{m_b}")
                        else:
                            z_by_year_fh.write("\tNA,NA,0,0")

                        # Write in dz_interval.{pre}.tsv ... z_mean, SE, depth_rep_b, depth_rep_b
                        if dz_by_interval_fh:
                            if year_count > 0:
                                if last_z[0] != 'NA':
                                    z_mean_interval = (z_mean + last_z[0])/2.0
                                    if z_mean_interval > zlow and z_mean_interval < zhigh:
                                        dz2 = z_mean - last_z[0]
                                        evar = (null_var[sample] + 1.0/float(m_a) + 1.0/float(m_b))/4.0 + last_z[1]
                                        dz_by_interval_fh.write(f"\t{dz2:.4f},{evar:.4f}")
                                    else:
                                        dz_by_interval_fh.write("\tNA,NA")
                                else:
                                    dz_by_interval_fh.write("\tNA,NA")
                            else:                            
                                year_count+=1
                            last_z = [z_mean, (null_var[sample] + 1.0/float(m_a) + 1.0/float(m_b))/4.0]
                    else:
                        z_by_year_fh.write("\tNA,NA,0,0")                     
                        if dz_by_interval_fh:
                            if year_count > 0:
                                dz_by_interval_fh.write("\tNA,NA")
                            year_count+=1
                        last_z = ["NA","NA"]

                z_by_year_fh.write("\n")
                if dz_by_interval_fh:
                    dz_by_interval_fh.write("\n")

def main():
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    null_var = load_null_variance_recalc()
    paired_samples = load_paired_samples()
    prefixes = sorted(list({ "_".join(k.split("_")[:-1]) for k in paired_samples.keys() }))
        
    calculate_z_dz(prefixes, paired_samples, null_var, m_and_z_in)

if __name__ == '__main__':
    main()
