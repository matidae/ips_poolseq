#!/usr/bin/env python3

#----------------------------------------------------------------------
# Recalculate null variance for remaining SNPs after filtering results from previous steps
#
# Inputs:
#   - null_variance_summary.recalc.tsv: recalculated null variance estimates per sample
#   - genic_m_and_z.filter.tsv: filtered allele frequencies and counts per SNP
# Output:
#   - Z.by.year.tsv — Z values per year
#   - DZ.by.interval.tsv — diff in Z values in consecutive year intervals 
#----------------------------------------------------------------------

import re
from math import sqrt, asin
from utils import load_paired_samples, load_null_variance_recalc


work_dir = "../../results/06_SNPs_stats"

# Input files
null_var_in = f"{work_dir}/null_variance_summary.recalc.tsv"  # Null variance data input
m_and_z_in = f"{work_dir}/genic_m_and_z.filter.tsv"           # Genic m and z data input

# Output files
#z_by_year_out = f"{work_dir}/Z.by.year.tsv"            # Z values per year
#dz_by_year_out = f"{work_dir}/DZ.by.interval.tsv"      # Diff in Z values in consecutive year intervals 

def calculate_dz_by_year(prefixes, paired_samples, null_var, m_and_z_in):
    minMAF = 0.05
    zlow = 2.0*asin(sqrt(minMAF))
    zhigh= 2.0*asin(sqrt(1.0-minMAF))
    min_depth = 100

    for pre in prefixes:
        paired_samples_same= {}
        # Get all samples from same country and season
        for key, val in paired_samples.items():
            if key.startswith(pre):                
                paired_samples_same[key] = val
        
        if len(paired_samples_same) > 1:
            vstat = {}
            for k in paired_samples_same.keys():
                vstat[k]=[0,0.0,0.0]

            # Output files for this prefix
            z_by_year_out = f"{work_dir}/Z.by.year.{pre}.tsv"
            dz_by_year_out = f"{work_dir}/DZ.by.interval.{pre}.tsv"
                       
            with open(m_and_z_in, "r") as m_and_z_fh, \
                open(z_by_year_out, "w") as z_by_year_fh, \
                open(dz_by_year_out, "w") as dz_by_year_fh :

                # Header
                header = m_and_z_fh.readline().strip().split("\t")
                z_by_year_fh.write("CHROM\tPOS\tREF\tALT")
                dz_by_year_fh.write("CHROM\tPOS\tREF\tALT")

                for sample in paired_samples_same:
                    z_by_year_fh.write(f"\t{sample}")
                for i in range(len(paired_samples_same) - 1):
                    dz_by_year_fh.write(f"\t{pre}_interval{i+1}")
                z_by_year_fh.write("\n")
                dz_by_year_fh.write("\n")

                for line in m_and_z_fh:
                    cols = line.strip().split('\t')
                    row_start = "\t".join(cols[:4])
                    z_by_year_fh.write(row_start)
                    dz_by_year_fh.write(row_start)
                    last = None
                    m_total = 0

                    # Iterate over each group of replicates
                    for sample in paired_samples_same:
                        rep_a = cols[paired_samples[sample][0]].split(",")
                        m_a = int(rep_a[0])
                        rep_b = cols[paired_samples[sample][1]].split(",")
                        m_b = int(rep_b[0])
                        # Get total replicates depth
                        m_total += (m_a + m_b)

                   # next(m_and_z_fh)

                        if m_a >= min_depth and m_b >= min_depth:
                            z_a = float(rep_a[1])                    
                            z_b = float(rep_b[1])
                            z_mean=(z_a + z_b)/2.0
                            var = (null_var[key] + 1.0/float(m_a) + 1.0/float(m_b))/4.0
                            std_error = sqrt(var)

                            if z_mean > zlow and z_mean < zhigh:
                                z_by_year_fh.write('\t'+str(z_mean)+','+str(std_error)+','+str(m_a)+','+str(m_b))
                            else:
                                z_by_year_fh.write("\tNA,NA,0,0")

                            if (last is not None and z_mean+last[0])/2.0 >zlow and (z_mean+last[0])/2.0 <zhigh: 
                                dz2 = z_mean - last[0]
                                evar = (null_var[key] + 1.0/float(m_a) + 1.0/float(m_b))/4.0 + last[1]
                                z_by_year_fh.write('\t', dz2, evar)
                            else:
                                z_by_year_fh.write('\tNA')

def main():
    null_var = load_null_variance_recalc()
    paired_samples = load_paired_samples()
    prefixes = sorted(list({ "_".join(k.split("_")[:-1]) for k in paired_samples.keys() }))
    calculate_dz_by_year(prefixes, paired_samples, null_var, m_and_z_in)

if __name__ == '__main__':
    main()
