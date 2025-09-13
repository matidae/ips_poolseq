#!/usr/bin/env python3

#----------------------------------------------------------------------
# Recalculate null variance for remaining SNPs after filtering results from previous steps
#
# Inputs:
#   - null_variance_summary.recalc.tsv: recalculated null variance estimates per sample
#   - genic_m_and_z.filter.tsv: filtered allele frequencies and counts per SNP
# Output:
#   - Z.by.year.{prefix}.tsv — Z values per year
#   - DZ.by.interval.{prefix}.tsv — diff in Z values in consecutive year intervals 
#----------------------------------------------------------------------

from math import sqrt, asin
from utils import load_paired_samples, load_null_variance_recalc


work_dir = "../../results/06_SNPs_stats"

# Input files
null_var_in = f"{work_dir}/null_variance_summary.recalc.tsv"  # Null variance data input
m_and_z_in = f"{work_dir}/genic_m_and_z.filter.2.tsv"           # Genic m and z data input

# Output files
#z_by_year_out = f"{work_dir}/Z.by.year.{prefix}.tsv"         # Z values per year
#dz_by_year_out = f"{work_dir}/DZ.by.interval.{prefix}.tsv"   # Diff in Z values in consecutive year intervals 

def calculate_z_dz(prefixes, paired_samples, null_var, m_and_z_in):
    minMAF = 0.05
    zlow = 2.0*asin(sqrt(minMAF))
    zhigh= 2.0*asin(sqrt(1.0-minMAF))
    min_depth = 50

    for pre in prefixes:
        paired_samples_set= {}
        paired_samples_year= {}
        # Get all samples from same country and season
        for key, val in paired_samples.items():
            if key.startswith(pre):                
                paired_samples_set[key] = val
                paired_samples_year[key] = key.split("_")[2]
       # if len(paired_samples_same) > 1:
            # Output files for this prefix
            z_by_year_out = f"{work_dir}/Z.by.year.{pre}.tsv"
            dz_by_year_out = f"{work_dir}/DZ.by.interval.{pre}.tsv"
                       
            with open(m_and_z_in, "r") as m_and_z_fh, \
                open(z_by_year_out, "w") as z_by_year_fh, \
                open(dz_by_year_out, "w") as dz_by_year_fh :

                # Header
                z_by_year_fh.write("CHROM\tPOS\tREF\tALT")
                dz_by_year_fh.write("CHROM\tPOS\tREF\tALT")
                for sample in paired_samples_set:
                    z_by_year_fh.write(f"\t{sample}")
                
                for i in range(len(paired_samples_set) - 1):
                    dz_by_year_fh.write(f"\t{pre}_dif_{i+1}")
                z_by_year_fh.write("\n")
                dz_by_year_fh.write("\n")

                next(m_and_z_fh)
                for line in m_and_z_fh:                    
                    cols = line.strip().split('\t')
                    row_start = "\t".join(cols[:4])
                    z_by_year_fh.write(row_start)
                    dz_by_year_fh.write(row_start)
                    last = None
                    m_total = 0

                    # Iterate over each group of replicates
                    for sample in paired_samples_set:
                        rep_a = cols[paired_samples[sample][0]].split(",")
                        m_a = int(rep_a[0])
                        rep_b = cols[paired_samples[sample][1]].split(",")
                        m_b = int(rep_b[0])
                        # Get total replicates depth
                        m_total += (m_a + m_b)

                        if m_a >= min_depth and m_b >= min_depth:
                            z_a = float(rep_a[1])                    
                            z_b = float(rep_b[1])
                            z_mean=(z_a + z_b)/2.0
                            var = (null_var[key] + 1.0/float(m_a) + 1.0/float(m_b))/4.0
                            std_error = sqrt(var)

                            if z_mean > zlow and z_mean < zhigh:
                                z_by_year_fh.write(f"\t{z_mean:.4f},{std_error:.4f},{m_a},{m_b}")                                
                            else:
                                z_by_year_fh.write("\tNA,NA,0,0")
                        else:
                            z_by_year_fh.write("\tNA,NA,0,0")

##############
                    
                        if (last is not None and z_mean+last[0])/2.0 >zlow and (z_mean+last[0])/2.0 <zhigh: 
                            dz2 = z_mean - last[0]
                            evar = (null_var[key] + 1.0/float(m_a) + 1.0/float(m_b))/4.0 + last[1]
                            dz_by_year_fh.write(f"\t{dz2:.4f},{evar:.4f}")
                        else:
                            dz_by_year_fh.write('\tNA')
                        last = (z_mean, var) if var is not None else None
                    z_by_year_fh.write("\n")
                    dz_by_year_fh.write("\n")


def main():
    null_var = load_null_variance_recalc()
    paired_samples = load_paired_samples()
    prefixes = sorted(list({ "_".join(k.split("_")[:-1]) for k in paired_samples.keys() }))    
    calculate_z_dz(prefixes, paired_samples, null_var, m_and_z_in)

if __name__ == '__main__':
    main()
