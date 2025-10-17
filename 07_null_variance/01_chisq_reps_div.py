#!/usr/bin/env python3

#----------------------------------------------------------------------
# Test SNP reliability based on replicate divergence using allele frequencies as (z-scores).
# Compute chi-square statistic and p-value based on null variance (p-values < 0.01 are unreliable). 
#
# Inputs:
#   - null_variance_summary.tsv: null variance estimates per sample
#   - samples_paired.tsv: sample names and replicate column indices
#   - genic_m_and_z.tsv: allele frequencies and counts per SNP
#
# Output:
#   - snpdev.m_and_z.tsv: SNP-level coverage, comparisons, chi-square stat, and p-value
#----------------------------------------------------------------------
import sys
sys.path.append("../utils") 
from math import sqrt, pi
from scipy.stats import chi2
from utils import load_paired_samples
import pandas as pd 

work_dir = "../../results/07_null_variance"
input_dir = "../../results/06_SNPs_stats"

# Input files
m_and_z_in = f"{input_dir}/genic_m_and_z.tsv"
null_var_in = f"{work_dir}/null_variance_summary.tsv"
# Ouput file
snpdev_out = f"{work_dir}/snpdev_m_and_z.tsv"
cov_bin_out = f"{work_dir}/snp_depth_div_bin.tsv"

# Parameters
#mincalls = 20  # Minimum number of SNPs per bulk
min_depth = 50 # Lowest count to include a sample

def load_null_variance(null_var_in):    
    df = pd.read_csv(null_var_in, sep='\t')
    null_var = dict(zip(df.iloc[:, 0], df.iloc[:, 5]))    
    return null_var

def aggregate_SNPs_by_coverage_bin(coverage_bin_stats, dcat, p_value):
    if dcat not in coverage_bin_stats:
        coverage_bin_stats[dcat] = [0, 0]
    coverage_bin_stats[dcat][0] += 1
    if p_value < 0.01:
        coverage_bin_stats[dcat][1] += 1    

def write_SNP_stats(snpdev_out, cols_mz, m_total, snpstats, p_value, is_valid):
    if is_valid:        
        snpdev_out.write(f"{cols_mz[0]}\t{cols_mz[1]}\t{m_total}\t{snpstats[0]}\t{snpstats[1]:.4f}\t{p_value:.6f}\n")
    else:
        snpdev_out.write(f"{cols_mz[0]}\t{cols_mz[1]}\t{m_total}\t{snpstats[0]}\t-99\t1.0\n")       

def main():
    null_var = load_null_variance(null_var_in)
    paired_samples = load_paired_samples()

    with open(m_and_z_in, "r") as m_and_z_fh, open(snpdev_out, "w") as snpdev_fh:
        snpdev_fh.write("CHROM\tPOS\tDepth\tNcomp\tChi2\tpval\n")
        next(m_and_z_fh)    
        coverage_bin_stats={}    

        for line in m_and_z_fh:
            cols = line.strip().split('\t')       
            snpstats=[0, 0.0]        
            m_total = 0

            # Iterate over each group of replicates
            for sample in paired_samples:
                rep_a = cols[paired_samples[sample][0]].split(",")
                m_a = int(rep_a[0])
                rep_b = cols[paired_samples[sample][1]].split(",")
                m_b = int(rep_b[0])
                # Get total replicates depth
                m_total += (m_a + m_b)

                #If both replicate pass minimum read depth threshold, then calculate replicate variance
                if m_a >= min_depth and m_b >= min_depth:
                    z_a = float(rep_a[1])
                    z_b = float(rep_b[1])
                    if (z_a + z_b) > 0 and (z_a + z_b) < 2 * pi: #sum of z values should be between 0 and 2 Pi
                        #Standardized allele freq. difference, adjusted by null variance and depth.
                        rx = (z_a - z_b)/ sqrt(null_var[sample] + 1/m_a + 1/m_b) # E[rx] = nullvar

                        snpstats[0]+=1 # Count valid comparisons
                        snpstats[1]+= rx*rx # Accumulate chi-squared-like statistic

            dcat=int(m_total/500)  # Bin SNP by total read coverage (e.g., 0,1,... for every 500 reads)
            p_value = 1.0 
            if snpstats[0]>0:  # If there were valid comparisons
                if snpstats[1]>0:
                    p_value=1.0 - chi2.cdf(snpstats[1], snpstats[0])  # Chi-squared p-value
                else:
                    p_value=1.0  # No divergence            
                # Write SNP-level statistics
                write_SNP_stats(snpdev_fh, cols, m_total, snpstats, p_value, is_valid=True)
                # Aggregate statistics by coverage bin
                aggregate_SNPs_by_coverage_bin(coverage_bin_stats, dcat, p_value)             
            else:  # No valid comparisons
                write_SNP_stats(snpdev_fh, cols, m_total, snpstats, p_value, is_valid=False)

    with open(cov_bin_out, "w") as cov_bin_fh:
        cov_bin_fh.write("depth_start\tSNP_count\tprop_pval_lt_0.01\n")
        for dcat in sorted(coverage_bin_stats):
            coverage_start = dcat * 500
            num_snps = coverage_bin_stats[dcat][0]
            prop = coverage_bin_stats[dcat][1] / num_snps
            cov_bin_fh.write(f"{coverage_start}\t{num_snps}\t{prop:.4f}\n")

if __name__ == "__main__":
    main()