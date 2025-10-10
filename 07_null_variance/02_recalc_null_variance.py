#!/usr/bin/env python3

#----------------------------------------------------------------------
# Recalculate null variance for remaining SNPs after filtering results from previous steps
#
# Inputs:
#   - null_variance_summary.tsv: null variance estimates per sample
#   - genic_m_and_z.tsv: allele frequencies and counts per SNP
#   - snpdev_m_and_z.tsv: SNP-level coverage, comparisons, chi-square stat, and p-value
# Output:
#   - genic_m_and_z.filter.tsv — Filtered SNP data, only for SNPs passing criteria.
#   - null_variance_summary.recalc.tsv — Estimated null variance per sample group
#----------------------------------------------------------------------

from math import sqrt, asin
from utils import load_paired_samples, load_depth_threshold

work_dir = "../../results/07_null_variance"
input_dir = "../../results/06_SNPs_stats"

# Input files
m_and_z_in = f"{input_dir}/genic_m_and_z.tsv"           # Genic m and z data input
null_var_in = f"{work_dir}/null_variance_summary.tsv"  # Null variance data input
snpdev_in = f"{work_dir}/snpdev_m_and_z.tsv"           # SNP deviation input file

# Output files
snp_filter_out = f"{work_dir}/genic_m_and_z.filter.tsv"      # Filtered SNP data output
null_var_recalc_out = f"{work_dir}/null_variance_summary.recalc.tsv"  # Recalculated null variance output

#Variables
min_depth_rep = 100 # minimum depth in each replicate
minMAF = 0.05
zlow = 2.0*asin(sqrt(minMAF))
zhigh= 2.0*asin(sqrt(1.0-minMAF))

#dev2 = 3.0 #filter for genotype data

# Filter SNPs by p-value from Chi-square test. Discard those with p-value<0.01
def filter_SNPs(snpdev_in, min_depth, max_depth):
    valid_SNPs = set()
    pval_min = 0.01
    with open(snpdev_in, "r") as snpdev_fh:
        next(snpdev_fh)
        for line in snpdev_fh:
            cols = line.strip().split('\t')
            pval = float(cols[5])
            depth = int(cols[2])
            if  pval >= pval_min and depth >= min_depth and depth <= max_depth: 
                #and abs(float(cols[8])) <= dev2 for genotype data
                valid_SNPs.add(f"{cols[0]}_{cols[1]}")    
    return valid_SNPs

# Process SNPs file, filter SNPs, and update stats
def process_snps(m_and_z_in, snp_filter_out, valid_snps, paired_samples, stats):    
    with open(m_and_z_in, 'r') as m_and_z_fh, open(snp_filter_out, 'w') as snp_filter_fh:
        header = next(m_and_z_fh)
        snp_filter_fh.write(header)
        for line in m_and_z_fh:
            cols = line.strip().split('\t')
            snp_key = f"{cols[0]}_{cols[1]}"
            if snp_key not in valid_snps:
                continue
            snp_filter_fh.write(line)            
            for sample, (idx_a, idx_b) in paired_samples.items():                
                m1, z1 = cols[idx_a].split(',')
                m2, z2 = cols[idx_b].split(',')
                m1, m2 = int(m1), int(m2)
                if m1 < min_depth_rep or m2 < min_depth_rep:
                    continue

                z1, z2 = float(z1), float(z2)
                z_mean = (z1 + z2) / 2.0
                if not (zlow < z_mean < zhigh):
                    continue

                dz2 = (z1 - z2) ** 2
                inv_depth = 1 / m1 + 1 / m2

                stats[sample]['count'] += 1
                stats[sample]['sum_dz2'] += dz2
                stats[sample]['sum_inv_depth'] += inv_depth

# Recalculate null variance for each sample and write to file
def calculate_null_variance(null_var_recalc_out, stats):    
    with open(null_var_recalc_out, 'w') as null_var_recalc_fh:
        null_var_recalc_fh.write("Sample\tNull_var\tn_SNP\tMean_DZ2\tMean_inv_depth\trdprop\n")
        for sample, data in stats.items():
            if data['count'] == 0:
                continue
            # Calculate means 
            mean_dz2 = data['sum_dz2'] / data['count']
            mean_inv_depth = data['sum_inv_depth'] / data['count']
            # Estimate null variance
            null_var = mean_dz2 - mean_inv_depth
            # Prop. of variance explained by read depth
            rd_prop = mean_inv_depth/mean_dz2
            null_var_recalc_fh.write(
                f"{sample}\t{null_var:.6f}\t{data['count']}\t{mean_dz2:.6f}\t"
                f"{mean_inv_depth:.6f}\t{rd_prop:.4f}\n")

def main():
    min_depth, max_depth = load_depth_threshold()
    valid_SNPs = filter_SNPs(snpdev_in, min_depth, max_depth)
    paired_samples = load_paired_samples()
    stats = {sample: {'count': 0, 'sum_dz2': 0.0, 'sum_inv_depth': 0.0} for sample in paired_samples}

    process_snps(m_and_z_in, snp_filter_out, valid_SNPs, paired_samples, stats)
    calculate_null_variance(null_var_recalc_out, stats)

if __name__ == "__main__":
    main()