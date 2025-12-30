#!/usr/bin/env python3

#-------------------------------------------------------------------------------
# Calculates read depth (m) and Fisher-Ford transformation (z) for allele frequencies
# Removes SNPs with depth values of outliers >1.5 median and <0.5 median
#
# Input: 
#   - genic_depth_stats.tsv: summary statistics (mean, median, quantiles, IQR, whiskers) 
#   - genic_readcounts.tsv: file with CHR, POS, REF, ALT, ADs* for all samples
# Output: 
#   - genic_m_and_z.tsv: tsv file with m and z values for each SNP where min_depth <= total_depth <= max_depth
#   - genic_m_and_z.excluded.tsv: tsv file with m and z values that where previous condition does not apply
#-------------------------------------------------------------------------------
import sys
import os
from math import sqrt, asin
sys.path.append("./utils")
from utils import parse_counts, load_depth_threshold

work_dir = "../results/06_SNPs_stats"
input_dir = "../results/04_varcalls"

# Input files
depth_stats_in = f"{input_dir}/genic_depth_stats.tsv"
readcounts_in = f"{input_dir}/genic_readcounts.tsv"
# Output files
m_and_z_out = f"{work_dir}/genic_m_and_z.tsv"
m_and_z_excluded_out = f"{work_dir}/genic_m_and_z.excluded.tsv"

min_depth, max_depth = load_depth_threshold()

# Calculates Fisher-Ford arcsine square root transformation of allele frequencies
def calculate_z_score(ref, alt):    
    m = ref + alt
    if m == 0:
        return 0, 'NA', 'NA'
    p = ref / m  # Calculates p as ref/(ref+alt)
    z = 2 * asin(sqrt(p)) # Fisher-Ford transform of p
    return m, f"{p:.4f}", f"{z:.4f}"

def main():
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    with open(readcounts_in, 'r') as readcounts_fh, open(m_and_z_out, 'w') as \
    m_and_z_fh, open(m_and_z_excluded_out, 'w') as m_and_z_excluded_fh:
        header = next(readcounts_fh)
        m_and_z_fh.write(header)
        for line in readcounts_fh:        
            cols = line.strip().split('\t')        
            # Parse sample read counts
            sample_counts = [parse_counts(x) for x in cols[4:]]            
            total_depth = sum(ref + alt for ref, alt in sample_counts)           
            # Filter by SNPs position depth sum           
            if min_depth <= total_depth <= max_depth:
                m_and_z_fh.write('\t'.join(cols[:4]))
                for ref, alt in sample_counts:
                    m, p, z = calculate_z_score(ref, alt)
                    m_and_z_fh.write(f'\t{m},{z}')
                m_and_z_fh.write('\n')
            else:
                m_and_z_excluded_fh.write('\t'.join(cols[:4]))                
                for ref, alt in sample_counts:
                    m, p, z = calculate_z_score(ref, alt)
                    m_and_z_excluded_fh.write(f'\t{m},{z}')
                m_and_z_excluded_fh.write('\n')

if __name__ == "__main__":
    main()
