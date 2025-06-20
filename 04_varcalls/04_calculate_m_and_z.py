#!/usr/bin/env python3

#----------------------------------------------------------------------
# Calculates read depth (m) and Fisher-Ford transformation (z) for allele frequencies
# Removes SNPs with depth values of outliers >1.5 median and <0.5 median
# Run: ./04_calculate_m_and_z.py
# In: tables of statistics calculated by ./02_get_SNPs_depth.py
# Out: tsv file with m and z values for each SNP
#----------------------------------------------------------------------

from math import sqrt, asin

min_depth = 5508 # 1/2 the median
max_depth = 16524 # 1.5 the median

wd = "../../results/04_varcalls"

def parse_counts(count_str):    
    ref, alt = map(int, count_str.split(','))
    return ref, alt

# Calculates z Fisher transformation
def calculate_z_score(ref, alt):    
    m = ref + alt
    if m == 0:
        return 0, 'NA', 'NA'
    p = ref / m  # Calculates p as ref/(ref+alt)
    z = 2 * asin(sqrt(p)) # Fisher transform of p
    return m, f"{p:.4f}", f"{z:.4f}"

with open(f"{wd}/genic_readcounts.tsv", 'r') as infile, open(f"{wd}/genic_m_and_z.tsv", 'w') as out:
    for line in infile:
        cols = line.strip().split('\t')
        # Write header
        if cols[0] == "#CHROM":
            out.write(line)
            continue
        # Parse sample read counts
        sample_counts = [parse_counts(x) for x in cols[4:]]
        total_depth = sum(ref + alt for ref, alt in sample_counts)

        # Filter by SNPs position depth sum
        if min_depth <= total_depth <= max_depth:
            out.write('\t'.join(cols[:4]))
            for ref, alt in sample_counts:
                m, p, z = calculate_z_score(ref, alt)
                out.write(f'\t{m},{z}')
            out.write('\n')


        



