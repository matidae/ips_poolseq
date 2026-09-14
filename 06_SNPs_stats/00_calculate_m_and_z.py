#!/usr/bin/env python3

#-------------------------------------------------------------------------------
# Calculates read depth (m) and Fisher-Ford transformation (z) for allele frequencies
# Removes per-sample depth values that are outliers (>1.5x or <0.5x that sample's own median)
# Filtering is done per-sample, not on the summed depth across all samples
#
# Input: 
#   - genic_depth_stats.tsv: summary statistics (mean, median, quantiles, IQR, whiskers) per sample
#   - genic_readcounts.tsv: file with CHR, POS, REF, ALT, ADs* for all samples
# Output: 
#   - genic_m_and_z.tsv: tsv file with m and z values for every SNP; per-sample values are 'NA'
#     where that sample's own depth at that position falls outside its [min_depth, max_depth] range
#-------------------------------------------------------------------------------
import sys
import os
from math import sqrt, asin
sys.path.append("./utils")
from utils import parse_counts, load_depth_threshold, load_config, log

cfg = load_config()

input_dir = cfg["VARCALL_RESULTS"]
work_dir = cfg["STATS_RESULTS"]

# Input files
depth_stats_in = f"{input_dir}/genic_depth_stats.tsv"
readcounts_in = f"{input_dir}/genic_readcounts.tsv"
# Output file
m_and_z_out = f"{work_dir}/genic_m_and_z.tsv"

min_max_by_sample = load_depth_threshold()  # dict of {sample_name: (min_depth, max_depth)}

# Calculates Fisher-Ford arcsine square root transformation of allele frequencies
# Returns NA for m/z if this sample's own depth at this SNP is outside its own [min_depth, max_depth]
def calculate_z_score(ref, alt, sample_name):
    m = ref + alt
    min_depth, max_depth = min_max_by_sample[sample_name]
    if m == 0 or not (min_depth <= m <= max_depth):
        return m, 'NA', 'NA'
    p = ref / m  # Calculates p as ref/(ref+alt)
    z = 2 * asin(sqrt(p)) # Fisher-Ford transform of p
    return m, f"{p:.4f}", f"{z:.4f}"

def main():
    log("=== Start: calculate m and z (per-sample depth masking) ===")
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    with open(readcounts_in, 'r') as readcounts_fh, open(m_and_z_out, 'w') as m_and_z_fh:
        header = next(readcounts_fh)
        # Get sample column names from header, same order as cols[4:] on every data line
        sample_names = header.strip().split('\t')[4:]        
        log(f"Loaded depth thresholds for {len(sample_names)} samples, matched against readcounts header")
        m_and_z_fh.write(header)
        n_rows = 0
        for line in readcounts_fh:        
            cols = line.strip().split('\t')            
            # Parse sample read counts
            sample_counts = [parse_counts(x) for x in cols[4:]]
            # No row-level filter: every SNP is kept, masking happens per-sample inside calculate_z_score
            m_and_z_fh.write('\t'.join(cols[:4]))
            for sample_name, (ref, alt) in zip(sample_names, sample_counts):
                m, p, z = calculate_z_score(ref, alt, sample_name)
                m_and_z_fh.write(f'\t{m},{z}')
            m_and_z_fh.write('\n')
            n_rows += 1
    log(f"Wrote {n_rows} SNPs to {m_and_z_out}")
    log("=== Done ===")

if __name__ == "__main__":
    main()