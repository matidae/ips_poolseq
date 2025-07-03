#!/usr/bin/env python3

#------------------------------------------------------------------------------
# Compute summary statistics of read depths per sample and overall
#
# Input: 
#   - genic_readcounts.tsv: file with CHR, POS, REF, ALT, ADs* for all samples
#   - intergenic_readcounts.tsv 
# Output: 
#   - [inter]genic_depth_stats.tsv: summary statistics (mean, median, quantiles, IQR, whiskers)
#   - [inter]genic_processed_depths.tsv: sum of alt+ref per SNP and sample 
#-------------------------------------------------------------------------------

import os
import pandas as pd
from utils import parse_counts

work_dir = "../../results/04_varcalls"
out_dir = "../../results/05_SNPs_depths"
# Input files
genic_counts_in = f"{work_dir}/genic_readcounts.tsv" 
intergenic_counts_in = f"{work_dir}/intergenic_readcounts.tsv" 
# Output files
genic_depth_stats_out = f"{out_dir}/genic_depth_stats.tsv"
intergenic_depth_stats_out = f"{out_dir}/intergenic_depth_stats.tsv"
genic_processed_depths_out = f"{out_dir}/genic_processed_depths.tsv"
intergenic_processed_depths_out = f"{out_dir}/intergenic_processed_depths.tsv"
    
def compute_sample_quantiles(counts_in, depth_stats_out, processed_depths_out):
    df = pd.read_csv(counts_in, sep="\t", header=0)
    # Extract sample column names
    sample_cols = df.columns[4:]
    # Iterate over cols and use apply to sum alleles depths
    for sample in sample_cols:
        df[sample] = df[sample].apply(lambda x: sum(parse_counts(x)))
    
    df["TOTAL"] = df[sample_cols].sum(axis=1)
    processed_df = df[["#CHROM", "POS", "REF", "ALT"] + list(sample_cols) + ["TOTAL"]]
    processed_df.to_csv(processed_depths_out, sep="\t", index=False)
    
    sample_cols = list(sample_cols) + ["TOTAL"]

    stats = df[sample_cols].describe().T
    stats["mean"] = stats["mean"].round(1)
    stats["std"] = stats["std"].round(1)

    stats["iqr"] = stats["75%"] - stats["25%"]    
    stats["whislo"] = (stats["25%"] - 1.5 * stats["iqr"]).clip(lower=0).round(1)
    stats["whishi"] = (stats["75%"] + 1.5 * stats["iqr"]).round(1)
       
    stats["median"] = df[sample_cols].median().round(1)
    stats["median_50"] = (stats["median"].astype(float)*0.5).round(1)
    stats["median_150"] = (stats["median"].astype(float)*1.5).round(1)

    # Cast columns back to int
    int_cols = stats.columns.difference(["mean", "std", "median", "whislo", "whishi","median_50", "median_150"])
    stats[int_cols] = stats[int_cols].astype(int)
    stats.index.name = "sample" 
    stats.to_csv(depth_stats_out, sep="\t")

def main():
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    compute_sample_quantiles(genic_counts_in, genic_depth_stats_out, genic_processed_depths_out)
    compute_sample_quantiles(intergenic_counts_in, intergenic_depth_stats_out, intergenic_processed_depths_out)    

if __name__ == "__main__":
    main()