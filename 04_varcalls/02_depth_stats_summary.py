#!/usr/bin/env python3

#------------------------------------------------------------------------------
# Compute summary statistics of read depths per sample and overall
#
# Input: 
#   - genic_readcounts.tsv: file with CHR, POS, REF, ALT, ADs* for all samples
#   - intergenic_readcounts.tsv 
# Output: 
#   - genic.depth.stats.tsv: summary statistics (mean, median, quantiles, IQR, whiskers) 
#   - intergenic.depth.stats.tsv
#-------------------------------------------------------------------------------

import pandas as pd
from utils import parse_counts

work_dir = "../../results/04_varcalls"

# Input files
genic_counts_in = f"{work_dir}/genic_readcounts.tsv" 
intergenic_counts_in = f"{work_dir}/intergenic_readcounts.tsv" 
# Output files
genic_depth_stats_out = f"{work_dir}/genic_depth_stats.tsv"
intergenic_depth_stats_out = f"{work_dir}/intergenic_depth_stats.tsv"
    
def compute_sample_quantiles(input_file, output_file):
    df = pd.read_csv(input_file, sep="\t", header=0)
    # Extract sample column names
    sample_cols = df.columns[4:]         
    # Iterate over cols and use apply to sum alleles depths
    for sample in sample_cols:
        df[sample] = df[sample].apply(lambda x: sum(parse_counts(x)))
    
    df["TOTAL"] = df[sample_cols].sum(axis=1)
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
    stats.to_csv(output_file, sep="\t")

def main():
    compute_sample_quantiles(genic_counts_in, genic_depth_stats_out)
    compute_sample_quantiles(intergenic_counts_in, intergenic_depth_stats_out)    

if __name__ == "__main__":
    main()