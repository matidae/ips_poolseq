#!/usr/bin/env python3

#----------------------------------------------------------------------
# Calculates basic statistics of readcounts for each sample and for the total sum
# Run: ./02_get_SNPs_depth.py
# In: genic_readcounts.tsv and intergenic_readcounts.tsv (output of 01_filter_vcf.sh)
# Out: tables of statistics (median, quantiles, IQR) for further calculations and for making boxplots
#----------------------------------------------------------------------

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

wd = "../../results/04_varcalls"

def sum_depth_alleles(i):
    depth = int(i.split(",")[0]) + int(i.split(",")[1])
    return(depth)

def compute_sample_quantiles(input_file, output_file):
    df = pd.read_csv(input_file, sep="\t", header=0)
    # Extract sample column names
    sample_cols = df.columns[4:]         
    # Iterate over cols and use apply to sum alleles depths
    for sample in sample_cols:
        df[sample] = df[sample].apply(sum_depth_alleles)
    
    df["TOTAL"] = df[sample_cols].sum(axis=1)
    sample_cols = list(sample_cols) + ["TOTAL"]

    stats = df[sample_cols].describe().T
    stats["mean"] = stats["mean"].map("{:.2f}".format)
    stats["std"] = stats["std"].map("{:.2f}".format)

    stats["iqr"] = stats["75%"] - stats["25%"]    
    stats["whislo"] = (stats["25%"] - 1.5 * stats["iqr"]).clip(lower=0)
    stats["whishi"] = stats["75%"] + 1.5 * stats["iqr"]
   
    stats["median"] = df[sample_cols].median().map("{:.1f}".format)
    stats["whislo"] = stats["whislo"].round(1)
    stats["whishi"] = stats["whishi"].round(1)

    # Cast int-type columns back to int
    int_cols = stats.columns.difference(["mean", "std", "median", "whislo", "whishi"])
    stats[int_cols] = stats[int_cols].astype(int)    
    stats.to_csv(output_file, sep="\t", float_format="%.2f")

compute_sample_quantiles(f"{wd}/genic_readcounts.tsv", f"{wd}/genic.snp_total_depths.stats.tsv")
compute_sample_quantiles(f"{wd}/intergenic_readcounts.tsv", f"{wd}/intergenic.snp_total_depths.stats.tsv")