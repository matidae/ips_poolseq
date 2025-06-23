#!/usr/bin/env python3

#----------------------------------------------------------------------
# Plot summary figures based on null variance and divergence calculations.
# In:
#   - nv.genic_m_and_z: estimated nucleotide variation values
#   - null_variance_summary.tsv: summary statistics for dz² and depth variance
#   - pairdz.byfreq.txt: dz differences across allele frequency bins
# Outputs:
#   - nv_per_sample.png: bar plot of null variance estimates per sample
#   - avg_zdiff_per_bin.png: replicate differences across allele frequency bins
#   - snp_counts_per_bin.png: SNP counts per allele frequency bin
#----------------------------------------------------------------------

import matplotlib.pyplot as plt
import pandas as pd


wd = "../../results/04_varcalls"
nv_file=f"{wd}/nv.genic_m_and_z"
infile2=f"{wd}/null_variance_summary.tsv"
freq_file = f"{wd}/pairdz.byfreq.txt"

#color = 'cornflowerblue'
color= "#009688"
# Barplot of null variance estimates per sample
def plot_null_variance_bar(nv_file):
    nv = pd.read_csv(nv_file, sep="\t", header=None, names=["sample", "nv"])
    plt.figure(figsize=(12,6))
    plt.style.use("ggplot")
    plt.bar(nv["sample"], nv["nv"], color=color)
    plt.xticks(rotation=90)
    plt.ylabel("Null variance estimate (nv_init)")
    plt.tight_layout()
    plt.savefig(f"{wd}/nv_per_sample.png")
    plt.close()

# Plot replicate dz differences and counts per allele frequency bin
def plot_dz_diff_by_freq(freq_file):
    df = pd.read_csv(freq_file, sep="\t", header=None, names=["pcat", "count", "avg_abs_diff"])
    plt.figure(figsize=(8, 5))
    plt.style.use("ggplot")
    plt.plot(df["pcat"], df["avg_abs_diff"], marker='o', color=color)
    plt.xlabel("Allele frequency bin (pcat)")
    plt.ylabel("Average absolute difference between replicates")    
    plt.grid(True)    
    plt.tight_layout()
    plt.savefig(f"{wd}/avg_zdiff_per_bin.png")
    plt.close()

    # Bar plot of SNP counts per bin
    plt.figure(figsize=(8, 5))
    plt.bar(df["pcat"], df["count"], color=color)
    plt.xlabel("Allele frequency bin (pcat)")
    plt.ylabel("Number of SNPs")
    plt.grid(axis='y')
    plt.tight_layout()
    plt.savefig(f"{wd}/snp_counts_per_bin.png")
    plt.close()

plot_null_variance_bar(nv_file)    
plot_dz_diff_by_freq(freq_file)




