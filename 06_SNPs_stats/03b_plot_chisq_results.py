#!/usr/bin/env python3

#----------------------------------------------------------------------
# Plot SNP reliability results based on replicate divergence analysis.
#
# Inputs:
#   - snpdev.m_and_z.tsv: SNP-level coverage, chi-square statistics, and p-values generated from replicate divergence testing.
# Outputs:
#   - pval_histogram.png: Histogram plot of p-values.
#----------------------------------------------------------------------

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np


work_dir = "../../results/06_SNPs_stats"
# Input files
input_file = f"{work_dir}/snpdev_m_and_z.tsv"
# Output files
pval_histogram_plot = f"{work_dir}/pval_histogram.png"
pval_vs_depth_plot = f"{work_dir}/pval_vs_depth.png"

color= "#009688"
plt.style.use("ggplot")

def load_data():
    return pd.read_csv(input_file, sep='\t')

def plot_pval_histogram(df, pval_histogram_plot):
    plt.figure(figsize=(8, 4))
    plt.hist(df['pval'], bins=100, color=color, edgecolor='white')
    plt.axvline(0.01, color='#8B0000', linestyle="--", linewidth=1.5 ,label="p = 0.01")
    ax = plt.gca()  # get current axis
    ax.xaxis.set_major_locator(MultipleLocator(0.1))  # set major ticks every 0.1
    plt.title("Distribution of chi-square p-values")
    plt.xlabel("p-value")
    plt.ylabel("SNPs count")
    plt.tight_layout()
    plt.savefig(pval_histogram_plot)
    plt.close()


def plot_depth_vs_pval(df, depth_pval_plot):
    # Avoid log(0) by replacing zero p-values with the smallest positive float
    df["pval"] = df["pval"].replace(0, np.nextafter(0, 1))
    # Calculate -log10(pval)
    df["-log10(pval)"] = -np.log10(df["pval"])
    df = df.copy() 

    plt.style.use("ggplot")
    plt.figure(figsize=(10, 6))
    plt.scatter(
        df["Depth"],
        df["-log10(pval)"],
        alpha=0.3,
        s=5,
        color=color,
        edgecolors='none'
    )
    plt.ylim(0,5)
    plt.axhline(y=2, color='#8B0000', linestyle='--', linewidth=0.5, label='p = 0.01')
    plt.xlabel("Depth (per SNP)")
    plt.ylabel("-log10(p-value)")
    plt.title("Chi2 p-value vs total depth of SNPs")
    plt.tight_layout()
    plt.savefig(pval_vs_depth_plot, dpi=300)
    plt.close()

def main():
    df = load_data()
    plot_pval_histogram(df, pval_histogram_plot)
    plot_depth_vs_pval(df, pval_vs_depth_plot)

if __name__ == "__main__":
    main()