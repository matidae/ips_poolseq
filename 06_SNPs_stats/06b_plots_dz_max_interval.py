#!/usr/bin/env python3

#----------------------------------------------------------------------
# Plots the distribution of AF changes (dz) between first and last year for each population
#
# Inputs:
#   - dz_max.{prefix}.tsv : per-population file with one dz value per SNP
#
# Outputs:
#   - dz_hist_{prefix}.png : histogram of dz values
#----------------------------------------------------------------------


import os
import matplotlib.pyplot as plt

work_dir = "../../results/06_SNPs_stats"
plot_dir = f"{work_dir}/plots"

color= "#009688"
plt.style.use("ggplot")


def load_dz_file(path):
    dz_vals = []
    with open(path) as fh:
        for line in fh:            
            dz_vals.append(float(line.rstrip()))
    return dz_vals

def plot_dz_distribution(prefix):
    path = f"{work_dir}/dz_max.{prefix}.tsv"
    dz_vals = load_dz_file(path)

    if not dz_vals:
        print(f"No data for {prefix}")
        return

    # Histogram
    plt.figure(figsize=(8, 5))
    plt.hist(dz_vals, bins=50, color=color, edgecolor="black", alpha=0.7)
    plt.axvline(0, color="red", linestyle="--")
    plt.title(f"dz distribution for {prefix}")
    plt.xlabel("dz")
    plt.ylabel("nSNPs")
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/dz_max_hist_{prefix}.png", dpi=300)
    plt.close()

def main():
    os.makedirs(plot_dir, exist_ok=True)
    # Find all dz_max.*.tsv files
    files = [f for f in os.listdir(work_dir) if f.startswith("dz_max.") and f.endswith(".tsv")]
    prefixes = [f.replace("dz_max.", "").replace(".tsv", "") for f in files]

    for pre in prefixes:
        plot_dz_distribution(pre)

if __name__ == "__main__":
    main()
