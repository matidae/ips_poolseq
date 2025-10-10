#!/usr/bin/env python3

#----------------------------------------------------------------------
# Plot summary figures based on null variance and divergence calculations.
# In:
#   - null_variance_summary.tsv: summary statistics for dz2 and depth variance
#   - dz2_by_pbin.tsv: dz differences across allele frequency bins
# Outputs:
#   - nv_per_sample.png: bar plot of null variance estimates per sample
#   - avg_zdiff_per_bin.png: replicate differences across allele frequency bins
#   - snp_counts_per_bin.png: SNP counts per allele frequency bin
#----------------------------------------------------------------------

import matplotlib.pyplot as plt
import pandas as pd

work_dir = "../../results/07_null_variance"

#Input files
null_var_in = f"{work_dir}/null_variance_summary.tsv"
dz2_bin_in=f"{work_dir}/dz2_by_pbin.tsv"

color= "#009688"
plt.style.use("ggplot")

# Barplot of null variance estimates per sample
def plot_null_variance_bar(nv_file):
    nv = pd.read_csv(nv_file, sep="\t", usecols=["Sample", "Null_var"])
        # Extract parts for sorting
    nv[["group", "season", "year"]] = nv["Sample"].str.extract(r"^([A-Z]+)_([EL])_(\d{4})")    
    nv["year"] = nv["year"].astype(int) # Convert year to integer     
    season_order = {"E": 0, "L": 1} # Define season order (E < L)
    nv["season_order"] = nv["season"].map(season_order)
    # Sort by group, then season (E before L), then year
    nv = nv.sort_values(by=["group", "season_order", "year"])
   # Create a new column combining group and season to detect changes
    nv["group_season"] = nv["group"] + "_" + nv["season"]
    # Prepare lists for plot data including gaps
    samples_with_gaps = []
    null_vars_with_gaps = []        
    prev_gs = None
    for idx, row in nv.iterrows():
        current_gs = row["group_season"]
        
        # Insert gap if group_season changed (and it's not the first entry)
        if prev_gs is not None and current_gs != prev_gs:
            samples_with_gaps.append("")    # empty label for gap
            null_vars_with_gaps.append(0)   # zero-height bar as gap                    
        samples_with_gaps.append(row["Sample"])
        null_vars_with_gaps.append(row["Null_var"])        
        
        prev_gs = current_gs
    plt.figure(figsize=(14, 6))    
    plt.bar(range(len(samples_with_gaps)), null_vars_with_gaps, color=color)
    plt.xticks(ticks=range(len(samples_with_gaps)), labels=samples_with_gaps, rotation=90)
    plt.ylabel("Null variance estimate")
    plt.tight_layout()
    plt.savefig(f"{work_dir}/nv_per_sample.png")
    plt.close()

# Plot replicate dz averge and SNP counts  per allele frequency bin
def plot_dz_diff_by_freq(freq_file):
    df = pd.read_csv(freq_file, sep="\t")
    plt.figure(figsize=(8, 5))    
    plt.plot(df["pbin"], df["avg_abs_diff"], marker='o', color=color)
    plt.xlabel("Allele frequency bin")
    plt.ylabel("Average absolute difference between replicates")    
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"{work_dir}/avg_zdiff_per_bin.png")
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.bar(df["pbin"], df["count"], color=color)
    plt.xlabel("Allele frequency bin")
    plt.ylabel("Number of SNPs")
    plt.grid(axis='y')
    plt.tight_layout()
    plt.savefig(f"{work_dir}/snp_counts_per_bin.png")
    plt.close()

def more_plots(null_var_in):
    df = pd.read_csv(null_var_in, sep="\t")    
    fig, ax = plt.subplots(figsize=(7, 6))
    # Scatter plot
    ax.scatter(
        df["Depth_var"], df["DZ2_mean"], color=color, s=20)
    # Label dots where DZ2_mean > 0.025
    offset = 0.0003
    for i, (_, row) in enumerate(df[df["DZ2_mean"] > 0.0255].iterrows()):
        ax.text(
            row["Depth_var"],
            row["DZ2_mean"] + (offset if i % 2 == 0 else -offset),  # Added offset 
            str(row["Sample"]), fontsize=10, ha="left", va="bottom")   
    ax.set_xlabel("Mean depth variance")
    ax.set_ylabel("Mean z_diff_sq ")

    # Diagonal line y=x
    lims = [
        min(ax.get_xlim()[0], ax.get_ylim()[0]),
        max(ax.get_xlim()[1], ax.get_ylim()[1])
    ]
    ax.plot(lims, lims, 'k--', alpha=0.5)    
    plt.tight_layout()
    fig.savefig(f"{work_dir}/dz2_vs_depth_var.png", dpi=300)
    plt.close()

def main():
    plot_null_variance_bar(null_var_in)
    plot_dz_diff_by_freq(dz2_bin_in)
    more_plots(null_var_in)

if __name__ == "__main__":
    main()
