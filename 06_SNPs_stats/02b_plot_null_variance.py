#!/usr/bin/env python3

#----------------------------------------------------------------------
# Plot summary figures based on null variance and divergence calculations.
# In:
#   - null_variance_summary.tsv: summary statistics for dz² and depth variance
#   - dz2_by_pbin.tsv: dz differences across allele frequency bins
# Outputs:
#   - nv_per_sample.png: bar plot of null variance estimates per sample
#   - avg_zdiff_per_bin.png: replicate differences across allele frequency bins
#   - snp_counts_per_bin.png: SNP counts per allele frequency bin
#----------------------------------------------------------------------

import matplotlib.pyplot as plt
import pandas as pd


work_dir = "../../results/06_SNPs_stats"
#Input files
null_var_in = f"{work_dir}/null_variance_summary.tsv"
dz2_bin_in=f"{work_dir}/dz2_by_pbin.tsv"

color= "#009688"
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
    plt.style.use("ggplot")
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
    plt.style.use("ggplot")
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

import seaborn as sns

# Load data
def more_plots(null_var_in):
    df = pd.read_csv(null_var_in, sep="\t")

    sns.set(style="whitegrid")

    fig, axs = plt.subplots(1, 2, figsize=(14, 6))

   # Scatter plot DZ2_mean vs ReadDepth_var (no change)
    sns.scatterplot(data=df, x="ReadDepth_var", y="DZ2_mean", hue="Sample", ax=axs[0], palette="tab20", legend=False)
    axs[0].set_title("DZ2_mean vs ReadDepth_var")
    axs[0].set_xlabel("Mean Read Depth Variance (ReadDepth_var)")
    axs[0].set_ylabel("Mean Squared Difference (DZ2_mean)")

    lims = [
        min(axs[0].get_xlim()[0], axs[0].get_ylim()[0]),
        max(axs[0].get_xlim()[1], axs[0].get_ylim()[1])
    ]
    axs[0].plot(lims, lims, 'k--', alpha=0.5)
    axs[0].set_xlim(lims)
    axs[0].set_ylim(lims)

    # Barplot without palette warning & fix tick labels rotation
    sns.barplot(data=df, x="Sample", y="RDprop", ax=axs[1])
    axs[1].set_title("Proportion Variance Explained by Read Depth (RDprop)")
    axs[1].set_xlabel("Sample")
    axs[1].set_ylabel("RDprop")
    plt.setp(axs[1].get_xticklabels(), rotation=45, ha='right')  # <-- correct way to rotate and align x labels


    plt.tight_layout()

    # Save figure as PNG
    fig.savefig("null_variance_plots.png", dpi=300)

    # If you want to see it locally, uncomment this:
    # plt.show()


def main():
    #plot_null_variance_bar(null_var_in)
    #plot_dz_diff_by_freq(dz2_bin_in)
    more_plots(null_var_in)

if __name__ == "__main__":
    main()



