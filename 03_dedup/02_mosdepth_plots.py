#!/usr/bin/env python3

#----------------------------------------------------------------------
# Make mosdepth plots per genomic window for the html map and dedup report.
#
# Input:
#   - $DEDUP_RESULTS/depth_metrics/$prefix.regions.bed.gz
# Output:
#   - $DEDUP_RESULTS/depth_metrics/$prefix.depth.500K_bins.png
#----------------------------------------------------------------------

import pandas as pd
import matplotlib.pyplot as plt
import os, sys
sys.path.append("./utils")
from utils import load_config, log

# Load paths from config file
cfg = load_config()
out_dir = cfg["DEDUP_RESULTS"]


def plot(filepath, outplot):
    df = pd.read_csv(filepath, sep='\t', header=None,
            names=["chrom", "start", "end", "depth"])

    # Sort by contig and start, keep only chromo
    df = df[df["chrom"].str.startswith("LG")]
    df.sort_values(by=["chrom", "start"], inplace=True)

    # Get  contig lengths
    contig_lengths = df.groupby("chrom")["end"].max().to_dict()
    contigs_sorted = sorted(contig_lengths.keys(), key=lambda x: int(x.replace("LG", "")))

    # Compute cumulative offsets
    offsets = {}
    current_offset = 0
    for contig in contigs_sorted:
        offsets[contig] = current_offset
        current_offset += contig_lengths[contig]

    # Add cumulative start positions
    df["cum_start"] = df.apply(lambda row: row["start"] + offsets[row["chrom"]], axis=1)

    ylim_max = 500
    fig, ax = plt.subplots(figsize=(12, 5))

    # Plot basic settings
    ax.set_xlim(0, df["cum_start"].max())
    ax.set_ylim(0, ylim_max)
    ax.set_xlabel("Chromosomes")
    ax.set_ylabel("Depth")    
    title = os.path.basename(filepath).split(".")[0]
    ax.set_title(title + " - avg depth per 500k windows")

    # Alternate background colors
    for i, contig in enumerate(contigs_sorted):
        start = offsets[contig]
        end = start + contig_lengths[contig]
        if i % 2 == 0:
            ax.axvspan(start, end, color=(0.9, 0.9, 0.9, 0.5))

    # Plot depth lines
    for contig in contigs_sorted:
        sub_df = df[df["chrom"] == contig]
        ax.plot(sub_df["cum_start"], sub_df["depth"], color="royalblue", linewidth=2)

    # Add contig boundary lines
    for contig in contigs_sorted:
        ax.axvline(x=offsets[contig], color="gray", linestyle="dashed", linewidth=0.8)

    # Add vertical line at end of last contig
    last_contig = contigs_sorted[-1]
    end_last_contig = offsets[last_contig] + contig_lengths[last_contig]
    ax.axvline(x=end_last_contig, color="gray", linestyle="dashed", linewidth=0.8)

    # X-axis labels at contig midpoints
    midpoints = [offsets[c] + contig_lengths[c] / 2 for c in contigs_sorted]
    ax.set_xticks(midpoints)
    ax.set_xticklabels(contigs_sorted, rotation=90, fontsize=8)

    # Horizontal dashed lines in y axis
    ax.axhline(y=400, color="black", linestyle="dashed", linewidth=1)
    ax.axhline(y=300, color="black", linestyle="dashed", linewidth=1)
    ax.axhline(y=200, color="black", linestyle="dashed", linewidth=1)
    ax.axhline(y=100, color="black", linestyle="dashed", linewidth=1)

    # Final layout
    plt.style.use('default')
    plt.tight_layout()
    plt.savefig(outplot, dpi=300)

def main():
    prefix_file = "../data/01_proc_reads/prefixes_new"
    with open(prefix_file) as f:
        prefixes = [line.strip() for line in f]
    
    log(f"=== Mosdepth plots for {len(prefixes)} samples ===")
    for prefix in prefixes:
        i = out_dir + f"/depth_metrics/{prefix}.regions.bed.gz"
        outplot = i.replace("regions.bed.gz", "depth.500K_bins.png")
        plot(i, outplot)
        log(f"done: {outplot}")    
    log(f"=== Mosdepth plots complete ===")
    
if __name__ == "__main__":
    main()

