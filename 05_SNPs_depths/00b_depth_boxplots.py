#!/usr/bin/env python3

#----------------------------------------------------------------------
# Makes boxplots grouped by country, region and season.
# Input:
#   - genic_processed_depths.tsv : table with depths (ref + alt) per SNP, per sample
# Output:
#   - boxplots of genic SNPs per region_time, and for total SNPs
#----------------------------------------------------------------------

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
sys.path.append("./utils")
from plot_style import apply_style, C, style_boxplot, add_year_bands, add_hlines

apply_style(font_size=15)

work_dir = "../results/05_SNPs_depths"
genic_depths_in = f"{work_dir}/genic_processed_depths.tsv"

top_ylim = 400


def plot_depth_boxplots(genic_depths_in, top_ylim):
    years = list(range(2015, 2025))
    reps  = ["a", "b"]

    df = pd.read_csv(genic_depths_in, sep="\t")
    sample_cols = df.columns[4:]

    samples_info = []
    for s in sample_cols:
        parts = s.split("_")
        if len(parts) == 3 and len(parts[1]) >= 2:
            region = parts[0]
            season = parts[1][0]
            rep    = parts[1][1]
            year   = int(parts[2])
            samples_info.append((s, region, season, rep, year))

    groups = {}
    for s, region, season, rep, year in samples_info:
        key = f"{region}_{season}"
        groups.setdefault(key, []).append((s, rep, year))

    for group, samples in groups.items():
        col_labels = [f"{y}{r}" for y in years for r in reps]
        plot_df = pd.DataFrame(index=df.index, columns=col_labels, dtype=float)

        for s, rep, year in samples:
            plot_df[f"{year}{rep}"] = df[s]

        fig, ax = plt.subplots(figsize=(18, 6))

        boxplot_data = []
        for col in col_labels:
            data = plot_df[col].dropna()
            boxplot_data.append(data.values if len(data) > 0 else [np.nan])

        x_pos = np.arange(1, len(col_labels) + 1)

        add_year_bands(ax, years)
        add_hlines(ax, [100, 200, 300])

        box = ax.boxplot(
            boxplot_data, positions=x_pos, patch_artist=True,
            widths=0.55, showfliers=False, zorder=3,
        )
        style_boxplot(box)

        # Median annotations
        for i, median_line in enumerate(box["medians"]):
            med = median_line.get_ydata()[0]
            if pd.notnull(med) and med > 0:
                ax.text(
                    x_pos[i], med + 7, str(int(round(med))),
                    ha="center", va="bottom", fontsize=11,
                    color=C["median"], fontweight="bold",
                )

        ax.set_xlim(0.2, len(col_labels) + 0.8)
        ax.set_ylim(-5, top_ylim)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(col_labels, rotation=45, ha="right")
        ax.set_ylabel("Sequencing depth", labelpad=8)
        ax.set_title(group.replace("_", " · "), pad=14)

        legend_handles = [
            mpatches.Patch(facecolor=C["teal_light"], edgecolor=C["teal_dark"], label="IQR"),
            plt.Line2D([0], [0], color=C["median"], linewidth=2, label="Median"),
        ]
        ax.legend(handles=legend_handles, loc="upper left")

        fig.tight_layout()
        fig.savefig(f"{work_dir}/boxplots/{group}_depth_boxplots.png")
        plt.close(fig)


def plot_total_depths(genic_depths_in):
    genic = pd.read_csv(genic_depths_in, sep="\t")
    data  = [genic["TOTAL"]]

    fig, ax = plt.subplots(figsize=(5, 6))

    bplot = ax.boxplot(
        data, tick_labels=["Genic SNPs"], patch_artist=True,
        widths=0.25, showfliers=False, zorder=3,
    )
    style_boxplot(bplot)

    med = data[0].median()
    if pd.notna(med):
        ax.text(
            1, med + 5, str(int(round(med))),
            ha="center", va="bottom", fontsize=11,
            color=C["median"], fontweight="bold",
        )

    ax.set_ylabel("Total depth per SNP", labelpad=8)
    ax.set_title("Total depth distribution", pad=12)

    fig.tight_layout()
    fig.savefig(f"{work_dir}/boxplots/total_depth_boxplots.png")
    plt.close(fig)


def main(genic_depths_in):
    os.makedirs(f"{work_dir}/boxplots", exist_ok=True)
    plot_depth_boxplots(genic_depths_in, top_ylim)
    plot_total_depths(genic_depths_in)


if __name__ == "__main__":
    main(genic_depths_in)