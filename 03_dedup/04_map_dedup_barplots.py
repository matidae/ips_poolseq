#!/usr/bin/env python3

#----------------------------------------------------------------------
# Generate barplots comparing mapped vs deduplicated coverage depths.
# organized by region, country, and time period, with bars grouped by year and replicate.
#
# Input:
#    - $DEDUP_RESULTS/summary_table.tsv
#
# Output:
#    - $DEDUP_RESULTS/barplots/barplot_{region}{country}_{season}_depths_100.png
#----------------------------------------------------------------------

import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch
sys.path.append("./utils")
from plot_style import apply_style, C
from utils import load_config, log

apply_style()

cfg = load_config()
out_dir = cfg["DEDUP_RESULTS"]
out_barplots = out_dir + "/barplots"

ips_genome_size = int(cfg["GENOME_SIZE"])

def load_data():
    df = pd.read_csv(out_dir + "/summary_table.tsv", sep="\t")  

    df["nodedup"] = df["Depth_nodedup"]
    df["dedup_all"] = df["Depth_dedup"]
    df["dedup_optical"] = df["Depth_dedup_alt"]
    df["raw_reads"] = (df["Raw_reads"] * 1e6 * 150) / ips_genome_size
    df["qc_reads"] = (df["QC_reads"] * 1e6 * df["Length"]) / ips_genome_size

    num_cols = ["raw_reads", "qc_reads", "dedup_optical", "dedup_all","nodedup"]
    df[num_cols] = df[num_cols].apply(pd.to_numeric, errors="coerce")

    # Extract Region, Country, Season, Year, Rep from Idn 
    df[["RegionCountry", "SeasonRep", "Year"]] = df["Idn"].str.split("_", expand=True)
    df["Region"] = df["RegionCountry"].str[0]
    df["Country"] = df["RegionCountry"].str[1:]
    df["Season"] = df["SeasonRep"].str[0]  # E or L    
    # Extract Rep (last char)
    df["Rep"] = df["SeasonRep"].str[-1]
    df["Plot_ID"] = df["RegionCountry"] + "_" + df["Year"].astype(str)
    
    return df

def coverage_barplots(df):
    # Melt coverage values for plotting
    df_melted = pd.melt(
        df,
        id_vars=["Year", "Season", "Rep", "RegionCountry"],
        value_vars=["raw_reads", "qc_reads", "dedup_optical", "dedup_all", "nodedup"],
        var_name="CoverageType",
        value_name="Coverage"
    )

    region_countries = df["RegionCountry"].unique()
    all_years = [str(y) for y in range(2015, 2026)]
    all_reps = sorted(df["Rep"].unique())
    all_seasons = ["E", "L"]
    all_coverage_types = ["raw_reads", "qc_reads", "dedup_optical", "dedup_all", "nodedup"]

    for rc in region_countries:
        for season in all_seasons:
            subset = df_melted[
                (df_melted["RegionCountry"] == rc) &
                (df_melted["Season"] == season)
            ]

            full_index = pd.MultiIndex.from_product(
                [all_years, all_reps, all_coverage_types],
                names=["Year", "Rep", "CoverageType"]
            )
            full_df = pd.DataFrame(index=full_index).reset_index()
            full_df["RegionCountry"] = rc
            full_df["Season"] = season

            merged = pd.merge(
                full_df, subset,
                on=["RegionCountry", "Season", "Year", "Rep", "CoverageType"],
                how="left"
            )
            merged["Coverage"] = merged["Coverage"].fillna(0)

            # Pivot for plotting
            pivot = merged.pivot_table(
                index=["Year", "Rep"],
                columns="CoverageType",
                values="Coverage"
            ).reset_index()
            pivot = pivot.sort_values(["Year", "Rep"])

            # Compute clustered x positions
            x_positions = []
            label_positions = []
            label_texts = []
            current_x = 0
            width = 0.3
            cluster_offset = 0.35  # distance between replicates
            group_gap = 0.5        # gap between Season/Year

            n_years = len(pivot["Year"].unique())
            n_bars = (pivot["dedup_all"] > 0).sum()
            fig_width = 20
            fig, ax = plt.subplots(figsize=(fig_width, 5))

            grouped = pivot.groupby(["Year"])
            for year, group in grouped:
                reps = len(group)
                positions = [current_x + i*cluster_offset for i in range(reps)]
                x_positions.extend(positions)
                # one label per year, centered between replicates
                label_positions.append((positions[0] + positions[-1]) / 2)
                if group[["raw_reads","qc_reads","nodedup","dedup_optical","dedup_all"]].sum().sum() > 0:
                    label_texts.append(f"{year[0]}")
                else:
                    label_texts.append(" ")
                current_x += reps*cluster_offset + group_gap

            # Background and reference lines
            ax.set_facecolor(C["grey_bg"])
            for threshold in range(50, 750, 50):
                lw = 0.5 if threshold % 100 == 0 else 0.4
                alpha = 0.8 if threshold % 100 == 0 else 0.6
                ax.axhline(threshold, color='grey', linestyle='--', linewidth=lw, alpha=alpha, zorder=0)

            # Red line at 100
            ax.axhline(100, color=C["rust"], linestyle='--', linewidth=1, zorder=1)

            # Plot bars
            bar_dict = {}
            bar_dict["raw_reads"]     = ax.bar(x_positions, pivot["raw_reads"],     width=width, color=C["pipe"][0])
            bar_dict["qc_reads"]      = ax.bar(x_positions, pivot["qc_reads"],      width=width, color=C["pipe"][1])
            bar_dict["nodedup"]       = ax.bar(x_positions, pivot["nodedup"],       width=width, color=C["pipe"][2])
            bar_dict["dedup_optical"] = ax.bar(x_positions, pivot["dedup_optical"], width=width, color=C["pipe"][3])
            bar_dict["dedup_all"]     = ax.bar(x_positions, pivot["dedup_all"],     width=width, color=C["pipe"][4])

            # Add optical dedup values as text inside bars
            for x, val in zip(x_positions, pivot["dedup_all"]):
                if val > 0:  # only label if coverage > 0
                    ax.text(x, val/2, f"{val:.0f}", ha="center", va="center", fontsize=14, rotation=90, weight="bold", color="white")

            # Conditional alpha per Season group
            group_alpha = pivot.groupby(["Year"])["dedup_all"].transform(lambda x: (x < 100).any())
            n = len(pivot)
            for i, patch in enumerate(ax.patches):
                rep_idx = i % n
                if group_alpha.iloc[rep_idx]:
                    patch.set_alpha(0.3)

            # Legend: fully opaque
            legend_handles = [
                Patch(facecolor=C["pipe"][0], label="Sequenced reads (est.)", alpha=1.0),
                Patch(facecolor=C["pipe"][1], label="QC-filtered reads (est.)", alpha=1.0),
                Patch(facecolor=C["pipe"][2], label="Aligned (pre-dedup)", alpha=1.0),
                Patch(facecolor=C["pipe"][3], label="After optical dedup", alpha=1.0),
                Patch(facecolor=C["pipe"][4], label="Depth excl. all flagged dups", alpha=1.0)
            ]
            ax.legend(handles=legend_handles, fontsize=14)

            # X-axis labels
            ax.set_xticks(label_positions)
            ax.set_xticklabels(label_texts, rotation=0, ha="center", fontsize=14, fontweight="bold")

            region = rc[0]
            country = rc[1:]
            ax.set_title(f"{rc}_{season} depths")
            ax.set_ylabel("Depth", fontsize=14)
            ax.tick_params(axis='y', labelsize=14)
            ax.set_ylim(0, 750)
            for spine in ax.spines.values():
                spine.set_visible(False)

            plt.tight_layout()
            plt.savefig(f"{out_barplots}/barplot_{region}{country}_{season}_depths_100.png")
            plt.close()


def main():
    log("=== Barplots start ===")
    if not os.path.exists(out_barplots):
        os.makedirs(out_barplots)
    df = load_data()
    log(f"loaded {len(df)} samples from summary table")
    coverage_barplots(df)
    log(f"done: {out_barplots}")
    log("=== Barplots complete ===")


if __name__ == "__main__":
    main()