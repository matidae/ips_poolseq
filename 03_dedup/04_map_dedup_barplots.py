#!/usr/bin/env python3

"""
Generate barplots comparing mapped vs deduplicated coverage depths.
organized by region, country, and time period, with bars grouped by year and replicate.

Input:
    - ../results/03_dedup/summary_table.tsv

Output:
    - ../results/03_dedup/barplots/{Region}{Country}_depths.png
"""


import matplotlib.pyplot as plt
import pandas as pd
import os
from matplotlib.patches import Patch

out_dir="../results/03_dedup"
out_barplots = out_dir + "/barplots"

ips_genome_size = 224977219

def load_data():
    df = pd.read_csv(out_dir + "/summary_table.tsv", sep="\t")  

    df["nodedup"] = df["Mean_depth_nodedup"]
    df["dedup_all"] = df["Mean_depth_dedup"]
    df["dedup_optical"] = df["Mean_depth_dedup_alt"]
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
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    # Melt coverage values for plotting
    df_melted = pd.melt(
        df,
        id_vars=["Year", "Season", "Rep", "RegionCountry"],
        value_vars=["raw_reads", "qc_reads", "dedup_optical", "dedup_all", "nodedup"],
        var_name="CoverageType",
        value_name="Coverage"
    )

    region_countries = df["RegionCountry"].unique()
    all_years = sorted(df["Year"].unique())
    all_reps = sorted(df["Rep"].unique())
    all_seasons = ["E", "L"]
    all_coverage_types = ["raw_reads", "qc_reads", "dedup_optical", "dedup_all", "nodedup"]

    for rc in region_countries:
        fig, ax = plt.subplots(figsize=(12, 5))
        subset = df_melted[df_melted["RegionCountry"] == rc]

        # Complete combinations to include missing year/season/rep/coverage
        full_index = pd.MultiIndex.from_product(
            [all_years, all_seasons, all_reps, all_coverage_types],
            names=["Year", "Season", "Rep", "CoverageType"]
        )
        full_df = pd.DataFrame(index=full_index).reset_index()
        full_df["RegionCountry"] = rc

        merged = pd.merge(
            full_df, subset,
            on=["RegionCountry", "Year", "Season", "Rep", "CoverageType"],
            how="left"
        )
        merged["Coverage"] = merged["Coverage"].fillna(0)

        # Pivot for plotting
        pivot = merged.pivot_table(
            index=["Year", "Season", "Rep"],
            columns="CoverageType",
            values="Coverage"
        ).reset_index()
        pivot = pivot.sort_values(["Year", "Season", "Rep"])

        # Compute clustered x positions
        x_positions = []
        label_positions = []
        label_texts = []
        current_x = 0
        cluster_offset = 0.25  # distance between replicates
        group_gap = 0.5        # gap between Season/Year

        grouped = pivot.groupby(["Year", "Season"])
        for (year, season), group in grouped:
            reps = len(group)
            # positions for each replicate
            positions = [current_x + i*cluster_offset for i in range(reps)]
            x_positions.extend(positions)
            # position for label (end at last replicate)
            label_positions.append(positions[-1])
            # only add label if any coverage > 0
            if group[["raw_reads","qc_reads","nodedup","dedup_optical","dedup_all"]].sum().sum() > 0:
                label_texts.append(f"{rc}_{season}_{year}")
            else:
                label_texts.append(" ")
            current_x += reps*cluster_offset + group_gap

        width = 0.2

        # Background and reference lines
        ax.set_facecolor('#E9E9E9')
        for threshold in range(50, 750, 50):
            lw = 0.5 if threshold % 100 == 0 else 0.4
            alpha = 0.8 if threshold % 100 == 0 else 0.6
            ax.axhline(threshold, color='grey', linestyle='--', linewidth=lw, alpha=alpha, zorder=0)

        # Red line at 150
        ax.axhline(100, color='red', linestyle='--', linewidth=1, zorder=1)

        # Plot bars
        bar_dict = {}
        bar_dict["raw_reads"] = ax.bar(x_positions, pivot["raw_reads"], width=width, color="sienna")
        bar_dict["qc_reads"] = ax.bar(x_positions, pivot["qc_reads"], width=width, color="chocolate")
        bar_dict["nodedup"] = ax.bar(x_positions, pivot["nodedup"], width=width, color="goldenrod")
        bar_dict["dedup_optical"] = ax.bar(x_positions, pivot["dedup_optical"], width=width, color="khaki")
        bar_dict["dedup_all"] = ax.bar(x_positions, pivot["dedup_all"], width=width, color="olive")

        # Add optical dedup values as text inside bars
        for x, val in zip(x_positions, pivot["dedup_all"]):
            if val > 0:  # only label if coverage > 0
                ax.text(x, val/2, f"{val:.0f}", ha="center", va="center", fontsize=6, rotation=90,  weight="bold", color="white")


        # Conditional alpha per Season group
        group_alpha = pivot.groupby(["Year","Season"])["dedup_all"].transform(lambda x: (x < 100).any())
        n = len(pivot)
        for i, patch in enumerate(ax.patches):
            rep_idx = i % n
            if group_alpha.iloc[rep_idx]:
                patch.set_alpha(0.3)

        # Legend: fully opaque
        legend_handles = [
            Patch(facecolor="sienna", label="Sequenced reads (est.)", alpha=1.0),
            Patch(facecolor="chocolate", label="QC-filtered reads (est.)", alpha=1.0),
            Patch(facecolor="goldenrod", label="Aligned reads", alpha=1.0),
            Patch(facecolor="khaki", label="Optical duplicates removed", alpha=1.0),
            Patch(facecolor="olive", label="All duplicates removed", alpha=1.0)
        ]
        ax.legend(handles=legend_handles)

        # X-axis labels
        ax.set_xticks(label_positions)
        ax.set_xticklabels(label_texts, rotation=45, ha="right", fontsize=7)

        region = rc[0]
        country = rc[1:]
        ax.set_title(f"{rc} depths")
        ax.set_ylabel("Depth")
        ax.set_ylim(0, 750)
        for spine in ax.spines.values():
            spine.set_visible(False)

        plt.tight_layout()
        plt.savefig(f"{out_barplots}/barplot_{region}{country}_depths_100.png", dpi=300)
        plt.close()



def main():
    if not os.path.exists(out_barplots):
        os.makedirs(out_barplots)
    df = load_data()
    coverage_barplots(df)


if __name__ == "__main__":
    main()

