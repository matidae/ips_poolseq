#!/usr/bin/env python3

"""
Generate barplots comparing mapped vs deduplicated coverage depths.
organized by region, country, and time period, with bars grouped by year and replicate.

Input:
    - ../../results/03_dedup/summary_table.tsv

Output:
    - ../../results/03_dedup/barplots/{Region}{Country}_{Time}.png
"""
#----------------------------------------------------------------------

import matplotlib.pyplot as plt
import pandas as pd
import os

out_dir="../results/03_dedup"
out_barplots = out_dir + "/barplots"

ips_genome_size = 224977219

def load_data():
    df = pd.read_csv(out_dir + "/summary_table.tsv", sep="\t")
    df["Cov_mapped"] = (df["Mapped_reads"] * 1000000 * df["Length"]) / ips_genome_size
    df["Cov_dedup"] = df["Mean_depth"]
    df["Plot_ID"] = df["Region"] + "_" + df["Country"] + "_" + df["Time"]
    return df

def coverage_barplots(df):

    df_melted = pd.melt(
        df,
        id_vars=["Year", "Rep", "Plot_ID"],
        value_vars=["Cov_mapped", "Cov_dedup"],
        var_name="CoverageType",
        value_name="Coverage"
    )
    plot_ids = df["Plot_ID"].unique()
    all_years = sorted(df["Year"].unique())
    all_reps = sorted(df["Rep"].unique())
    all_coverage_types = ["Cov_mapped", "Cov_dedup"]

    for plot_id in plot_ids:
        fig, ax = plt.subplots(figsize=(10, 4))

        # Subset and full merge
        subset = df_melted[df_melted["Plot_ID"] == plot_id]
        full_index = pd.MultiIndex.from_product(
            [all_years, all_reps, all_coverage_types],
            names=["Year", "Rep", "CoverageType"]
        )
        full_df = pd.DataFrame(index=full_index).reset_index()
        full_df["Plot_ID"] = plot_id
        merged = pd.merge(
            full_df, subset,
            on=["Plot_ID", "Year", "Rep", "CoverageType"],
            how="left"
        )
        merged["Coverage"] = merged["Coverage"].fillna(0)

        # Pivot for bar plotting
        pivot = merged.pivot_table(
            index=["Year", "Rep"],
            columns="CoverageType",
            values="Coverage"
        ).reset_index()

        pivot = pivot.sort_values(["Year", "Rep"])
        labels = pivot["Year"].astype(str) + "_" + pivot["Rep"]
        x = range(len(pivot))
        width = 0.35

        ax.bar([i - width/2 for i in x], pivot["Cov_mapped"], width=width, label="No_dedup", color="peru")
        ax.bar([i + width/2 for i in x], pivot["Cov_dedup"], width=width, label="Dedup", color="darkkhaki")

        country = plot_id.split("_")[1]
        region = plot_id.split("_")[2]
        time = plot_id.split("_")[0]
        ax.set_title(f"{region}{country}_{time}")
        ax.set_ylabel("Depth")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_ylim(0, 400)
        ax.axhline(100, color='gray', linestyle='--')
        ax.axhline(200, color='gray', linestyle='--')
        ax.axhline(300, color='gray', linestyle='--')
        ax.legend()

        plt.tight_layout()    
        plt.savefig(f"{out_barplots}{region}{country}_{time}.png", dpi=300)
        plt.close()

def main():
    if not os.path.exists(out_barplots):
        os.makedirs(out_barplots)
    df = load_data()
    coverage_barplots(df)


if __name__ == "__main__":
    main()

