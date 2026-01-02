#!/usr/bin/env python3

"""
Generate barplots comparing mapped vs deduplicated coverage depths.
organized by region, country, and time period, with bars grouped by year and replicate.

Input:
    - ../results/03_dedup/summary_table.tsv

Output:
    - ../results/03_dedup/barplots/{Region}{Country}_{Time}.png
"""


import matplotlib.pyplot as plt
import pandas as pd
import os, sys

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
    df["RegionCountrySeason"] = df["RegionCountry"] + "_" + df["Season"]
    # Extract Rep (last char)
    df["Rep"] = df["SeasonRep"].str[-1]
    df["Plot_ID"] = df["RegionCountry"] + "_" + df["Year"].astype(str)
    
    return df

def coverage_barplots(df):
    df_melted = pd.melt(
        df,
        id_vars=["Year", "Rep", "RegionCountrySeason"],
        value_vars=["raw_reads", "qc_reads", "dedup_optical", "dedup_all", "nodedup"],
        var_name="CoverageType",
        value_name="Coverage"
        )

    region_countries = df["RegionCountrySeason"].unique()
    all_years = sorted(df["Year"].unique())
    all_reps = sorted(df["Rep"].unique())    
    all_coverage_types=["raw_reads", "qc_reads", "dedup_optical", "dedup_all", "nodedup"]

    for rc in region_countries:
        fig, ax = plt.subplots(figsize=(12, 5))
        subset = df_melted[df_melted["RegionCountrySeason"] == rc]

        # Create full index to include missing year/rep/coverage combinations
        full_index = pd.MultiIndex.from_product(
            [all_years, all_reps, all_coverage_types],
            names=["Year", "Rep", "CoverageType"]
        )
        full_df = pd.DataFrame(index=full_index).reset_index()
        full_df["RegionCountrySeason"] = rc

        merged = pd.merge(
            full_df, subset,
            on=["RegionCountrySeason", "Year", "Rep", "CoverageType"],
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

        labels = pivot["Year"].astype(str) + "_" + pivot["Rep"]
        x = range(len(pivot))
        width = 0.7

        ax.set_facecolor('#E9E9E9') 
        # Add horizontal reference lines
        for threshold in range(50, 750, 50):
            if threshold % 100 == 0:
                # Thick line every 100
                ax.axhline(threshold, color='grey', linestyle='--', linewidth=0.5, alpha=0.8, zorder=0)
            else:
                # Thin line every 50
                ax.axhline(threshold, color='grey', linestyle='--', linewidth=0.4, alpha=0.6, zorder=0)

        ax.bar([i for i in x], pivot["raw_reads"], width=width, label="Raw reads (est.)", color="sienna")
        ax.bar([i for i in x], pivot["qc_reads"], width=width, label="QC reads (est.)", color="chocolate")
        ax.bar([i for i in x], pivot["nodedup"], width=width, label="Mapped reads", color="goldenrod")
        ax.bar([i for i in x], pivot["dedup_optical"], width=width, label="Optical dedup.", color="khaki")
        ax.bar([i for i in x], pivot["dedup_all"], width=width, label="Full dedup.", color="olive")

        region = rc[0]
        country = rc[1:]
        ax.set_title(f"{rc} depths")
        ax.set_ylabel("Depth")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_ylim(0, 750)
        ax.legend()
        for spine in ax.spines.values():
          spine.set_visible(False)

        plt.tight_layout()
        plt.savefig(f"{out_barplots}/barplot_{region}{country}_depths.png", dpi=300)        
        plt.close()


def main():
    if not os.path.exists(out_barplots):
        os.makedirs(out_barplots)
    df = load_data()
    coverage_barplots(df)


if __name__ == "__main__":
    main()

