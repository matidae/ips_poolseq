#!/usr/bin/env python3

#----------------------------------------------------------------------
# Makes boxplots grouped by country, region and season.
# Input: 
#   - [inter]genic_processed_depths.tsv : table with depths (ref + alt) per SNP, per sample
# Output: 
#   - boxplots of genic SNPs per region_time, and for total SNPs
#----------------------------------------------------------------------

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

work_dir = "../../results/05_SNPs_depths"

#Input files
genic_depths_in = f"{work_dir}/genic_processed_depths.tsv"
intergenic_depths_in = f"{work_dir}/intergenic_processed_depths.tsv"

face_color = "#009688"
median_color = "#795548"
plt.style.use('ggplot')

def plot_depth_boxplots(genic_depths_in):
    years = list(range(2015, 2025))
    reps = ['a', 'b']

    df = pd.read_csv(genic_depths_in, sep="\t")
    sample_cols = df.columns[4:]

    samples_info = []
    for s in sample_cols:
        parts = s.split('_')
        if len(parts) == 3 and len(parts[1]) >= 2:
            region = parts[0]
            season = parts[1][0]
            rep = parts[1][1]
            year = int(parts[2])
            samples_info.append((s, region, season, rep, year))

    groups = {}
    for s, region, season, rep, year in samples_info:
        key = f"{region}_{season}"
        groups.setdefault(key, []).append((s, rep, year))

    for group, samples in groups.items():
        col_labels = [f"{y}{r}" for y in years for r in reps]
        plot_df = pd.DataFrame(index=df.index, columns=col_labels, dtype=float)

        for s, rep, year in samples:
            col_name = f"{year}{rep}"
            plot_df[col_name] = df[s]

        plt.figure(figsize=(16, 6))
        boxplot_data = []
        for col in col_labels:
            data = plot_df[col].dropna()
            boxplot_data.append(data.values if len(data) > 0 else [np.nan])

        x_pos = np.arange(1, len(col_labels) + 1)
        for y in [100, 200, 300]:
            plt.axhline(y=y, color='gray', linestyle='--', linewidth=0.8)

        box = plt.boxplot(
            boxplot_data, positions=x_pos, patch_artist=True,widths=0.6,
            medianprops=dict(color=median_color), showfliers=False
        )        
        # Add median
        for i, median_line in enumerate(box["medians"]):            
            y_data = median_line.get_ydata()
            median_value = y_data[0]
            if pd.notnull(median_value):
            # x_pos[i] is the position of the box
                plt.text(
                    x_pos[i] -0.6, median_value, f"{int(round(median_value))}",  
                    va="center", ha="left", fontsize=10
                )
        # Apply colors
        for patch in box['boxes']:
            patch.set(facecolor=face_color)        
        for median in box['medians']:
            median.set(color=median_color)
        
        plt.ylim(-5, 400)
        plt.xticks(ticks=x_pos, labels=col_labels, rotation=45, fontsize=12)
        plt.title(f"{group}")
        plt.ylabel("Depth")
        plt.tight_layout()
        plt.savefig(f"{work_dir}/boxplots/{group}_depth_boxplots.png")
        plt.close()

def plot_total_depths(genic_depths_in, intergenic_depths_in):
    genic = pd.read_csv(genic_depths_in, sep="\t")
    intergenic = pd.read_csv(intergenic_depths_in, sep="\t")
    data = [genic["TOTAL"], intergenic["TOTAL"]]
    labels = ["Genic", "Intergenic"]
    plt.figure(figsize=(8, 6))
    plt.style.use('ggplot')

    bplot = plt.boxplot(
        data, tick_labels=labels, patch_artist=True, widths=0.2, showfliers=False)

    # Set colors for boxes
    for patch in bplot['boxes']:
        patch.set_facecolor(face_color)
        
    # Set colors for medians
    for median_line in bplot['medians']:
        median_line.set_color(median_color)        

    # Add median value labels
    for i, d in enumerate(data):
        median_val = d.median()
        if pd.notna(median_val):
            median_val_int = int(round(median_val))
            plt.text(
                i+0.83, median_val, f"{median_val_int}",
                ha='center', va='bottom', fontsize=10                
            )
    plt.ylabel("TOTAL depth per SNP")
    plt.tight_layout()
    plt.savefig(f"{work_dir}/boxplots/total_depth_boxplots.png")
    plt.close()

def main(genic_depths_in, intergenic_depths_in):
    boxplot_dir = f"{work_dir}/boxplots"
    if not os.path.exists(boxplot_dir):
        os.makedirs(boxplot_dir)
    plot_depth_boxplots(genic_depths_in)
    plot_total_depths(genic_depths_in, intergenic_depths_in)

if __name__ == "__main__":
    main(genic_depths_in, intergenic_depths_in)
