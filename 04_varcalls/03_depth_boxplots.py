#!/usr/bin/env python3

#----------------------------------------------------------------------
# Makes boxplots grouped by country, region and season.
# Run: ./03_depth_boxplots.py
# In: tables of statistics calculated by ./02_get_SNPs_depth.py
# Out: boxplots of genic, intergenic, and total sum of SNPs
#----------------------------------------------------------------------

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

wd = "../../results/04_varcalls"
face_color = "#009688"
edge_color = "#3C3C3C"
median_color = "#795548"

def make_boxplots(stats_file):
    # Load the stats file
    df = pd.read_csv(stats_file, sep="\t", index_col=0)
    df.index = df.index.str.strip() 
    df = df.iloc[:-1, :]    
    # List samples and prefixes for grouping    
    samples = [s.strip() for s in df.index.tolist()]
    group_prefixes = sorted(set([i[:6] for i in samples]))
    grouped_samples = [[s for s in samples if s.startswith(prefix)] for prefix in group_prefixes]
    # Precompute boxplot data for all samples
    boxplot_data_all = {}
    for sample in samples:
        row = df.loc[sample]
        # Store data in dict for quick lookup
        boxplot_data_all[sample] = {
            'med': row["50%"],
            'q1': row["25%"],
            'q3': row["75%"],
            'whislo': row["whislo"],
            'whishi': row["whishi"],
            'fliers': [],
            'label': sample
        }

    # Define fixed order of x-axis labels
    years = list(range(2015, 2025))
    reps = ['a', 'b']
    fixed_labels = [f"{year}_{rep}" for year in years for rep in reps]     
    
    for prefix, samples_group in zip(group_prefixes, grouped_samples):        
        label_to_sample = {}
        actual_labels = []
        # Make new labels to sort and group properly
        for sample in samples_group:            
            parts = sample.split("_")  # ['SFIN', 'Ea', '2015']
            rep = parts[1][-1]         # 'a'
            year = parts[2]            # '2015'
            simple_label = f"{year}_{rep}"  # '2015_a'
            label_to_sample[simple_label] = sample
        # Build boxplot data in fixed order, insert null values if data for year is missing
        boxplot_data = []        
        for label in fixed_labels:
            
            if label in label_to_sample:
                sample = label_to_sample[label]
                boxplot_data.append(boxplot_data_all[sample])
            else:                
                boxplot_data.append({
                    'med': None, 'q1': None, 'q3': None, 'whislo': None, 'whishi': None, 'fliers': [], 'label': label                    
                })
            actual_labels.append(label)

        # Create a new figure for each group
        fig, ax = plt.subplots(figsize=(20, 6))
        ax.axhline(100, color='gray', linestyle='--')
        ax.axhline(200, color='gray', linestyle='--')
        ax.axhline(300, color='gray', linestyle='--')
        plt.style.use("ggplot")
        # Draw boxplot for the group        
        positions = list(range(1, len(boxplot_data) + 1))
        bp = ax.bxp(boxplot_data, positions=positions, showmeans=False, showfliers=False, manage_ticks=True, patch_artist=True)
        ax.set_xticks(positions)
        # Set new styles for the boxplots        
        for box in bp['boxes']:
            box.set(facecolor=face_color, alpha=0.7)
            box.set(edgecolor=edge_color, linewidth=1.5) 
        for median in bp['medians']:
            median.set(color=median_color, linewidth=2) 
        for whisker in bp['whiskers']:
            whisker.set(color=edge_color, linewidth=1.5)
        for cap in bp['caps']:
            cap.set(color=edge_color, linewidth=1.5)
        for pos, box in zip(positions, boxplot_data):
            if box['med'] is not None:
                ax.text(pos - 0.28, box['med'] -5 , f"{box['med']:.0f}", ha='right', va='bottom', fontsize=14, color='black')
        ax.set_xticklabels(actual_labels, rotation=60, fontsize=14)        
        ax.set_title(prefix)        
        ax.set_ylabel("Depth", fontsize=14)        
        plt.ylim(-10, 400)
        plt.tight_layout()        
        plot_name = prefix + ".SNPs_depth.png"       
        plt.savefig(f"{wd}/{plot_name}", dpi=300, bbox_inches='tight')
        plt.close()

def make_total_boxplots(genic_stats, intergenic_stats):
    genic = pd.read_csv(genic_stats, sep="\t", header=0).tail(1)
    intergenic = pd.read_csv(intergenic_stats, sep="\t", header=0).tail(1)
    combined = pd.concat([genic, intergenic])
    combined.index = ['TOTAL_genic', 'TOTAL_intergenic']

    box_data = []
    for idx, row in combined.iterrows():
        box_data.append({
            'label': idx,
            'med': row['50%'],
            'q1': row['25%'],
            'q3': row['75%'],
            'whislo': row['whislo'],
            'whishi': row['whishi'],
            'fliers': []  
        })
    # Create the plot
    fig, ax = plt.subplots(figsize=(6, 5))
    bp = ax.bxp(box_data, showfliers=False, patch_artist=True)
    ax.set_ylabel("Depth")
    ax.set_title("Total SNPs depth (all samples combined)")
    # Set new styles for the boxplots
    positions = list(range(1, len(box_data) + 1))
    for box in bp['boxes']:
        box.set(facecolor=face_color, alpha=0.7)
        box.set(edgecolor=edge_color, linewidth=1.5) 
    for median in bp['medians']:
        median.set(color=median_color, linewidth=2) 
    for whisker in bp['whiskers']:
        whisker.set(color=edge_color, linewidth=1.5)
    for cap in bp['caps']:
        cap.set(color=edge_color, linewidth=1.5)
    for pos, box in zip(positions, box_data):
        if box['med'] is not None:
            ax.text(pos - 0.1, box['med'] -5 , f"{box['med']:.0f}", ha='right', va='bottom', fontsize=14, color='black')
    plt.style.use("ggplot")
    plt.tight_layout()    
    plt.savefig(f"{wd}/total_SNPs_depth.png", dpi=300, bbox_inches='tight')
    plt.close()
    
make_boxplots(f"{wd}/genic.snp_total_depths.stats.tsv")
make_boxplots(f"{wd}/intergenic.snp_total_depths.stats.tsv")
make_total_boxplots(f"{wd}/genic.snp_total_depths.stats.tsv", f"{wd}/intergenic.snp_total_depths.stats.tsv")
