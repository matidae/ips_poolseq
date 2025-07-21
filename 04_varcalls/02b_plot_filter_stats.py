#!/usr/bin/env python3

#----------------------------------------------------------------------
# Makes barplots of counts of missing GT per sample (labels those with more 1% SNPs missing GT)
# Input: 
#   - missing_GT_stats.tsv : table with counts and proportion of missing GT per sample
# Output: 
#   - barplot of percent of missing SNPs missing GT per sample
#----------------------------------------------------------------------

import matplotlib.pyplot as plt

work_dir = "../../results/04_varcalls"
# Input files
gt_in = f"{work_dir}/missing_GT_stats.tsv" 
# Output files
gt_barplot = f"{work_dir}/missing_GT_barplot.png" 

def main(gt_in, gt_barplot):
    # Load data
    samples = []
    proportions = []
    counts = []

    with open(gt_in) as gt_fh:
        for line in gt_fh:
            parts = line.strip().split()
            pct_missing = float(parts[2]) * 100        
            samples.append(parts[0])
            counts.append(int(parts[1]))
            proportions.append(pct_missing)

    # Barplot
    plt.style.use('ggplot')
    plt.figure(figsize=(12, 6))
    bars = plt.bar(samples, proportions, color="#009688")
    plt.xticks(rotation=90, fontsize=8)
    plt.ylabel("Proportion of SNPs missing genotype per sample")
    plt.ylim(0, max(proportions) + 5)

    # Add labels to bars with pct of missing GT > 1%
    for i, bar in enumerate(bars):
        if bar.get_height() >= 1: 
            height = bar.get_height()
            label = f"{counts[i]}%\n({proportions[i]:.1f})"
            plt.text(bar.get_x() + bar.get_width() / 2, height + 0.1, label,
                        ha='center', va='bottom', fontsize=8, rotation=0, color="black")
    # Save plot
    plt.tight_layout()
    plt.savefig(gt_barplot, dpi=300)
    plt.close()

if __name__ == "__main__":
    main(gt_in, gt_barplot)