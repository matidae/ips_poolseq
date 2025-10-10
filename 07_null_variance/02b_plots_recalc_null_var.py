#!/usr/bin/env python3

#----------------------------------------------------------------------
# Recalculate null variance for remaining SNPs after filtering results from previous steps
#
# Inputs:
#   - null_variance_summary.tsv: null variance estimates per sample
#   - genic_m_and_z.tsv: allele frequencies and counts per SNP
#   - snpdev_m_and_z.tsv: SNP-level coverage, comparisons, chi-square stat, and p-value
# Output:
#   - genic_m_and_z.filter.tsv — Filtered SNP data, only for SNPs passing criteria.
#   - null_variance_summary.filter.tsv — Estimated null variance per sample group
#----------------------------------------------------------------------

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

work_dir = "../../results/07_null_variance/"
plot_dir = "../../results/07_null_variance/plots"


# Input files
null_var_recalc_in = f'{work_dir}/null_variance_summary.recalc.tsv'
# Output files
null_var_boxplot = f"{plot_dir}/02b_null_var_boxplot.png"


def plot_boxplot_null_var_by_group(df, plot_file):
    # Extract group prefix from sample name (prefix before first underscore)
    df['group'] = df['Sample'].apply(lambda x: x.split('_')[0])
    
    plt.figure(figsize=(10, 6))
    sns.boxplot(x='group', y='Null_var', data=df, palette='Set2', legend=False, hue='group')
    plt.title('Null variance distribution by sample group')
    plt.xlabel('Sample group')
    plt.ylabel('Null variance')
    plt.tight_layout()
    plt.savefig(plot_file)
    plt.close()

def main():
    # Load data
    df = pd.read_csv(null_var_recalc_in, sep='\t')
    plot_boxplot_null_var_by_group(df, null_var_boxplot)

if __name__ == '__main__':
    main()
