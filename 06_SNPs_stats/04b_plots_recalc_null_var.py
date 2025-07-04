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
work_dir = "../../results/06_SNPs_stats"

def plot_histogram_mean_inv_depth(df, output_file='plot1_mean_inv_depth_hist.png'):
    plt.figure(figsize=(8, 5))
    sns.histplot(df['mean_inv_depth'], bins=20, kde=True, color='skyblue')
    plt.title('Distribution of Mean Inverse Depth per Sample')
    plt.xlabel('Mean Inverse Depth')
    plt.ylabel('Number of Samples')
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    print(f'Saved histogram plot to {output_file}')

def plot_boxplot_null_var_by_group(df, output_file='plot2_null_var_boxplot.png'):
    # Extract group prefix from sample name (prefix before first underscore)
    df['group'] = df['sample'].apply(lambda x: x.split('_')[0])
    
    plt.figure(figsize=(10, 6))
    sns.boxplot(x='group', y='null_var', data=df, palette='Set2')
    plt.title('Null Variance Distribution by Sample Group')
    plt.xlabel('Sample Group')
    plt.ylabel('Null Variance')
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    print(f'Saved boxplot to {output_file}')

def plot_scatter_null_var_vs_mean_inv_depth(df, output_file='plot4_null_var_vs_mean_inv_depth.png'):
    plt.figure(figsize=(8, 6))
    sizes = (df['n_snp'] / df['n_snp'].max()) * 300  # scale sizes for visibility
    
    scatter = plt.scatter(df['mean_inv_depth'], df['null_var'],
                          s=sizes, alpha=0.7, c=df['rdprop'], cmap='viridis')
    
    plt.colorbar(scatter, label='Read Depth Proportion (rdprop)')
    plt.title('Null Variance vs Mean Inverse Depth (point size ~ # SNPs)')
    plt.xlabel('Mean Inverse Depth')
    plt.ylabel('Null Variance')
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    print(f'Saved scatter plot to {output_file}')

def main():
    # Load data
    df = pd.read_csv(f'{work_dir}/null_variance_summary.filtered.tsv2', sep='\t')
    
    plot_histogram_mean_inv_depth(df)
    plot_boxplot_null_var_by_group(df)
    plot_scatter_null_var_vs_mean_inv_depth(df)

if __name__ == '__main__':
    main()
