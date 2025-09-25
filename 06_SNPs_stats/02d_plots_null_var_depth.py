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
work_dir_05 = "../../results/05_SNPs_depths"
work_dir_01 = "../../results/01_proc_reads"

color= "#009688"
plt.style.use("ggplot")

def plot_nullvar_nSNPs(df_nullvar, output_file='nullvar_nSNPs_plot.png'):
    plt.figure(figsize=(8,6))
    plt.scatter(df_nullvar['nSNPs'], df_nullvar['nullvar'],
                           alpha=0.7, color=color)
    labels = df_nullvar['sample']    

    for i, label in enumerate(labels):
        plt.text(df_nullvar['nSNPs'][i], df_nullvar['nullvar'][i], label, fontsize=6)
    
    plt.xlabel('n_SNPs')
    plt.ylabel('nullvar')
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def plot_ne_genomes_nSNPs(df_nullvar,  output_file='nSNPs_ne_genomes_plot.png'):
    plt.figure(figsize=(8,6))
    plt.scatter(df_nullvar['nSNPs'], df_nullvar['ne_diploid'],
                           alpha=0.7, color=color)
    labels = df_nullvar['sample']    

    for i, label in enumerate(labels):
        plt.text(df_nullvar['nSNPs'][i], df_nullvar['ne_diploid'][i], label, fontsize=6)
    
    plt.xscale('log')
    plt.xlabel('number of SNPs')
    plt.ylabel('effective number of diploid genomes')
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()


def plot_depths_nSNPs(df_nullvar, output_file='depths_nSNPs_plot.png'):
    plt.figure(figsize=(8,6))
    plt.scatter(df_nullvar['nSNPs'], df_nullvar['depth'],
                           alpha=0.7, color=color)
    labels = df_nullvar['sample']    
    for i, label in enumerate(labels):
        plt.text(df_nullvar['nSNPs'][i], df_nullvar['depth'][i], label, fontsize=6)        
    
    plt.ylabel('median depth')
    plt.xlabel('number of SNPS')
    plt.savefig(output_file)
    plt.close()

def plot_ne_genomes_depths(df_nullvar,  output_file='ne_genomes_depth_plot.png'):
    plt.figure(figsize=(8,6))
    plt.scatter(df_nullvar['ne_diploid'], df_nullvar['depth'],
                           alpha=0.7, color=color)
    labels = df_nullvar['sample']    

    for i, label in enumerate(labels):
        plt.text(df_nullvar['ne_diploid'][i], df_nullvar['depth'][i], label, fontsize=6)
    
    plt.xscale('log')
    plt.xlabel('depth')
    plt.ylabel('effective number of diploid genomes')
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()


def plot_depths_qubit(df_qubit, output_file='depths_qubit_plot.png'):
    plt.figure(figsize=(8,6))        
    
    plt.scatter(df_qubit['qubit'], df_qubit['depth'],
                           alpha=0.7)
    labels = df_qubit['sample_rep']
    for i, label in enumerate(labels):
        plt.text(df_qubit['qubit'][i], df_qubit['depth'][i], label, fontsize=6)
           
    plt.ylabel(' depth')
    plt.xlabel('qubit')
    plt.savefig(output_file)
    plt.close()


def main():
    # Load null variance data
    df_nullvar = pd.read_csv(f'{work_dir}/null_variance_summary.tsv', sep='\t', \
                             skiprows=1, usecols=[0, 1, 5])    
    df_nullvar.columns = ['sample', 'nSNPs', 'nullvar']
    # Add number estimated diploid genomes
    df_nullvar["ne_diploid"] = 1 / df_nullvar["nullvar"]    

    # Load depth data 
    df_depths = pd.read_csv(f'{work_dir_05}/genic_depth_stats.tsv', sep='\t', \
                            skiprows=1, skipfooter=1, engine='python' ,usecols=[0, 12])
    df_depths.columns = ['sample_rep', 'depth']
    # Add column with main sample name, removing replicate id
    df_depths["sample"] = df_depths["sample_rep"].str.replace(r"([A-Z])([ab])_", r"\1_", regex=True)
    
    #Add average depth of replicates
    grouped = df_depths.groupby("sample")["depth"]
    df_depths_avg = grouped.mean().reset_index()
    df_nullvar = df_nullvar.merge(df_depths_avg, on="sample", how="left")
    
    # Load qubit data
    df_qubit = pd.read_csv(f'{work_dir_01}/sample_qubit', sep='\t', usecols=[0, 1])
    df_qubit.columns = ['qubit', 'sample_rep']
    df_qubit = df_qubit.merge(df_depths, on="sample_rep", how="left").drop(columns="sample") 
        
    plot_nullvar_nSNPs(df_nullvar)
    plot_ne_genomes_nSNPs(df_nullvar)
    plot_depths_nSNPs(df_nullvar)
    plot_ne_genomes_depths(df_nullvar)
    plot_depths_qubit(df_qubit)


if __name__ == '__main__':
    main()
