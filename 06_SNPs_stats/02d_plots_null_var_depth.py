#!/usr/bin/env python3

#----------------------------------------------------------------------
#  Multiple plots comparing depths, nSNPs, nullvar, n estimated genomes and qubits
#
# Inputs:
#   - null_variance_summary.tsv: null variance estimates per sample
#   - genic_depth_stats.tsv: depths stats for each replicate
#   - sample_qubit: qubits value for each replicate
# Output:
#   - nullvar_nSNPs_plot.png: null variance VS number of SNPs that pass depth filter
#   - ne_genomes_nSNPs_plot.png: n estimated genomes VS nSNPs after filter
#   - depths_nSNPs_plot.png:    average of the median depth of replicates VS nSNPs
#   - ne_genomes_depth_plot.png: n estimated genomes VS average of median depth
#   - depths_qubit_plot_rep_a.png: qubits VS depth of  replicate a
#   - depths_qubit_plot_rep_b.png: qubits VS depth of  replicate b
#----------------------------------------------------------------------

import os
import pandas as pd
import matplotlib.pyplot as plt

work_dir = "../../results/06_SNPs_stats/plots"

# Input files
nullvar_in = "../../results/06_SNPs_stats/null_variance_summary.tsv"
depths_in = "../../results/05_SNPs_depths/genic_depth_stats.tsv"
qubit_in = "../../results/01_proc_reads/sample_qubit"


color= "#009688"
plt.style.use("ggplot")

def plot_nullvar_nSNPs(df_nullvar, out_file):
    plt.figure(figsize=(8,6))
    plt.scatter(df_nullvar['nSNPs'], df_nullvar['nullvar'],
                           alpha=0.7, color=color)
    labels = df_nullvar['sample']    

    for i, label in enumerate(labels):
        plt.text(df_nullvar['nSNPs'][i], df_nullvar['nullvar'][i], label, fontsize=6)
    
    plt.xscale('log')
    plt.xlabel('n_SNPs')
    plt.ylabel('nullvar')
    plt.tight_layout()
    plt.savefig(out_file)
    plt.close()


def plot_ne_genomes_nSNPs(df_nullvar,  out_file):
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
    plt.savefig(out_file)
    plt.close()


def plot_depths_nSNPs(df_nullvar, out_file):
    plt.figure(figsize=(8,6))
    plt.scatter(df_nullvar['nSNPs'], df_nullvar['depth'],
                           alpha=0.7, color=color)
    labels = df_nullvar['sample']    
    for i, label in enumerate(labels):
        plt.text(df_nullvar['nSNPs'][i], df_nullvar['depth'][i], label, fontsize=6)        
    
    plt.ylabel('depth')
    plt.xlabel('number of SNPS')
    plt.tight_layout()
    plt.savefig(out_file)
    plt.close()


def plot_ne_genomes_depths(df_nullvar,  out_file):
    plt.figure(figsize=(8,6))
    plt.scatter(df_nullvar['ne_diploid'], df_nullvar['depth'],
                           alpha=0.7, color=color)
    labels = df_nullvar['sample']    

    for i, label in enumerate(labels):
        plt.text(df_nullvar['ne_diploid'][i], df_nullvar['depth'][i], label, fontsize=6)
    
    plt.xlabel('depth')
    plt.ylabel('effective number of diploid genomes')
    plt.tight_layout()
    plt.savefig(out_file)
    plt.close()


def plot_depths_qubit(df_qubit, out_file, rep):
    df_qubit_rep = df_qubit[df_qubit['sample_rep'].str.contains(rep + "_")]
    df_qubit_rep = df_qubit_rep.reset_index(drop=True)

    plt.figure(figsize=(8,6))    
    plt.scatter(df_qubit_rep['qubit'], df_qubit_rep['depth'],
                           alpha=0.7, color=color)
    labels = df_qubit_rep['sample_rep']
    for i, label in enumerate(labels):
        plt.text(df_qubit_rep['qubit'][i], df_qubit_rep['depth'][i], label, fontsize=6)
           
    plt.ylabel('depth')
    plt.xlabel('qubit')
    plt.tight_layout()
    plt.savefig(out_file)    
    plt.close()


def main():
    # Make directory if it doesn't exist 
    os.makedirs(work_dir, exist_ok=True)

    # Load null variance data
    df_nullvar = pd.read_csv(nullvar_in, sep='\t', usecols=[0, 1, 5])    
    df_nullvar.columns = ['sample', 'nSNPs', 'nullvar']
    # Add number estimated diploid genomes
    df_nullvar["ne_diploid"] = 1 / df_nullvar["nullvar"]

    # Load depth data 
    df_depths = pd.read_csv(depths_in, sep='\t', \
                             skipfooter=1, engine='python' ,usecols=[0, 12])
    df_depths.columns = ['sample_rep', 'depth']
    # Add column with main sample name, removing replicate id
    df_depths["sample"] = df_depths["sample_rep"].str.replace(r"([A-Z])([ab])_", r"\1_", regex=True)
    
    #Add average depth of replicates
    grouped = df_depths.groupby("sample")["depth"]
    df_depths_avg = grouped.mean().reset_index()
    df_nullvar = df_nullvar.merge(df_depths_avg, on="sample", how="left")
    
    # Load qubit data
    df_qubit = pd.read_csv(qubit_in, sep='\t')    
    df_qubit.columns = ['qubit', 'sample_rep']
    
    df_qubit = df_qubit.merge(df_depths, on="sample_rep", how="left")\
        .drop(columns="sample").dropna()
    
    # Call each plot function
    plot_nullvar_nSNPs(df_nullvar, f'{work_dir}/nullvar_nSNPs_plot.png')
    plot_ne_genomes_nSNPs(df_nullvar, f'{work_dir}/ne_genomes_nSNPs_plot.png')
    plot_depths_nSNPs(df_nullvar, f'{work_dir}/depths_nSNPs_plot.png')
    plot_ne_genomes_depths(df_nullvar, f'{work_dir}/ne_genomes_depth_plot.png')
    plot_depths_qubit(df_qubit, f'{work_dir}/depths_qubit_plot_rep_a.png', 'a')
    plot_depths_qubit(df_qubit, f'{work_dir}/depths_qubit_plot_rep_b.png', 'b')


if __name__ == '__main__':
    main()
