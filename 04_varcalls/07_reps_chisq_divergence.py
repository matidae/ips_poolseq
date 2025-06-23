#!/usr/bin/env python3

#----------------------------------------------------------------------
# Test SNP reliability based on replicate divergence using allele frequencies as (z-scores).
# Compute chi-square statistic and p-value based on null variance. p-values  < 0.01 are unreliable. 
#
# Inputs:
#   - nv.genic_m_and_z: null variance estimates per sample
#   - paired_samples.txt: sample names and replicate column indices
#   - genic_m_and_z.tsv: allele frequencies and counts per SNP
#
# Output:
#   - snpdev.genic.m_z.txt: SNP-level coverage, comparisons, chi-square stat, and p-value
##----------------------------------------------------------------------


from math import sqrt, asin
from scipy.stats import chi2

wd = "../../results/04_varcalls"
# Input files
null_var_file = f"{wd}/nv.genic_m_and_z"
paired_file = f"{wd}/paired_samples.txt"
genic_m_and_z_file = f"{wd}/genic_m_and_z.tsv" # Genic m and z file

# Ouput file
out_file = f"{wd}/snpdev.genic.m_z.txt"

# Parameters
mincalls = 20  # Minimum number of SNPs per bulk
min2 = 50 # Lowest count to include a sample

def load_null_variance(filename):
    null_var = {}
    with open(filename ,"r") as f:
        for line in f:
            sample, var = line.strip().split('\t')
            null_var[sample] = float(var)
    return null_var

def load_paired_samples(filename):
    xpositions = {}
    with open(filename ,"r") as f:
        for line in f:
            sample, rep_a, rep_b = line.strip().split('\t') 
            xpositions[sample] = [int(rep_a), int(rep_b)]
    return xpositions

def aggregate_SNPs_by_coverage_bin(stats, dcat, tail):
    if dcat not in stats:
        stats[dcat] = [0, 0]
    stats[dcat][0] += 1
    if tail < 0.01:
        stats[dcat][1] += 1    

def write_SNP_stats(outfile, cols_mz, m_total, snpstats, tail, is_valid):
    if is_valid:        
        outfile.write(f"{cols_mz[0]}\t{cols_mz[1]}\t{m_total}\t{snpstats[0]}\t{snpstats[1]:.4f}\t{tail:.4f}\n")
    else:
        outfile.write(f"{cols_mz[0]}\t{cols_mz[1]}\t{snpstats[0]}\t-99\t1.0\n")
       # print(f"full fail {cols_mz[0]}\t{cols_mz[1]}\t{snpstats[0]}\t-99\t1.0")       

null_var = load_null_variance(null_var_file)
paired_samples = load_paired_samples(paired_file)

with open(genic_m_and_z_file, "r") as f_mz, open(out_file, "w") as outfile:
    outfile.write("CHROM\tPOS\tCOV\tNcomp\tChi2\tpval\n")
    next(f_mz)    
    stats={}    

    for line_mz in f_mz:
        cols_mz = line_mz.strip().split('\t')       
        chrom_mz, pos_mz = cols_mz[0], cols_mz[1]

        snpstats=[0,0.0]        
        m_total = 0

        # Iterate over each group of replicates
        for sample in paired_samples:
            rep_a = cols_mz[paired_samples[sample][0]].split(",")
            m1 = int(rep_a[0])
            rep_b = cols_mz[paired_samples[sample][1]].split(",")
            m2 = int(rep_b[0])
            m_total += (m1 + m2)

            #If both replicate pass minimum read depth threshold, then calculate variance
            if m1>=min2 and m2>=min2:
                z1 = float(rep_a[1])
                z2 = float(rep_b[1])
                if (z1 + z2) > 0 and (z1 + z2) < 4.0*asin(1.0): #sum of z values should be between 0 and 2 Pi
                    rx = (z1-z2)/ sqrt(null_var[sample] + 1.0/float(m1) + 1.0/float(m2)) # E[rx] = nullvar
                    snpstats[0]+=1 # Count valid comparisons
                    snpstats[1]+= rx*rx # Accumulate chi-squared-like statistic

        dcat=int(m_total/500)  # Bin SNP by total coverage (e.g., 0,1,... for every 500 reads)
        tail = 1.0
        if snpstats[0]>0:  # If there were valid comparisons
            if snpstats[1]>0:
                tail=1.0-chi2.cdf(snpstats[1], snpstats[0])  # Chi-squared p-value
            else:
                tail=1.0  # No divergence            
            # Write SNP-level statistics
            write_SNP_stats(outfile, cols_mz, m_total, snpstats, tail, is_valid=True)
            # Aggregate statistics by coverage bin
            aggregate_SNPs_by_coverage_bin(stats, dcat, tail)             
        else:  # No valid comparisons
            write_SNP_stats(outfile, cols_mz, m_total, snpstats, tail, is_valid=False)

for dcat in stats:    
    print(f"{dcat*500}\t{stats[dcat][0]}\t{stats[dcat][1]/stats[dcat][0]:.4f}")
