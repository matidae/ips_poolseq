#!/usr/bin/env python3

#-------------------------------------------------------------------------------
# Calculates null variance from transformed allele frequencies (z-scores) in paired samples (reps).
# Filters out low-quality SNPs, removes outliers and low depth SNPs.
#
# Input: tables of m,z values calculated by ./04_calculate_m_and_z.py
# Output: 
#   - dz2_by_pbin.tsv : histogram of mean reps dz2 changes by allele frequencies bins
#   - null_variance_summary.tsv : summary statistics per sample 
#    (mean dz2 between reps, mean read depth var, prop dz2 explained by read depth, null variance)
#-------------------------------------------------------------------------------

from math import sqrt, asin, sin
from utils import load_paired_samples

work_dir = "../../results/06_SNPs_stats"

# Input files
m_and_z_in = f"{work_dir}/genic_m_and_z.tsv"

# Output files
dz2_bin_out=f"{work_dir}/dz2_by_pbin.tsv"
summary_out = f"{work_dir}/null_variance_summary.tsv"

def calculate_null_variance(samples_reps, m_and_z_in):
    min_read_depth = 100 # Minimum read depth
    min_MAF = 0.05 # Minimum minor allele frequency
    # Transformed AF boundaries
    zlow = 2.0 * asin(sqrt(min_MAF)) 
    zhigh = 2.0 * asin(sqrt(1.0 - min_MAF))

    stats = {sample: [0, 0.0, 0.0] for sample in samples_reps} # Statistics for samples
    changedist = [[0, 0.0] for _ in range(51)] # Bins for dz2 changes by AF.

    with open(m_and_z_in, "r", newline="") as m_and_z_fh:
        next(m_and_z_fh)
        for line in m_and_z_fh:            
            cols = line.strip().split('\t')
            for sample_name, (idx_a, idx_b) in samples_reps.items():
                m_z_a = cols[idx_a].split(",")
                m_z_b = cols[idx_b].split(",")                
                if 'NA' in m_z_a or 'NA' in m_z_b:
                    continue
                # Load variables m (depth) and z (transformed AF)
                m_a, z_a = int(m_z_a[0]), float(m_z_a[1])
                m_b, z_b = int(m_z_b[0]), float(m_z_b[1])
                 # Compute mean z-value of the two replicates
                z_mean = (z_a + z_b) / 2.0
                # Skip low coverage pairs
                if m_a < min_read_depth or m_b < min_read_depth:                    
                    continue               
                # Skip pairs outside the threshold
                if not (zlow < z_mean < zhigh):                    
                    continue
                
                # Increment count of valid data points
                stats[sample_name][0] += 1
                # Accumulate sum of squared differences between reps
                stats[sample_name][1] += (z_a - z_b) ** 2 
                # Accumulate sum of inverse read counts
                stats[sample_name][2] += (1.0 / m_a) + (1.0 / m_b)  
                
                # Calculate p transformation from z_mean using sine squared function
                p = sin(z_mean / 2.0) ** 2
                # Calculate category index for histogram binning using symmetric binning
                pcat = min(int(p * 100), int((1 - p) * 100))                
                # Increment count of SNPs in this bin
                changedist[pcat][0] += 1
                # Accumulate absolute difference in z for this bin
                changedist[pcat][1] += abs(z_a - z_b)

    return stats, changedist

def write_output(stats, changedist, samples_reps, dz2_bin_out, summary_out):
    with open(dz2_bin_out, 'w') as dz2_bin_fh, open(summary_out, 'w') as summary_fh:
        summary_fh.write("Sample\tCount\tDZ2_mean\tDepth_var\tDepth_prop_var\tNull_var\n")
        for sample_name in samples_reps:            
            n = int(stats[sample_name][0])
            # Calculate mean of squared difference of z (dz2_ mean)
            dz2_mean = stats[sample_name][1] / n            
            # Calculate mean read depth variance (sum of inverse read counts)
            readdepth_var = stats[sample_name][2] / n
            # Proportion of total variance explained by read depth variation
            rdprop = readdepth_var / dz2_mean
            # Calculate initial nucleotide variation estimate (nv_init)
            null_var = dz2_mean - readdepth_var
            
            summary_fh.write(f"{sample_name}\t{n}\t{dz2_mean:.6f}\t{readdepth_var:.6f}\t{rdprop:.6f}\t{null_var:.6f}\n")
            
        dz2_bin_fh.write("pbin\tcount\tavg_abs_diff\n")
        # Process the changedist histogram symmetrical bins
        for pbin in range(51):
            count = changedist[pbin][0]
            total_abs_diff = changedist[pbin][1]    
            # Only process bins with at least one data point
            if count > 0:
                # Calculate average absolute difference in this bin
                avg_abs_diff = total_abs_diff / float(count)                        
                # Write pcat, count, and average absolute difference to output file out1
                dz2_bin_fh.write(f"{pbin}\t{count}\t{avg_abs_diff:.6f}\n")       

def main():
    samples_reps = load_paired_samples()
    stats, changedist = calculate_null_variance(samples_reps, m_and_z_in)
    write_output(stats, changedist, samples_reps, dz2_bin_out, summary_out)

if __name__ == "__main__":
    main()
