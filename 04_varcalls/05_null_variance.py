#!/usr/bin/env python3

#----------------------------------------------------------------------
# Estimates null variance (v) between pair of replicates  using AF data.
# Removes outliers and low depth SNPs
# Run: ./05_null_variance.py
# In: tables of m,z values calculated by ./04_calculate_m_and_z.py
# Out: paired_samples.txt (paired sample indices)
#      paired_samples.no_replicate.txt (samples without replicates)
#      pairdz.byfreq.txt (histogram of changes by allele frequency)
#      nv.genic_m_and_z (initial nucleotide variation estimate per sample)
#      null_variance_summary.tsv (summary statistics per sample)
#----------------------------------------------------------------------
import csv
from random import shuffle
from math import sqrt, asin, sin

# Files input and outputs
work_dir = "../../results/04_varcalls"
m_and_z_in = f"{work_dir}/genic_m_and_z.tsv"
paired_out = f"{work_dir}/paired_samples.txt"
norep = f"{work_dir}/paired_samples.no_replicate.txt"
out_hist=f"{work_dir}/pairdz.byfreq.txt"
out_nv=f"{work_dir}/nv.genic_m_and_z"
out_summary = f"{work_dir}/null_variance_summary.tsv"

# Generate the file paired_samples.txt
def paired_samples(m_and_z_in, paired_out, norep):
    sample_dict = {}
    xpositions = {} # maps sample to column idx of replicates
    # Open genic_m_and_z to get the header
    with open(m_and_z_in, "r") as m_and_z_fh, open(paired_out, "w") as paired_fh, open(norep, "w") as out_norep:
        reader = csv.reader(m_and_z_fh, delimiter="\t")
        header = next(reader)
        for idx, sample in enumerate(header[4:], start=4):  # start in col 4:
            location, season_rep, year = sample.split("_")
            season = season_rep[0]  # E or L
            rep = season_rep[1]     # a or b
            key = (location, season, year)
            if key not in sample_dict:
                    sample_dict[key] = {}
            sample_dict[key][rep] = idx             
        # Write paired samples to file
        for (location, season, year), reps in sorted(sample_dict.items()):
            if "a" in reps and "b" in reps:
                name = f"{location}_{season}_{year}"
                paired_fh.write(f"{name}\t{reps['a']}\t{reps['b']}\n")
                xpositions[name] = [int(reps['a']), int(reps['b'])]
            else:
                # Write non paired samples to file
                out_norep.write(f"{location}_{season}_{year}: {reps}\n")
    return xpositions

def calculate_null_variance(xpositions, infile):
    min2 = 100 #Minimum read depth
    minMAF = 0.05 #Minimum allele frequency

    #Transformed AF boundaries
    zlow = 2.0 * asin(sqrt(minMAF)) 
    zhigh = 2.0 * asin(sqrt(1.0 - minMAF))

    stats = {sample: [0, 0.0, 0.0] for sample in xpositions}
    changedist = [[0, 0.0] for _ in range(51)] #Bins for divergence magnitude per AF.

    with open(infile, "r") as f:
        for line_idx, line in enumerate(f):
            cols = line.strip().split('\t')     
            if line_idx==0:
                    continue        
            for sample_name, (idx_a, idx_b) in xpositions.items():                             
                count_z_a = cols[idx_a].split(",")
                count_z_b = cols[idx_b].split(",")
                if 'NA' in count_z_a or 'NA' in count_z_b:
                    continue
                count_a = int(count_z_a[0])
                count_b = int(count_z_b[0])

                z_a = float(count_z_a[1])                        
                z_b = float(count_z_b[1])
                
                # Compute mean z-value of the two replicates
                z_bar = (z_a + z_b) / 2.0

                # Skip low coverage pairs
                if count_a < min2 or count_b < min2:
                    continue            
                # Skip pairs outside the threshold
                if not (zlow < z_bar < zhigh):
                    continue

                # Update statistics for this sample:
                stats[sample_name][0] += 1                # Increment count of valid data points
                stats[sample_name][1] += (z_a - z_b) ** 2 # Accumulate sum of squared differences between replicates
                stats[sample_name][2] += (1.0 / count_a) + (1.0 / count_b)  # Accumulate sum of inverse read counts
                # Calculate p-value transformation from z_bar using sine squared function
                p = sin(z_bar / 2.0) ** 2

                # Calculate category index for histogram binning using symmetric binning
                pcat = min(int(p * 100), int((1 - p) * 100))

                # Update histogram bin counters in changedist
                changedist[pcat][0] += 1                   # Increment count of SNPs in this bin
                changedist[pcat][1] += abs(z_a - z_b)     # Accumulate absolute difference in z-scores for this bin
    return stats, changedist

def write_output(stats, changedist):
    with open(out_hist, 'w') as hist_f, open(out_nv, 'w') as nv_f, open(out_summary, 'w') as summary_f:    
        # your final output loops here, replace out1.write with out1_f.write etc.
        for sample_name in xpositions:
            # Calculate mean squared difference of z (mean dz2)
            mean_dz2 = stats[sample_name][1] / float(stats[sample_name][0])
            n = stats[sample_name][0]
            # Calculate mean read depth variance (sum of inverse read counts)
            readdepth_var = stats[sample_name][2] / float(stats[sample_name][0])
            
            rdprop = readdepth_var / mean_dz2
            null_var = mean_dz2 - readdepth_var

            # Print summary statistics to console
            summary_f.write(f"{sample_name}\t{n}\t{mean_dz2:.6f}\t{readdepth_var:.6f}\t{rdprop:.6f}\t{null_var:.6f}\n")            
            
            # Calculate initial nucleotide variation estimate (nv_init)
            nv_init = mean_dz2 - readdepth_var
            
            # Write sample name and nv_init value to output file out2
            nv_f.write(f"{sample_name}\t{nv_init}\n")        

        # Process the changedist histogram bins (only up to 50)
        for pcat in range(51):
            count = changedist[pcat][0]
            total_abs_diff = changedist[pcat][1]    
            # Only process bins with at least one data point
            if count > 0:
                # Calculate average absolute difference in this bin
                avg_abs_diff = total_abs_diff / float(count)        
                # Write pcat, count, and average absolute difference to output file out1
                hist_f.write(f"{pcat}\t{count}\t{avg_abs_diff}\n")

xpositions = paired_samples(m_and_z_in, paired_out, norep)
stats, changedist = calculate_null_variance(xpositions, m_and_z_in)
write_output(stats, changedist)

