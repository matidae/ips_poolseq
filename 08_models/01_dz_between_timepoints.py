#!/usr/bin/env python3

#----------------------------------------------------------------------
# Calculates the values of z for each year_season and the difference (dz) between consecutive periods
#
# Inputs:
#   - genic_m_and_z.filter.tsv: filtered allele frequencies and counts per SNP
# Output:
#   - z_year.{prefix}.tsv — Z values per year_season
#   - dz_interval.{prefix}.tsv — diff in Z values in consecutive year intervals 
#----------------------------------------------------------------------

import sys
sys.path.append("./utils")
import os

work_dir = "../results/08_models"

# Input files

# Output files
#dz_timepoint_out = f"{work_dir}/dz_timepoint.{prefix}.tsv"  # ziff in z values in consecutive intervals 


def get_prefixes(work_dir):
    # look for files like z_mean.<prefix>.tsv
    prefixes = []
    for fname in os.listdir(work_dir):
        if fname.startswith("z_mean.") and fname.endswith(".tsv"):
            # Extract prefix
            prefix = fname.split(".")[1]
            prefixes.append(prefix)
    return sorted(prefixes)

def calculate_dz(prefix):
    z_mean_file = f"{work_dir}/z_mean.{prefix}.tsv"
    dz_file = f"{work_dir}/dz_interval.{prefix}.tsv"

    with open(z_mean_file) as z_mean_fh, open(dz_file, "w") as dz_timepoint_fh:
        header = z_mean_fh.readline().strip().split("\t")
        samples = header[4:]  # skip CHROM,POS,REF,ALT

        # Prepare dz header
        dz_timepoint_fh.write("\t".join(header[:4]))
        for i in range(len(samples)-1):
            s1 = "_".join(samples[i].split("_")[1:])   # 'L_2015'
            s2 = "_".join(samples[i+1].split("_")[1:]) # 'E_2018'
            dz_timepoint_fh.write(f"\t{s1}-{s2}")
        dz_timepoint_fh.write("\n")

        for line in z_mean_fh:
            cols = line.strip().split("\t")
            dz_timepoint_fh.write("\t".join(cols[:4]))
            for i in range(len(samples)-1):
                a = cols[4 + i].split(",")
                b = cols[4 + i + 1].split(",")
                
                if "NA" in (a[0], b[0]):
                    dz_timepoint_fh.write("\tNA,NA")
                    continue
                
                z_a = float(a[0])
                z_b = float(b[0])
                se_a = float(a[1])
                se_b = float(b[1])

                dz = z_b - z_a
                dz_var = se_a**2 + se_b**2
                dz_timepoint_fh.write(f"\t{dz:.4f},{dz_var:.4f}")
            dz_timepoint_fh.write("\n")

def main():
    prefixes = get_prefixes(work_dir)
    for pre in prefixes:
        calculate_dz(pre)

if __name__ == '__main__':
    main()
