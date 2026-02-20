#!/usr/bin/env python3

#----------------------------------------------------------------------
# Characterizes temporal AFs dynamics of SNPs under a fluctuating model;
# SNPs are minor-allele polarized, and dz is computed  between consecutive intervals. 
# Interval types are inferred from the sign of mean dz  and used to compute a 
# "greenscore" per SNP measuring if it follows the population-avg fluctuation.
#
# Inputs:
#    thinned_file    : tsv file with one SNP per gene 
#    z_file          : tsv file with z-transformed AFs
#
# Outputs:
#   - _deltaZ_per_snp.tsv      : dz per SNP per interval
#   - _mean_deltaZ_ttests.tsv  : mean dz per interval with t-tests
#   - _correlations.tsv        : Pearson correlation matrix between intervals
#   - _greenscores.tsv         : greenscore and rank per SNP
#   - _interval_types.tsv      : interval types (green/yellow) inferred from mean dz
#----------------------------------------------------------------------

import math
import glob
from scipy.stats import pearsonr, ttest_1samp
import csv, os, sys

work_dir = "../results/08_models"
thinned_snps_dir = "../results/08_models/snp_to_gene"
out_dir = "../results/08_models/analysis"

# Input files
thinned_snps_in = sorted(glob.glob(f"{thinned_snps_dir}/tests_*_fluctuating.thinned.tsv"))

# Output files
#tests_fdr_out = [f.replace(".tsv", "_FDR.tsv") for f in tests_in]

prefixes = ["_".join(os.path.basename(f).split("_")[1:3]) for f in thinned_snps_in]

#Read thinned SNP file, returning a set
def read_thinned_file(thinned_file):    
    rsnps = set()
    with open(thinned_file, "r") as f:
        for line in f:
            cols = line.strip().split("\t")            
            rsnps.add(cols[1] + "_" + cols[2]) # chromosome_position as key
    return rsnps


# Reads the Z-transformed allele frequency file 
def read_z_file(z_file, rsnps):
    dznext = {}
    dz_per_snp = {}
    snp_zbars = {}
    with open(z_file, "r") as f:
        header = f.readline().strip().split("\t")
        # Extract timepoint columns automatically
        timepoint_cols = header[4:]  # all columns after CHROM, POS, REF, ALT
        years = [col for col in timepoint_cols]

        # initialize dznext for each interval
        dznext = {year: [] for year in years[:-1]}

        for line in f:
            cols = line.strip().split("\t")
            key = cols[0] + "_" + cols[1]  # chromosome + position

            if key not in rsnps:
                continue

            zlist = []
            valid = 0
            for j in range(4, len(cols)):
                val = cols[j].split(",")[0]
                if val != "NA":
                    zlist.append(float(val))
                    valid += 1

            if valid < len(years):
                continue  # skip SNPs with missing timepoints

            # Minor allele polarization
            zbar = sum(zlist) / len(zlist)
            # Store mean z before polarization
            snp_zbars[key] = zbar
            if zbar > math.pi / 2:
                zlist = [math.pi - z for z in zlist]

            # Compute dz for each interval
            dz_list = []
            for j in range(len(zlist) - 1):
                dz = zlist[j + 1] - zlist[j]
                dznext[years[j]].append(dz)
                dz_list.append(dz)
            dz_per_snp[key] = dz_list

    return dznext, dz_per_snp, years, snp_zbars

# Computes mean dz, t-tests, correlations, and saves output files.
def analyze_and_save(dznext, dz_per_snp, years, prefix, snp_zbars, greenscores, interval_types):
    # - Save dz per SNP per interval
    dz_file = prefix + "_deltaZ_per_snp.tsv"
    with open(dz_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["SNP"] + [f"{years[i]}_to_{years[i+1]}" for i in range(len(years)-1)])
        for snp, dzlist in dz_per_snp.items():
            writer.writerow([snp] + dzlist)

    # - Save mean dz and t-tests
    mean_file = prefix + "_mean_deltaZ_ttests.tsv"
    with open(mean_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["interval", "mean_dz", "t_stat", "p_value"])
        for i, year in enumerate(years[:-1]):
            mean_dz = sum(dznext[year]) / len(dznext[year])
            tval, pval = ttest_1samp(dznext[year], 0.0)
            label = f"{year}_to_{years[i+1]}"
            writer.writerow([label, mean_dz, tval, pval])

    # - Correlation matrix
    corr_file = prefix + "_correlations.tsv"
    with open(corr_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        intervals = [f"{years[i]}_to_{years[i+1]}" for i in range(len(years)-1)]
        writer.writerow([""] + intervals)
        #writer.writerow([""] + years[:-1])
        n = len(years) - 1
        # Pre-compute upper triangle only
        corr_matrix = [[1.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                r, p = pearsonr(dznext[years[i]], dznext[years[j]])
                corr_matrix[i][j] = round(r, 6)
                corr_matrix[j][i] = round(r, 6)  # mirror

        for i, interval in enumerate(intervals):
            writer.writerow([interval] + corr_matrix[i])

    # - Green score per SNP
    greenscore_file = prefix + "_greenscores.tsv"
    with open(greenscore_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["SNP", "greenscore", "mean_z", "rank"])
        sorted_snps = sorted(greenscores.items(), key=lambda x: x[1], reverse=True)
        for rank, (snp, gs) in enumerate(sorted_snps, start=1):
            zbar = snp_zbars[snp]  # mean_z stored during read_z_file
            writer.writerow([snp, round(gs, 6), round(zbar, 6), rank])

    # Correlation greenscore vs mean allele frequency
    gs_vals = [gs for _, gs in sorted_snps]
    zbar_vals = [snp_zbars[snp] for snp, _ in sorted_snps]
    r, p = pearsonr(gs_vals, zbar_vals)

    # Save interval types to file
    interval_file = prefix + "_interval_types.tsv"
    with open(interval_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["interval", "type", "mean_dz"])
        for i, year in enumerate(years[:-1]):
            mean_dz = sum(dznext[year]) / len(dznext[year])
            label = f"{year}_to_{years[i+1]}"
            writer.writerow([label, interval_types[i], round(mean_dz, 6)])

# Infers interval types from mean dz, and scores SNPs by how consistently follow population-avg fluctuation.
def compute_greenscores(dznext, dz_per_snp, years):
    # Infer interval types
    interval_types = []
    for year in years[:-1]:
        mean_dz = sum(dznext[year]) / len(dznext[year])
        interval_types.append("g" if mean_dz > 0 else "y")
    # Score each SNP
    greenscores = {}
    for snp, dz_list in dz_per_snp.items():
        greenscore = 0.0
        for j, dz in enumerate(dz_list):
            greenscore += dz if interval_types[j] == "g" else -dz
        greenscores[snp] = greenscore
    return greenscores, interval_types

def main():
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    for pre in prefixes:
        thinned_file = f"{thinned_snps_dir}/tests_{pre}_FDR_fluctuating.thinned.tsv"
        output_prefix = f"{out_dir}/{pre}_fluctuating"
        z_file = f"{work_dir}/z_year.{pre}.tsv"
        rsnps = read_thinned_file(thinned_file)
        dznext, dz_per_snp, years, snp_zbars = read_z_file(z_file, rsnps)
        greenscores, interval_types = compute_greenscores(dznext, dz_per_snp, years)
        analyze_and_save(dznext, dz_per_snp, years, output_prefix, snp_zbars, greenscores, interval_types)

if __name__ == "__main__":
    main()
