#!/usr/bin/env python3

import re
import csv
from math import sqrt, asin
from utils import load_paired_samples


work_dir = "../../results/06_SNPs_stats"

# Input files
null_var_in = f"{work_dir}/null_variance_summary.recalc.tsv"  # Null variance data input
m_and_z_in = f"{work_dir}/genic_m_and_z.filter.tsv"           # Genic m and z data input

# Output files
z_by_year_out = f"{work_dir}/Z.by.year.tsv"            # Filtered SNP data output
dz_by_year_out = f"{work_dir}/DZ.by.interval.tsv"      # Recalculated null variance output


def parse_groupings_from_genic_file(m_and_z_in): 
    # Read header only
    with open(m_and_z_in, 'r') as m_and_z_fh:
        header =  m_and_z_fh.readline().strip().split('\t')
    
    # Skip the first 4 columns: CHROM, POS, REF, ALT
    sample_columns = header[4:]

    grouped = {}

    for col in sample_columns:
        # Match pattern like: LAAU_Ea_2023 or SFIN_Lb_2018
        match = re.match(r"([A-Z]+)_(E|L)(a|b)_(\d{4})", col)
        group, season, rep, year = match.groups()

        if group not in grouped:
            grouped[group] = {}
        if season not in grouped[group]:
            grouped[group][season] = {}
        if rep not in grouped[group][season]:
            grouped[group][season][rep] = []

        grouped[group][season][rep].append(col)
        print(grouped)

    return grouped

def load_null_variances(null_var_file):
    null_vars = {}
    with open(null_var_file) as null_var_fh:
        reader = csv.DictReader(null_var_fh, delimiter="\t")
        for row in reader:
            null_vars[row["sample"]] = float(row["null_var"])
    return null_vars

    
def main():
    grouped = parse_groupings_from_genic_file(f'{m_and_z_in}')
    print(grouped['SFIN']['E'])

if __name__ == '__main__':
    main()
