#!/usr/bin/env python3

#----------------------------------------------------------------------
# Adjusts p-values from likelihood ratio tests for multiple comparisons
# using the Benjamini-Hochberg FDR correction.
#
# Inputs:
#   - tests_{prefix}.tsv — evolutionary models tests results per SNP
#
# Output:
#   - tests_{prefix}_fdr.tsv — add to input file FDR results
#----------------------------------------------------------------------

import glob, os
import pandas as pd
from statsmodels.stats.multitest import multipletests

work_dir = "../../results/06_SNPs_stats"

# Input files
tests_in = sorted(glob.glob(f"{work_dir}/tests_*.tsv"))

# Output files
tests_fdr_out = [f.replace(".tsv", "_FDR.tsv") for f in tests_in]

def compute_fdr(df, columns, alpha=0.05):
    # Apply Benjamini–Hochberg FDR correction to selected p-value columns
    for col in columns:        
        pvals_fdr = multipletests(df[col].values, method="fdr_bh", alpha=alpha)[1]
        df[f"{col}_FDR"] = pvals_fdr
    return df

def main():
    # Process each input/output file pair
    for infile, outfile in zip(tests_in, tests_fdr_out):
        df = pd.read_csv(infile, sep="\t")        
        if df.empty:
            continue
        # Compute FDR for pval_1df and pval_df
        df = compute_fdr(df, ["pval_1df", "pval_df"])
        # Save with FDR results
        df.to_csv(outfile, sep="\t", index=False)

if __name__ == "__main__":
    main()