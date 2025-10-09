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

work_dir = "../../results/06_SNPs_stats"

import glob, os
import pandas as pd
from statsmodels.stats.multitest import multipletests

def main():
    for filepath in sorted(glob.glob(f"{work_dir}/tests_*.tsv")):
        prefix = os.path.basename(filepath).split(".")[0].replace("tests_", "")        
        df = pd.read_csv(filepath, sep="\t")
        if not df.empty:        
            # Compute FDR for pval_1df and pval_df
            for col in ["pval_1df", "pval_df"]:
                pvals = df[col].values
                pvals_fdr = multipletests(pvals, method="fdr_bh", alpha=0.05)[1]
                df[f"{col}_FDR"] = pvals_fdr            

            # Save with FDR results
            fdr_out = f"{work_dir}/tests_{prefix}_fdr.tsv"
            df.to_csv(fdr_out, sep="\t", index=False)

if __name__ == "__main__":
    main()