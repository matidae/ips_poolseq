#!/usr/bin/env python3

#----------------------------------------------------------------------
# 
#
# Inputs:
#   - 
# Output:
#   - 
#----------------------------------------------------------------------

from pathlib import Path
import numpy as np
from math import sqrt, asin, log
from scipy.stats import multivariate_normal, chi2
from utils import load_paired_samples 


work_dir = "../../results/06_SNPs_stats"

# Input files
ne_in = f"{work_dir}/Ne_estimates.tsv"
#z_by_year_in = f"{work_dir}/Z.by.year.{prefix}.tsv"   # Z values per year

# Output files


minMAF = 0.05
z_low = 2.0 * asin(sqrt(minMAF))
z_high = 2.0 * asin(sqrt(1.0 - minMAF))


def get_years_and_ne(prefix, ne_in):
    with open(ne_in) as ne_fh:
        next(ne_fh)
        for line in ne_fh:
            cols = line.strip().split("\t")
            if cols[0] == prefix:
                n_years = int(cols[3])
                ne      = float(cols[9])
                return (n_years, ne)
    return (None, None)


def tests(z_by_year_in, n_years, var_drift):
    
#    with open(ne_tests_out, "w") as ne_tests_fh:        
#       ne_tests_fh.write("\t".join([
#           "CHROM", "POS", "REF", "ALT", "tests","mu_null", "LL0", "mu98", "mu21", "LL1", 
#           "LRT_1df", "pval_1df", "LL_sat", "LRT_df", "pval_df" ]) + "\n")
        with open(z_by_year_in) as z_by_year_fh:
            cols = z_by_year_fh.readline().strip().split("\t")
            if len(cols) > 5:
                header = cols[4:]  # All years columns
                header_years = [int(h.split("_")[-1]) for h in header]
                
                for row in z_by_year_fh:
                    cols = row.rstrip().split("\t")
                    z = []
                    se = []
                    years = []
                    
                    for idx, value in enumerate(cols[4:]):
                        if not "NA" in value.split(",")[0]:
                            z.append(float(value.split(",")[0]))
                            se.append(float(value.split(",")[1]))
                            years.append(header_years[idx])                  
                    n_z = len(z)                    
                    if n_z >= n_years:
                        z_mean = sum(z)/n_z
                        if z_mean > z_low and z_mean < z_high:
                            # Drift model (null hypothesis)
                            # Initialize an empty covariance matrix with zeros
                            cov_matrix = np.zeros((n_z, n_z))                            
                            # Get the difference between each year to first year
                            t = np.array(years) - years[0]
                            # Compute covariance matrix as: cov[i, j] = var_drift * min(t[i], t[j])
                            cov_matrix = var_drift * np.minimum.outer(t, t)
                            # Add error variance se^2 to diagonal
                            np.fill_diagonal(cov_matrix, np.diag(cov_matrix) + np.square(se))

                            # Inverse of covariance matrix
                            inv_cov_matrix = np.linalg.inv(cov_matrix)
                            # Transform z into a numpy array for linalg operations
                            z_array = np.array(z)
                            # Create vector of ones representing the null model of drift
                            null_vector = np.ones(len(z))
                            # Maximum likelihood estimate for the null hypothesis (drift model)
                            mu0 = (null_vector @ inv_cov_matrix @ z_array) / (null_vector @ inv_cov_matrix @ null_vector)
                            LL0 = log_likelihood_drift(mu0, z_array, cov_matrix)
                            LL01 = LL_b(mu0, z_array, cov_matrix)
                            print("LL0:  ", LL0)

                            

def log_likelihood_drift(mu0, z_array, cov_matrix):
    mvec = [mu0] * len(z_array)
    probd = multivariate_normal.pdf(z_array, mvec, cov_matrix)
    if probd > 0 :
        return log(probd)
    else:
        return -99999

def log_likelihood_selection(mu_params, z_array, cov_matrix, years):
    mu_start, mu_end = mu_params
    total_years = years.max() - years.min()
    # Expected mean vector: linear interpolation between mu_start and mu_end
    mu_vector = mu_start + (mu_end - mu_start) * (years-years.min())/total_years
    # Compute the log-likelihood using the multivariate normal distribution
    prob_density = multivariate_normal.pdf(z_array, mean=mu_vector, cov=cov_matrix)
    
    if prob_density > 0:
        return log(prob_density)
    else:
        return -99999



def main():
    paired_samples = load_paired_samples()    
    prefixes = sorted(list({ "_".join(k.split("_")[:-1]) for k in paired_samples.keys()}))    

    # Process for each prefix
    for pre in prefixes:
        n_years, ne = get_years_and_ne(pre, ne_in)
        z_by_year_in = f"{work_dir}/z_by_year.{pre}.tsv"
       # ne_tests_out = f"{work_dir}/Ne_tests_{pre}.tsv"
        if ne is not None:
            var_drift = 1.0 / (2.0 * ne) if ne > 0 else 1.0 / (2.0 * 1.0)
            z_path=Path(z_by_year_in)
            if z_path.is_file():
                tests(z_by_year_in, n_years, var_drift)

if __name__ == "__main__":
    main()