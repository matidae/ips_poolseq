#!/usr/bin/env python3

#----------------------------------------------------------------------
# Tests SNP allele frequency trajectories for signatures of selection.
# Uses three models for AF data across years: null(drift), directional, fluctating
# Calculates  likelihood ratios (LRT) and p-values 
#
# Inputs:
#   - z_year.{prefix}.tsv — z-scores per SNP and per year
#
# Output:
#   - tests_{prefix}.tsv — likelihoods, LRTs, and p-values per SNP
#----------------------------------------------------------------------  

import sys
sys.path.append("../utils")
from pathlib import Path
import numpy as np
from math import sqrt, asin, log
from scipy.stats import multivariate_normal, chi2
from utils import load_paired_samples 


work_dir = "../../results/08_evolutionary_dynamics"

# Input files
ne_in = f"{work_dir}/Ne_estimates.tsv" # Years of data and Ne estimates
#z_year_in = f"{work_dir}/z_year.{prefix}.tsv"   # z values per year

# Output files
#tests_out = f"{work_dir}/tests_{prefix}.tsv" # Output file with mu, LRT, p-values for each SNP


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
                ne = float(cols[9])
                return n_years, ne
    return None, None


def run_selection_models(z_year_in, n_years, var_drift, tests_out):    
    with open(tests_out, "w") as tests_fh:
        tests_fh.write("\t".join([
           "CHROM", "POS", "REF", "ALT", "n_years","mu_null", "LL0", "mu_start", "mu_end", "LL1", 
           "LRT1", "pval_LRT1", "LL2", "LRT2", "pval_LRT2" ]) + "\n")
        with open(z_year_in) as z_year_fh:
            cols = z_year_fh.readline().strip().split("\t")
            if len(cols) > 5:
                header = cols[4:]  # All years columns
                header_years = [int(h.split("_")[-1]) for h in header]
                
                for row in z_year_fh:
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
                            # Get the covariance matrix
                            cov_matrix = build_covariance_matrix(years, var_drift, se)
                            # Inverse of covariance matrix
                            inv_cov_matrix = np.linalg.inv(cov_matrix)
                            # Transform z into a numpy array for linalg operations
                            z_array = np.array(z)
                            # Create vector of ones representing the null model of drift
                            null_vector = np.ones(len(z))
                            # Maximum likelihood estimate for the null hypothesis (drift model)
                            mu0 = (null_vector @ inv_cov_matrix @ z_array) / (null_vector @ inv_cov_matrix @ null_vector)
                            LL0 = log_likelihood_null(mu0, z_array, cov_matrix)                            

                            # Design matrix X for linear model (columns: coefficients for mu_start and mu_end)
                            years = np.array(years)
                            total_years = years[-1] - years[0]
                            coef_start = 1.0 - (years - years[0]) / total_years
                            coef_end   = (years - years[0]) / total_years
                            X = np.column_stack((coef_start, coef_end)) #(n, 2) matrix                           
                            
                            # Multiply X transposed by inv_cov_matrix to get the weighted transpose ( X^T * C^-1 )
                            Xt_weighted = X.T @ inv_cov_matrix   #(2, n) matrix
                            
                            # Compute the "normal equation" matrix: (X^T C^-1 X)
                            normal_matrix = Xt_weighted @ X             #(2, 2) matrix
                               
                            # Compute the weighted response vector: (X^T C^-1 z)
                            weighted_response = Xt_weighted @ z_array   #(2, ) matrix                            

                            # Solve for bhat: coefficients (mu_start, mu_end) of the linear model
                            # Equivalent to: bhat = (X^T C^-1 X)^(-1) * (X^T C^-1 z)
                            # bhat = np.linalg.inv(normal_matrix) @ weighted_response
                            bhat = np.linalg.solve(normal_matrix, weighted_response)
                            
                            LL1 = log_likelihood_selection(bhat, z_array, cov_matrix, years)                            

                            # Likelihood ratio test (LRT) for selection
                            LRT_selection = 2.0 * (LL1 - LL0)
                            pvalue_selection = 1.0 - chi2.cdf(LRT_selection, df=1)

                            # Fluctuating (alternative) model: mean vector = observed z
                            # In multivariate normal terms, this is the "perfect fit" log-likelihood.
                            LL2 = multivariate_normal.logpdf(z_array, z_array, cov_matrix)
                            # Likelihood ratio test against null
                            LRT_fluctuating = 2.0 * (LL2 - LL0)
                            #Calculate df as number of extra parameters in fluctuating model compared to null model
                            df = len(z_array) - 1
                            pvalue_fluctuating = 1.0 - chi2.cdf(LRT_fluctuating, df=df)

                            tests_fh.write(
                                "\t".join([*cols[:4], str(n_z), f"{mu0:.6f}", f"{LL0:.6f}", 
                                           f"{bhat[0]:.6f}", f"{bhat[1]:.6f}", 
                                           f"{LL1:.6f}", f"{LRT_selection:.4e}", f"{pvalue_selection:.4e}", 
                                           f"{LL2:.6f}", f"{LRT_fluctuating:.4e}", f"{pvalue_fluctuating:.4e}"]) + "\n")


def build_covariance_matrix(years, var_drift, se):               
    # Get the difference between each year to first year
    t = np.array(years) - years[0]
    # Compute covariance matrix as: cov[i, j] = var_drift * min(t[i], t[j])
    cov_matrix = var_drift * np.minimum.outer(t, t)
    # Add error variance se^2 to diagonal
    np.fill_diagonal(cov_matrix, np.diag(cov_matrix) + np.square(se))
    return cov_matrix


def log_likelihood_null(mu0, z_array, cov_matrix):
    mvec = [mu0] * len(z_array)
    probd = multivariate_normal.logpdf(z_array, mvec, cov_matrix)
    return probd


def log_likelihood_selection(mu_params, z_array, cov_matrix, years):
    mu_start, mu_end = mu_params
    total_years = years.max() - years.min()
    # Expected mean vector: linear interpolation between mu_start and mu_end
    mu_vector = mu_start + (mu_end - mu_start) * (years-years.min())/total_years
    # Compute the log-likelihood using the multivariate normal distribution
    prob_density = multivariate_normal.logpdf(z_array, mean=mu_vector, cov=cov_matrix)
    return prob_density
    

def main():
    paired_samples = load_paired_samples()    
    prefixes = sorted(list({ "_".join(k.split("_")[:-1]) for k in paired_samples.keys()}))    

    # Process for each prefix
    for pre in prefixes:
        n_years, ne = get_years_and_ne(pre, ne_in)
        z_year_in = f"{work_dir}/z_year.{pre}.tsv"
        tests_out = f"{work_dir}/tests_{pre}.tsv"
        if ne is not None:
            var_drift = 1.0 / (2.0 * ne) if ne > 0 else 1.0 / (2.0 * 1.0)
            z_path=Path(z_year_in)
            if z_path.is_file():
                run_selection_models(z_year_in, n_years, var_drift, tests_out)


if __name__ == "__main__":
    main()