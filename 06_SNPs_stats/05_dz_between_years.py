#!/usr/bin/env python3

#----------------------------------------------------------------------
# Recalculate null variance for remaining SNPs after filtering results from previous steps
#
# Inputs:
#   - null_variance_summary.recalc.tsv: recalculated null variance estimates per sample
#   - genic_m_and_z.filter.tsv: filtered allele frequencies and counts per SNP
# Output:
#   - Z.by.year.tsv — Z values per year
#   - DZ.by.interval.tsv — diff in Z values in consecutive year intervals 
#----------------------------------------------------------------------

import re
from math import sqrt, asin
from utils import load_paired_samples, load_null_variance_recalc


work_dir = "../../results/06_SNPs_stats"

# Input files
null_var_in = f"{work_dir}/null_variance_summary.recalc.tsv"  # Null variance data input
m_and_z_in = f"{work_dir}/genic_m_and_z.filter.tsv"           # Genic m and z data input

# Output files
z_by_year_out = f"{work_dir}/Z.by.year.tsv"            # Z values per year
dz_by_year_out = f"{work_dir}/DZ.by.interval.tsv"      # Diff in Z values in consecutive year intervals 

def calculate_dz_by_year(prefixes, paired_samples, null_var, m_and_z_in):
    minMAF = 0.05
    zlow = 2.0*asin(sqrt(minMAF))
    zhigh= 2.0*asin(sqrt(1.0-minMAF))
    
    with open(m_and_z_in, "r") as m_and_z_fh, open(z_by_year_out, "w") as z_by_year_fh:
        #z_by_year_fh.write("CHROM\tPOS\tREF\tALT")
        next(m_and_z_fh)
        #print(paired_samples)        
        for pre in prefixes:
            paired_samples_select = {}
            # Get all samples from same country and season
            for key, val in paired_samples.items():
                if key.startswith(pre):                
                    paired_samples_select[key] = val
            # 
            # Iterate over each group of replicates
            print(paired_samples_select)
            if len(paired_samples_select) > 1:
                vstat = {}
                for k in paired_samples_select.keys():
                    vstat[k]=[0,0.0,0.0]
                c = 0
                last = []
                z_mean = -1000
                for key, val in paired_samples_select.items():                    
                    rep_a_idx = paired_samples_select[key][0]
                    rep_b_idx = paired_samples_select[key][1]                
                    for line in m_and_z_fh:
                        cols = line.strip().split('\t')                    
                        rep_a = cols[rep_a_idx].split(",")
                        rep_b = cols[rep_b_idx].split(",")
                                                                
                        m_a = int(rep_a[0])                    
                        m_b = int(rep_b[0])
                    
                        if m_a <= 100 and m_b <= 100:
                            z_a = float(rep_a[1])                    
                            z_b = float(rep_b[1])
                            z_mean=(z_a + z_b)/2.0
                            SE = sqrt(null_var[key] + 1.0/float(m_a) + 1.0/float(m_b))/4.0  
                            if z_mean > zlow and z_mean < zhigh:
                                vstat[key][0]+=1
                                vstat[key][1]+=(null_var[key] + 1.0/float(m_a) + 1.0/float(m_b))/4.0
                                vstat[key][2]+=(1.0/float(m_a) + 1.0/float(m_b))/4.0

                        # Get total replicates depth                        
                        if c == 0:
                            last = [ z_mean, (null_var[key] + 1.0/float(m_a) + 1.0/float(m_b))/4.0 ]
                            print(last)
                            c+=1
                        else:
                            if (z_mean+last[0])/2.0 >zlow and (z_mean+last[0])/2.0 <zhigh: 
                                dz2 = z_mean - last[0]
                                evar = (null_var[key] + 1.0/float(m_a) + 1.0/float(m_b))/4.0 + last[1]
        for line in m_and_z_fh:
            cols = line.strip().split('\t')            
            m_total = 0
            # Iterate over each group of replicates
            for sample in paired_samples:
                rep_a = cols[paired_samples[sample][0]].split(",")
                m_a = int(rep_a[0])
                rep_b = cols[paired_samples[sample][1]].split(",")
                m_b = int(rep_b[0])
                # Get total replicates depth
                m_total += (m_a + m_b)
    
def main():
    null_var = load_null_variance_recalc()
    paired_samples = load_paired_samples()
    prefixes = sorted(list({ "_".join(k.split("_")[:-1]) for k in paired_samples.keys() }))
    calculate_dz_by_year(prefixes, paired_samples, null_var, m_and_z_in)

if __name__ == '__main__':
    main()
