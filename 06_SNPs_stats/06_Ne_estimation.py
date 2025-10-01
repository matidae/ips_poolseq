#!/usr/bin/env python3

#----------------------------------------------------------------------
# Estimates Ne 
#
# Inputs:
#   - Z.by.year.{prefix}.tsv - Z values per year
# Output:
#   - Z.by.year.{prefix}.tsv - Z values per year
#   - dz.max_interval.{prefix}.tsv - diff in Z values in consecutive year intervals 
#----------------------------------------------------------------------

from math import sqrt, asin

work_dir = "../../results/06_SNPs_stats"

# Input files
# Z.by.year.{prefix}.tsv - Z values per year (loaded dynamically)

# Output files
# dz.max_interval.{prefix}.tsv - diff in Z values in largest intervals (written dyamically)

minMAF = 0.05
z_low = 2.0*asin(sqrt(minMAF))
z_high= 2.0*asin(sqrt(1.0-minMAF))

def get_samples_years(prefixes):    
    prefixes_ne = {}
    for pre in prefixes:
        with open(f"{work_dir}/z_by_year.{pre}.tsv") as pre_fh:
            header = pre_fh.readline().strip().split("\t")
            ncols = len(header)            
            if ncols > 5:                
                start_year = header[4].split("_")[-1]
                last_year = header[-1].strip().split("_")[-1]
                prefixes_ne[pre] = [start_year, last_year]    
    return(prefixes_ne)

def calc(prefix):    
    dzlist = []
    meanEvar =0
    meanDz2 = 0    
    z_by_interval_in = f"{work_dir}/z_by_year.{prefix}.tsv"
    with open(z_by_interval_in, 'r') as z_by_interval_fh:
        next(z_by_interval_fh)
        n_row = 0
        for row in z_by_interval_fh:
            cols = row.strip().split("\t")
            nrow += 1
            # Get first and last year z values
            z_first =cols[4].split(",")[0]
            z_last = cols[len(cols)-1].split(",")[0]

            if z_first != "NA"  and z_last != "NA":
                z_first = float(z_first)
                z_last = float(z_last)
                z_mean = (z_first + z_last)/2.0
                if z_mean > z_low and z_mean < z_high:
                    dz = z_first - z_last                    
                    if n_row%2:
                        dz = -dz
                    # Get first and last year SE values and estimate variance
                    se_first = float(cols[4].split(",")[1])
                    se_last = float(cols[-1].split(",")[1])
                    evar = se_first**2 + se_last **2
                    meanEvar += evar
                    meanDz2 += (dz**2)
                    dzlist.append([dz, evar])
    return(dzlist, meanDz2, meanEvar)


def main():    
    paired_samples = load_paired_samples()
    prefixes = sorted(list({ "_".join(k.split("_")[:-1]) for k in paired_samples.keys() }))
    samples = get_samples_years(prefixes)    
    for key, value in samples.items():
        prefix = key 
        dzlist, meanDz2, meanEvar = calc(prefix)

        nx = len(dzlist)
        if nx > 0:
            dzlist.sort()
            qlist= {}
            for i in range(1, 20):
                qlist[i] = dzlist[int(i*nx/20)][0]
            
            IQR = qlist[15] - qlist[5]
            Q25 = -IQR/2
            start_year = value[0]
            last_year = value[1]
            t = int(value[1]) - int(value[0])
            print(prefix, start_year, last_year)
            print("IQR=",IQR,"Q25=", Q25)
            print("SNPs included ",nx)
            meanEvar=meanEvar/float(nx) 
            meanDz2=meanDz2/float(nx) 
            print("cases",nx,"mean dz2, Estimation var",meanDz2,meanEvar)
            x1 = (IQR/1.34896)**2.0 - meanEvar 
            print("IQR Ne",t/(2*x1))
            x1 = meanDz2 - meanEvar 
            print("Ez2 Ne",t/(2*x1))
        else:
            print(f"No SNPs passed filters for {prefix}")
        
        print("-"*30)


if __name__ == '__main__':
    main()