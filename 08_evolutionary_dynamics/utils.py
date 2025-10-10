import csv

depths_dir = "../../results/05_SNPs_depths"
stats_dir = "../../results/06_SNPs_stats"
null_var_dir = "../../results/07_null_variance"

def parse_counts(ad_field):    
    ref, alt = map(int, ad_field.split(','))
    return ref, alt

def load_paired_samples():
    samples_reps = {}
    with open(f"{stats_dir}/samples_paired.tsv") as fh:
        for line in fh:
            sample, idx1, idx2 = line.strip().split('\t')
            samples_reps[sample] = [int(idx1), int(idx2)]
    return samples_reps

def load_depth_threshold():
    depth_threshold = f"{depths_dir}/genic_depth_stats.tsv"
    with open(depth_threshold ,"r") as f:
        for line in f:
            pass
        min_depth, max_depth = map(float, line.strip().split('\t')[-2:])
        return min_depth, max_depth

def load_null_variance_recalc():
    null_var = {}
    null_var_file = f"{null_var_dir}/null_variance_summary.recalc.tsv"    
    with open(null_var_file, newline='') as null_var_fh:
        next(null_var_fh)
        null_var_read = csv.reader(null_var_fh, delimiter='\t')
        for row in null_var_read:
            null_var[row[0]] = float(row[1])
    return null_var
