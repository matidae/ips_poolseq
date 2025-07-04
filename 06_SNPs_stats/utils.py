depths_dir = "../../results/05_SNPs_depths"
stats_dir = "../../results/06_SNPs_stats"

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
        return map(float, line.strip().split('\t')[-2:])
