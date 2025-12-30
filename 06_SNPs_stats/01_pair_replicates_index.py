#!/usr/bin/env python3

#-------------------------------------------------------------------------------
# Creates a file with the column index of sample replicates
#
# Input: 
#   - genic_m_and_z.tsv: table of m, z values
# Output: 
#   - paired_samples.tsv : indexes of paired samples
#   - paired_samples.no_replicate.tsv : samples without replicates
#-------------------------------------------------------------------------------
import csv

work_dir = "../results/06_SNPs_stats"

# Input files
m_and_z_in = f"{work_dir}/genic_m_and_z.tsv"

# Output files
paired_out = f"{work_dir}/samples_paired.tsv"
not_paired_out = f"{work_dir}/samples_no_replicates.tsv"


# Generate the file paired_samples.txt
def paired_samples(m_and_z_in, paired_out, norep):
    sample_dict = {}
    xpositions = {} # maps sample to column idx of replicates
    # Open genic_m_and_z to get the header
    with open(m_and_z_in, "r") as m_and_z_fh, open(paired_out, "w") as paired_fh, \
        open(not_paired_out, "w") as not_paired_fh:
        reader = csv.reader(m_and_z_fh, delimiter="\t")
        header = next(reader)
        for idx, sample in enumerate(header[4:], start=4):  # start in col 4:
            location, season_rep, year = sample.split("_")
            season = season_rep[0]  # E or L
            rep = season_rep[1]     # a or b
            key = (location, season, year)
            if key not in sample_dict:
                    sample_dict[key] = {}
            sample_dict[key][rep] = idx             
        # Write paired samples to file
        for (location, season, year), reps in sorted(sample_dict.items()):
            if "a" in reps and "b" in reps:
                name = f"{location}_{season}_{year}"
                paired_fh.write(f"{name}\t{reps['a']}\t{reps['b']}\n")
                xpositions[name] = [int(reps['a']), int(reps['b'])]
            else:
                # Write non paired samples to file
                not_paired_fh.write(f"{location}_{season}_{year}: {reps}\n")
    return xpositions

def main():
    paired_samples(m_and_z_in, paired_out, not_paired_out)

if __name__ == "__main__":
    main()