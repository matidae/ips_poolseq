#!/usr/bin/env python3

"""
Script that  generates a YAML metadata file for PoolSeq samples. 

Each sample entry in the YAML file includes:
- id, country, region, season, year, replicate, list of fastq as R1, R2 pairs

Input:
- prefix_file: file containing sample prefixes
- filelist_file:  file containing FASTQ file paths

Output:
- metadata.yaml: a YAML file with structured sample metadata
"""
 
import yaml
from collections import defaultdict

work_dir = "../../results/01_proc_reads"

# Input files
prefix_file = "../../data/01_proc_reads/prefixes"
filelist_file = "../../data/01_proc_reads/filelist"
output_yaml = f"{work_dir}/metadata.yaml"

# Read prefixes
with open(prefix_file) as f:
    prefixes = [x.strip() for x in f if x.strip()]

# Read filelist
with open(filelist_file) as f:
    fastq_files = [x.strip() for x in f if x.strip()]

# Dictionary: {prefix: {"R1": [...], "R2": [...]}}
samples = defaultdict(lambda: {"R1": [], "R2": []})

# Match each FASTQ to its prefix
for fq in fastq_files:
    for prefix in prefixes:
        if prefix in fq:
            if fq.endswith("_1.fq.gz") or fq.endswith("_R1.fastq.gz") or "_1.fastq.gz" in fq:
                samples[prefix]["R1"].append(fq)
            elif fq.endswith("_2.fq.gz") or fq.endswith("_R2.fastq.gz") or "_2.fastq.gz" in fq:
                samples[prefix]["R2"].append(fq)
            break  # stop checking once matched

# Build YAML structure
yaml_data = {}
samples_sorted = sorted(samples.keys(), key=lambda x: (x.split('_')[0], int(x.split('_')[2])))

for prefix in samples_sorted:
    if "FIN" in prefix:
        country = prefix.split("_")[0][1:]  # FIN from SFIN
        region = prefix.split("_")[0][0]  # S from SFIN        
    else:
        country = prefix.split("_")[0][2:]  # AU from LAAU
        region = prefix.split("_")[0][:2]  # LA from LAAU
        
    season_replicate = prefix.split("_")[1]
    season = season_replicate[0]
    replicate = season_replicate[-1]
    year = int(prefix.split("_")[2])

    # Combine R1/R2 pairs by sorting
    runs = []
    for r1 in sorted(samples[prefix]["R1"]):
        base = r1.rsplit("_1", 1)[0]  # try to find the matching _2
        r2 = [r for r in samples[prefix]["R2"] if r.startswith(base)]
        if r2:
            runs.append({"R1": r1, "R2": r2[0]})

    yaml_data.setdefault(country, {}).setdefault(region, {}).setdefault(season, {}).setdefault(year, {}).setdefault(replicate, {})
    # Add this sample
    yaml_data[country][region][season][year][replicate] = {
        "fastq": runs
    }

# Write YAML to file
with open(output_yaml, "w") as out:
    yaml.safe_dump(yaml_data, out, sort_keys=False, width=120)


