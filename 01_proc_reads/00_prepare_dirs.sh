#!/usr/bin/env bash

#----------------------------------------------------------------------
# Generates metadata files with prefixes and paths to reads
# Also creates the directory structure for downstream analysis
#
# Input: 
#   - reads in fq.gz format located in $reads_path
# Output: 
#    - prefixes: file with the prefixes of all samples 
#    - filelist: list of all fastq files 
#    - $out_dir/$prefix (directory created for each sample)
#----------------------------------------------------------------------

reads_path="../data/00_raw_reads"
out_dir="../data/01_proc_reads"

# Input
# List of files from: find "$reads_path/" -name "*.fq.gz"

# Output
out_prefixes="$out_dir/prefixes"
out_readslist="$out_dir/filelist"

# Create out directory if it doesn't exist
mkdir -p "$out_dir"

# Create a file with a list of all the reads paths
find "$reads_path/" -type f -name "*.fq.gz" | sort -V  > "$out_readslist"

prefixes_array=($(awk -F'/' '{print $(NF-1)}' "$out_readslist" | sort -u))

# Create a file with the list of all the sample prefixes 
printf "%s\n" "${prefixes_array[@]}" | sort -t_ -k1,1 -k2.1,2.1 -k3,3 -k2.2,2.2 > "$out_prefixes"

# Create a directory to process the reads of each sample
while read -r prefix; do 
   mkdir -p "$out_dir/$prefix"
done < "$out_prefixes"
