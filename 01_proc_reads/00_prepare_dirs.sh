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
prefixes="$out_dir/prefixes"
filelist="$out_dir/filelist"

# Create out directory if it doesn' t exist
mkdir -p "$out_dir"

# Create a file with a list of all the reads paths
find "$reads_path/" -type f -name "*.fq.gz" > "$filelist"

prefixes=($(awk -F'/' '{print $(NF-1)}' "$filelist" | sort -Vu))

# Create a file with the list of all the sample prefixes 
printf "%s\n" "${prefixes[@]}" > "$prefixes"

# Create a directory to process the reads of each sample
while read -r prefix; do 
   mkdir -p "$out_dir/$prefix"
done < "$prefixes"
