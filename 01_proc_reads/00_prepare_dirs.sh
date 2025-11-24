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
#    - $out_path/$prefix (directory created for each sample)
#----------------------------------------------------------------------

reads_path="../data/00_raw_reads"
out_path="../data/01_proc_reads"

#Create out directory if it doesn' t exist
mkdir -p "$out_path"

#Create a file with a list of all the reads paths
find "$reads_path/" -type f -name "*.fq.gz" > "$out_path/filelist"

prefixes=($(awk -F'/' '{print $(NF-1)}' "$out_path/filelist" | sort -Vu))

#Create a file with the list of all the sample prefixes 
printf "%s\n" "${prefixes[@]}" > "$out_path/prefixes"

#Create a directory to process the reads of each sample
while read -r prefix; do 
   mkdir -p "$out_path/$prefix"
done < "$out_path/prefixes"
