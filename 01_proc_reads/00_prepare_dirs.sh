#!/usr/bin/env bash

#----------------------------------------------------------------------
# Generates files  with prefixes and paths to samples reads and creates
# the directory structure for downstream analysis
#
# Output: 
#    - prefixes: file with the prefixes of all samples 
#    - filelist: list of all fastq files 
#----------------------------------------------------------------------

reads_path="../../data/00_raw_reads/X204SC25022243*/01.RawData/*/*.fq.gz"
out_path="./data/01_proc_reads"

#Create a file with the list of all the prefixes that will be used
for i in $reads_path; do 
    prefix=$(basename "$(dirname "$i")")
    printf "%s\n" "$prefix"    
done | sort | uniq > "$out_path/prefixes"

#Create a directory for each dataset
while read -r prefix; do 
    mkdir -p "$out_path"/$prefix
done < "$out_path/prefixes"

#Create a list of all original files in the dataset
for i in $reads_path; do 
    printf "%s\n" "$i"
done > "$out_path/filelist"
