#!/usr/bin/env bash

#----------------------------------------------------------------------
# Setup directory structure and files with prefixes and paths
# Run: ./00_prepare_dirs.sh
#----------------------------------------------------------------------

##Create a file with the list of all the prefixes that will be used
for i in $(ls ./data/00_raw_reads/X204SC25022243*/01.RawData/*/*.fq.gz); do 
    prefix=$(echo $i | cut -f5 -d'/')  
    echo $prefix
done | sort | uniq > ./data/01_proc_reads/prefixes

#Create a directory for each dataset
while read -r prefix; do 
    mkdir -p ./data/01_proc_reads/$prefix
done < ./data/01_proc_reads/prefixes

#Create a list of all original files in the dataset
for i in $(ls ./data/00_raw_reads/X204SC25022243*/01.RawData/*/*.fq.gz); do 
    echo $i
done > ./data/01_proc_reads/filelist
