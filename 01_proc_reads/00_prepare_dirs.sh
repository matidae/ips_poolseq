#!/usr/bin/env bash

#----------------------------------------------------------------------
# Setup directory structure and files with prefixes and paths
# Run: ./00_prepare_dirs.sh
#----------------------------------------------------------------------

##Create a file with the list of all the prefixes that will be used
for i in $(ls ips_poolseq/X204SC25022243*/01.RawData/*/*.fq.gz); do 
    prefix=$(echo $i | cut -f4 -d'/')  
    echo $prefix
done | sort | uniq > ips_reads/prefixes

#Create a directory for each dataset
while read -r prefix; do 
    mkdir -p ./ips_reads/$prefix
done < ips_reads/prefixes

#Create a list of all original files in the dataset
for i in $(ls ips_poolseq/X204SC25022243*/01.RawData/*/*.fq.gz); do 
    echo $i
done > ips_reads/filelist