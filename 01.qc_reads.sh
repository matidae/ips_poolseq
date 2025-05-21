#!/usr/bin/env bash

#--------------
# Setup directory structure and files with prefixes and paths
#---------------

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

#--------------
# Fastp - trim and discard (no dedup)
#--------------

#Bulk QC processing of datasets with fastp.
while read -r prefix; do
    files=($(grep $prefix ./ips_reads/filelist))
    numfiles=$(grep -c $prefix ./ips_reads/filelist)    
    
    html=$(echo ${files[0]} | cut -f5 -d'/' | sed 's/_..fq.gz/.qc_report.html/')
    json=$(echo ${files[0]} | cut -f5 -d'/' | sed 's/_..fq.gz/.qc_report.json/')
    out1=$(echo ${files[0]} | cut -f5 -d'/' | sed 's/.fq/.qc.fq/')
    out2=$(echo ${files[1]} | cut -f5 -d'/' | sed 's/.fq/.qc.fq/')
    
    fastp -i ${files[0]} \
        -I ${files[1]} \
        -o ./ips_reads/$prefix/$out1 \
        -O ./ips_reads/$prefix/$out2 \
        -w 16 \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --length_required 50 \
        --qualified_quality_phred 20 \
        -j ./ips_reads/$prefix/$json \
        -h ./ips_reads/$prefix/$html
    
    if [ "$numfiles" -eq 4 ]; then
        html2=$(echo ${files[2]} | cut -f5 -d'/' | sed 's/_..fq.gz/.qc_report.html/')
        json2=$(echo ${files[0]} | cut -f5 -d'/' | sed 's/_..fq.gz/.qc_report.json/')
        out3=$(echo ${files[2]} | cut -f5 -d'/' | sed 's/.fq/.qc.fq/')
        out4=$(echo ${files[3]} | cut -f5 -d'/' | sed 's/.fq/.qc.fq/')
        fastp -i ${files[2]} \
            -I ${files[3]} \
            -o ./ips_reads/$prefix/$out3 \
            -O ./ips_reads/$prefix/$out4 \
            -w 16 \
            --detect_adapter_for_pe \
            --trim_poly_g \
            --length_required 50 \
            --qualified_quality_phred 20 \
            -j ./ips_reads/$prefix/$json2 \
            -h ./ips_reads/$prefix/$html2
    fi 
done < ips_reads/prefixes

