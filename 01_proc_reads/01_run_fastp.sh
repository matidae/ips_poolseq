#!/usr/bin/env bash

#----------------------------------------------------------------------
# Fastp - remove adapters, trim by quality and polyG end (no dedup)
# Run: ./01_run_fastp.sh
# In: original raw reads
# Out: qc reads, json and html reports
#----------------------------------------------------------------------

#Bulk QC processing of datasets with fastp.
while read -r prefix; do
    files=($(grep $prefix ./data/01_proc_reads/filelist))
    numfiles=$(grep -c $prefix ./data/01_proc_reads/filelist)    
    
    html=$(echo ${files[0]} | cut -f6 -d'/' | sed 's/_..fq.gz/.qc_report.html/')
    json=$(echo ${files[0]} | cut -f6 -d'/' | sed 's/_..fq.gz/.qc_report.json/')
    out1=$(echo ${files[0]} | cut -f6 -d'/' | sed 's/.fq/.qc.fq/')
    out2=$(echo ${files[1]} | cut -f6 -d'/' | sed 's/.fq/.qc.fq/')
    
    fastp -i ${files[0]} \
        -I ${files[1]} \
        -o ./data/01_proc_reads/$prefix/$out1 \
        -O ./data/01_proc_reads/$prefix/$out2 \
        -w 16 \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --length_required 50 \
        --qualified_quality_phred 20 \
        -j ./data/01_proc_reads/$prefix/$json \
        -h ./data/01_proc_reads/$prefix/$html
    
    if [ "$numfiles" -eq 4 ]; then
        html2=$(echo ${files[2]} | cut -f6 -d'/' | sed 's/_..fq.gz/.qc_report.html/')
        json2=$(echo ${files[2]} | cut -f6 -d'/' | sed 's/_..fq.gz/.qc_report.json/')
        out3=$(echo ${files[2]} | cut -f6 -d'/' | sed 's/.fq/.qc.fq/')
        out4=$(echo ${files[3]} | cut -f6 -d'/' | sed 's/.fq/.qc.fq/')
        fastp -i ${files[2]} \
            -I ${files[3]} \
            -o ./data/01_proc_reads/$prefix/$out3 \
            -O ./data/01_proc_reads/$prefix/$out4 \
            -w 16 \
            --detect_adapter_for_pe \
            --trim_poly_g \
            --length_required 50 \
            --qualified_quality_phred 20 \
            -j ./data/01_proc_reads/$prefix/$json2 \
            -h ./data/01_proc_reads/$prefix/$html2
    fi 
done < ./data/01_proc_reads/prefixes
