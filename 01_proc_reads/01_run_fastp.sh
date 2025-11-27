#!/usr/bin/env bash

#------------------------------------------------------------------------------
# Runs fastp and outputs QC reports, remove adapters, trim by quality and 
# removes polyG end (no dedup).
# 
# Input:
#   - filelist: list with FASTQ reads files for processing
#   - prefixes: list of sample prefixes
#   - FASTQ files for processing
# Output:
#   - *.qc.fq          : processed FASTQ files
#   - *.qc_report.json : fastp JSON report
#   - *.qc_report.html : fastp HTML report
#------------------------------------------------------------------------------

out_dir="../data/01_proc_reads"

# Input
prefixes="$out_dir/prefixes"
filelist="$out_dir/filelist"
# FASTQ files for processing : $prefix.fq.gz 

# Output
# fastp processed FASTQ files :$prefix.fq.qc.gz 

# Bulk QC processing of datasets with fastp.
while read -r prefix; do
    files=($(grep "$prefix" "$filelist"))
    n=${#files[@]}

    for ((i=0; i<n; i+=2)); do
        R1="${files[$i]}"
        R2="${files[$i+1]}"

        base1=$(basename "$R1")
        base2=$(basename "$R2")

        out1="${base1/.fq.gz/.qc.fq}"
        out2="${base2/.fq.gz/.qc.fq}"

        html="${base1/_..fq.gz/.qc_report.html}"
        json="${base1/_..fq.gz/.qc_report.json}"

        fastp \
            -i "../$R1" \
            -I "../$R2" \
            -o "$out_dir/$prefix/$out1" \
            -O "$out_dir/$prefix/$out2" \
            -w 16 \
            --detect_adapter_for_pe \
            --cut_tail \
            --cut_tail_mean_quality 20 \
            --trim_poly_g \
            --trim_poly_x \
            --length_required 50 \
            --qualified_quality_phred 20 \
            -j "$out_dir/$prefix/$json" \
            -h "$out_dir/$prefix/$html"

    done
done < "$prefixes"
