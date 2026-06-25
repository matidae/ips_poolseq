#!/usr/bin/env bash

#------------------------------------------------------------------------------
# Runs fastp and outputs QC reports, remove adapters, trim by quality and 
# removes polyG end.
# 
# Input:
#   - $PROC_DIR/filelist : list with FASTQ reads files for processing
#   - $PROC_DIR/prefixes : list of sample prefixes
#   - $RAW_DIR/*.fq.gz   : FASTQ files for processing
# Output:
#   - $PROC_DIR/$prefix/*.qc.fq          : processed FASTQ files
#   - $PROC_DIR/$prefix/*.qc_report.json : fastp JSON report
#   - $PROC_DIR/$prefix/*.qc_report.html : fastp HTML report
#------------------------------------------------------------------------------

source ../utils/paths.sh
set -euo pipefail

out_dir="$PROC_DIR"

threads=10 #Max fastp can use
jobs=4

export out_dir
export threads
export -f log

# Input
prefixes="$out_dir/prefixes"
filelist="$out_dir/filelist"

# Output
# fastp processed FASTQ files : $prefix.fq.qc.gz 

# Bulk QC processing of datasets with fastp.
pairs=()
while read -r prefix; do
    files=($(grep "$prefix" "$filelist"))
    n=${#files[@]}

    for ((i=0; i<n; i+=2)); do
        # Store as triplet: prefix, R1, R2
        pairs+=("$prefix|${files[i]}|${files[i+1]}")
    done
done < "$prefixes"

fastp_run(){
    IFS='|' read -r prefix R1 R2 <<< "$1"

    base1=$(basename "$R1")
    base2=$(basename "$R2")

    out1="${base1/.fq.gz/.qc.fq.gz}"
    out2="${base2/.fq.gz/.qc.fq.gz}"

    html="${base1/_1.fq.gz/.qc_report.html}"
    json="${base1/_1.fq.gz/.qc_report.json}"

    fastp \
        -i "$R1" \
        -I "$R2" \
        -o "$out_dir/$prefix/$out1" \
        -O "$out_dir/$prefix/$out2" \
        -w "$threads" \
        --detect_adapter_for_pe \
        --cut_tail \
        --cut_tail_mean_quality 20 \
        --trim_poly_g \
        --trim_poly_x \
        --length_required 50 \
        --qualified_quality_phred 20 \
        -j "$out_dir/$prefix/$json" \
        -h "$out_dir/$prefix/$html"

    log "done: $out_dir/$prefix/$json"
}

export -f fastp_run
log "=== Fastp QC start ==="
parallel -j "$jobs" fastp_run ::: "${pairs[@]}"
log "=== Fastp QC complete ==="
