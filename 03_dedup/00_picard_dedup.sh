#!/usr/bin/env bash

#----------------------------------------------------------------------
# Mark duplicates and remove optical ones (flag: DT:Z:SQ) 
# using samtools (for Illumina NovaSeq X PE genomic reads)
# 
# Input: 
#       - $prefix.sort.bam
# Output: 
#       - $out_dir/$prefix.dedup.sort.bam
#       - $results_dir/metrics/$prefix.dup_metrics.txt
#----------------------------------------------------------------------
set -euo pipefail

in_dir="../data/02_mappings"
out_dir="../data/03_dedup"
results_dir="../results/03_dedup"   

threads=4
parallel=10

mkdir -p "$out_dir"
mkdir -p "$results_dir/metrics"

dedup_optical() {
    bam="$1"
    prefix=$(basename "$bam" .sort.bam )
    tmpdir="$out_dir/tmp_$prefix"
    mkdir -p "$tmpdir"

    # Mark duplicates with different flags 
    # DT:Z:LB (library duplicate) or DT:Z:SQ (sequencing duplicate)
    picard -Xmx16g MarkDuplicates \
        I="$bam" \
        O="$out_dir/$prefix.markdup.bam" \
        M="$results_dir/metrics/$prefix.dup_metrics.txt" \
        ASSUME_SORT_ORDER=coordinate \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
        VALIDATION_STRINGENCY=LENIENT \
        TAGGING_POLICY=All \
        TMP_DIR="$tmpdir" \
        CREATE_INDEX=false ;

        # Remove optical duplicates (DT:Z:SQ)
        samtools view -@ "$threads" -h "$out_dir/$prefix.markdup.bam" | grep -v "DT:Z:SQ" \
        | samtools view -@ "$threads" -b -o "$out_dir/$prefix.dedup.sort.bam" -

        # Create bam index
        samtools index -@ "$threads" "$out_dir/$prefix.dedup.sort.bam"

        # Remove tmp files
        rm "$out_dir/$prefix.markdup.bam"
        rm -r "$tmpdir"
}

export out_dir results_dir threads 
export -f dedup_optical

parallel --halt soon,fail=1 -j $parallel dedup_optical ::: "$in_dir"/*.sort.bam