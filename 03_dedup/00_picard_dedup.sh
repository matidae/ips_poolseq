#!/usr/bin/env bash

#----------------------------------------------------------------------
# Mark duplicates and remove optical ones (flag: DT:Z:SQ) 
# using samtools (for Illumina NovaSeq X PE genomic reads)
# 
# Input: 
#   - $MAPPINGS_DIR/$prefix.sort.bam
# Output: 
#   - $DEDUP_DIR/$prefix.dedup.sort.bam
#   - $DEDUP_RESULTS/metrics/$prefix.dup_metrics.txt
#----------------------------------------------------------------------

source ./utils/paths.sh
set -euo pipefail

#Working dirs
in_dir="$MAPPINGS_DIR"
out_dir="$DEDUP_DIR"
results_dir="$DEDUP_RESULTS"
prefixes="$PREFIXES"

threads=4
jobs=10

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
        CREATE_INDEX=false

        # Remove optical duplicates (DT:Z:SQ)
        samtools view -@ "$threads" -h "$out_dir/$prefix.markdup.bam" | (grep -v "DT:Z:SQ" || true) \
        | samtools view -@ "$threads" -b -o "$out_dir/$prefix.dedup.sort.bam" -

        # Create bam index
        samtools index -@ "$threads" "$out_dir/$prefix.dedup.sort.bam"

        # Remove tmp files
        rm "$out_dir/$prefix.markdup.bam"
        rm -r "$tmpdir"
        log "done: $out_dir/$prefix.dedup.sort.bam"
}

export out_dir results_dir threads
export -f dedup_optical
export -f log
 
# Run deduplication in parallel
log "=== Deduplication start ==="
while read -r p; do echo "$in_dir/$p.sort.bam"; done < "$prefixes" \
    | parallel -j "$jobs" dedup_optical
log "=== Deduplication complete ==="