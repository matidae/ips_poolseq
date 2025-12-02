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


in_dir="../data/02_mappings"
out_dir="../data/03_dedup"
results_dir="../results/03_dedup"   

threads=4
parallel=10

mkdir -p "$out_dir"
mkdir -p "$results_dir/metrics"

export out_dir results_dir threads

parallel -j $parallel '
    bam="{}"
    prefix=$(basename "$bam" .sort.bam );
    # Mark duplicates with different flags (do not remove)
    picard -Xmx16g MarkDuplicates \
        I={} \
        O="$out_dir/$prefix.markdup.bam" \
        M="$results_dir/metrics/$prefix.dedup_metrics.txt" \
        ASSUME_SORT_ORDER=coordinate \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
        VALIDATION_STRINGENCY=LENIENT \
        TAGGING_POLICY=All \
        CREATE_INDEX=false ;        

        # Remove optical duplicates (DT:Z:SQ)
        samtools view -h "$out_dir/$prefix.markdup.bam" | grep -v "DT:Z:SQ" \
        | samtools view -b -o "$out_dir/$prefix.dedup.sort.bam" -

        # Create bam index
        samtools index -@ $threads "$out_dir/$prefix.dedup.sort.bam"
' ::: "$in_dir"/*.sort.bam

