#!/usr/bin/env bash

#----------------------------------------------------------------------
# Generate mapping and deduplication summary metrics for PoolSeq samples
# using samtools, mosdepth, and Picard
#
# Input:
#   - $DEDUP_DIR/$prefix.dedup.sort.bam
#   - $MAPPINGS_DIR/$prefix.sort.bam
#   - $PROC_RESULTS/all_poolseq_report.tsv
# Output:
#   - $DEDUP_RESULTS/summary_table.tsv
#----------------------------------------------------------------------

source ../utils/paths.sh
set -euo pipefail

dedup_dir="$DEDUP_DIR"
mappings_dir="$MAPPINGS_DIR"
proc_dir="$PROC_RESULTS"

out_dir="$DEDUP_RESULTS"
out_samtools="$out_dir/samtools_metrics"
out_insert_size="$out_dir/insertion_size_metrics"
out_depth_mapping="$out_dir/depth_mappings" # Depth of mappings
out_depth_dedup="$out_dir/depth_metrics" # Depth after deduplication
out_depth_opt="$out_dir/depth_metrics_alt"   # Depth after optical deduplication only
out_temp="$out_dir/temp"

prefix_file="$PREFIXES"
samples=$(cat "$prefix_file")

jobs=10
threads=6

mkdir -p "$out_samtools" "$out_insert_size" "$out_depth_mapping" "$out_depth_dedup" "$out_depth_opt" "$out_temp"

log "=== Mapping and deduplication metrics start ==="

# Run multiple samtools alignment metrics tools

for prefix in $samples; do
    i="$dedup_dir/$prefix.dedup.sort.bam"
    samtools flagstat "$i" > "$out_samtools/$prefix.samtools.flagstat"
    samtools view -c -F 4 "$mappings_dir/$prefix.sort.bam" > "$out_samtools/$prefix.samtools.before_dedup"
done
log "done: samtools flagstat and read counts"

# Calculate depth metrics with mosdepth for all mappings and deduplicated bams
for prefix in $samples; do echo "$prefix"; done | \
    parallel -j "$jobs" "mosdepth --by 500000 -t $threads --no-per-base $out_depth_dedup/{} $dedup_dir/{}.dedup.sort.bam"
log "done: mosdepth depth after deduplication"

for prefix in $samples; do echo "$prefix"; done | \
    parallel -j "$jobs" "mosdepth --by 500000 -F 772 -t $threads --no-per-base $out_depth_opt/{} $dedup_dir/{}.dedup.sort.bam"
log "done: mosdepth depth after optical deduplication"

for prefix in $samples; do echo "$prefix"; done | \
    parallel -j "$jobs" "mosdepth --by 500000 -t $threads --no-per-base $out_depth_mapping/{} $mappings_dir/{}.sort.bam"
log "done: mosdepth depth before deduplication"

# Run picard CollectInsertSizeMetrics to obtain the histogram plot
for prefix in $samples; do echo "$prefix"; done | \
    parallel -j $jobs "picard CollectInsertSizeMetrics \
        I=$dedup_dir/{}.dedup.sort.bam \
        O=$out_insert_size/{}.dedup.insertion_metrics.txt \
        H=$out_insert_size/{}.dedup.insertion_metrics_histogram.pdf"

log "done: picard insert size metrics"

 for i in "$out_insert_size"/*.pdf;  do
    magick convert -density 300 "$i" -quality 100 "${i%.pdf}.png" ; 
done
log "done: insert size histogram plots"


# Using () subshell to avoid bash code in output
# Extract the mean depth from mosdepth (first mapping, no deduplication)
(for i in $samples; do
    tail -n1 $out_depth_mapping/$i.mosdepth.summary.txt | cut -f4
done > "$out_temp/exact_mean_coverage_nodedup")

# Extract the mean depth from mosdepth (all deduplication)
(for i in $samples; do
    tail -n1 $out_depth_dedup/$i.mosdepth.summary.txt | cut -f4
done > "$out_temp/exact_mean_coverage_dedup")

# Extract the mean depth from mosdepth (optical only deduplication)
(for i in $samples; do
    tail -n1 $out_depth_opt/$i.mosdepth.summary.txt | cut -f4
done > "$out_temp/exact_mean_coverage_dedup_alt")

# Get a list of how many reads are mapping
(for i in $samples; do 
    awk '{printf "%.1f\n" , $1/1e6}'  "$out_samtools/$i.samtools.before_dedup";
done > "$out_temp/mapped_reads_before_dedup")


# Get a list of how many reads remain mapped after deduplication
(for i in $samples; do
    # extract primary mapped reads
    mapped=$(grep "primary mapped" "$out_samtools/$i.samtools.flagstat" | awk '{print $1}')
    # extract primary duplicates
    dup=$(grep "primary duplicates" "$out_samtools/$i.samtools.flagstat" | awk '{print $1}')
    # calculate usable reads
    usable=$((mapped - dup))
    # print both in millions
    awk -v m="$mapped" -v u="$usable" 'BEGIN{printf "%.1f\t%.1f\n", m/1e6, u/1e6}'
done > "$out_temp/mapped_reads_after_dedup")
log "done: depth and read count summaries"

# Prepare the table with data from the previous report (Before, After, GC, Length)
awk 'OFS="\t" {print $1, $8, $2, $6, $5}'  "$proc_dir/all_poolseq_report.tsv" | \
    grep -Ff "$prefix_file" > "$out_temp/half_previous_table"

# Write column names for summary table
out_table="$out_dir/summary_table.tsv"
[ -f "$out_table" ] && out_table="$out_dir/summary_table_new.tsv"

echo -e "Idn\tRaw_reads\tQC_reads\tGC\tLength\tMapped_reads\tDedup_optical\tDedup_all\tDepth_dedup\tDepth_dedup_alt\tDepth_nodedup" > "$out_table"
# Output a summary table with all the data needed for the html report
paste "$out_temp/half_previous_table" \
    "$out_temp/mapped_reads_before_dedup" \
    "$out_temp/mapped_reads_after_dedup"  \
    "$out_temp/exact_mean_coverage_dedup" \
    "$out_temp/exact_mean_coverage_dedup_alt" \
    "$out_temp/exact_mean_coverage_nodedup" \
    >> "$out_table"
log "done: $out_table"

log "=== Mapping and deduplication metrics complete ==="