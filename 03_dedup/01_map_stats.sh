#!/usr/bin/env bash

#----------------------------------------------------------------------
# Generate mapping and deduplication summary metrics for PoolSeq samples
# using samtools, mosdepth, and Picard
#
# Input:
#   - $DEDUP_DIR/$prefix.dedup.sort.bam
#   - $MAPPINGS_DIR/$prefix.sort.bam
#   - $PROC_RESULTS/all_poolseq_report.tsv
#   - $PREFIXES : sample prefix list
# Output:
#   - $DEDUP_RESULTS/summary_table.tsv
#   - $DEDUP_RESULTS/samtools_metrics/$prefix.samtools.flagstat
#   - $DEDUP_RESULTS/samtools_metrics/$prefix.samtools.before_dedup
#   - $DEDUP_RESULTS/depth_metrics/$prefix.mosdepth.*
#   - $DEDUP_RESULTS/depth_metrics_alt/$prefix.mosdepth.*
#   - $DEDUP_RESULTS/depth_mappings/$prefix.mosdepth.*
#   - $DEDUP_RESULTS/insertion_size_metrics/$prefix.dedup.insertion_metrics.txt
#   - $DEDUP_RESULTS/insertion_size_metrics/$prefix.dedup.insertion_metrics_histogram.pdf
#   - $DEDUP_RESULTS/insertion_size_metrics/$prefix.dedup.insertion_metrics_histogram.png
#----------------------------------------------------------------------

source ./utils/paths.sh
set -euo pipefail

dedup_dir="$DEDUP_DIR"
mappings_dir="$MAPPINGS_DIR"
proc_dir="$PROC_RESULTS"

out_dir="$DEDUP_RESULTS"
out_samtools="$out_dir/samtools_metrics"
out_insert_size="$out_dir/insertion_size_metrics"
out_depth_mapping="$out_dir/depth_mappings"
out_depth_dedup="$out_dir/depth_metrics"
out_depth_opt="$out_dir/depth_metrics_alt"
out_temp="$out_dir/temp"

prefix_file="$PREFIXES"
samples=$(cat "$prefix_file")

jobs=10
threads=6

mkdir -p "$out_samtools" "$out_insert_size" "$out_depth_mapping" "$out_depth_dedup" "$out_depth_opt" "$out_temp"

log "=== Mapping and deduplication metrics start ==="

# Run samtools alignment metrics
for prefix in $samples; do
    i="$dedup_dir/$prefix.dedup.sort.bam"
    samtools flagstat "$i" > "$out_samtools/$prefix.samtools.flagstat"
    samtools view -c -F 4 "$mappings_dir/$prefix.sort.bam" > "$out_samtools/$prefix.samtools.before_dedup"
done
log "done: samtools flagstat and read counts"

# mosdepth: skip samples whose output already exists under any naming convention
for prefix in $samples; do echo "$prefix"; done | \
    parallel -j "$jobs" \
    "ls $out_depth_dedup/{}.mosdepth.summary.txt 2>/dev/null | grep -q . || \
     mosdepth --by 500000 -t $threads $out_depth_dedup/{} $dedup_dir/{}.dedup.sort.bam"
log "done: mosdepth depth after deduplication"

for prefix in $samples; do echo "$prefix"; done | \
    parallel -j "$jobs" \
    "ls $out_depth_opt/{}.mosdepth.summary.txt 2>/dev/null | grep -q . || \
     mosdepth --by 500000 -F 772 -t $threads $out_depth_opt/{} $dedup_dir/{}.dedup.sort.bam"
log "done: mosdepth depth after optical deduplication"

#for prefix in $samples; do echo "$prefix"; done | \
    parallel -j "$jobs" \
    "ls $out_depth_mapping/{}.mosdepth.summary.txt 2>/dev/null | grep -q . || \
     mosdepth --by 500000 -t $threads $out_depth_mapping/{} $mappings_dir/{}.sort.bam"
log "done: mosdepth depth before deduplication"

# Calculate percentage of reads mapped
(for i in $samples; do
    raw=$(cat "$out_samtools/$i.samtools.before_dedup")
    mapped=$(grep "primary mapped" "$out_samtools/$i.samtools.flagstat" | awk '{print $1}')
    awk -v r="$raw" -v m="$mapped" 'BEGIN{printf "%.1f\n", (m/r)*100}'
done > "$out_temp/pct_mapped")
log "done: percentage of reads mapped"

# Picard: skip samples whose output already exists under any naming convention
for prefix in $samples; do echo "$prefix"; done | \
    parallel -j $jobs \
    "ls $out_insert_size/{}*.insertion_metrics.txt 2>/dev/null | grep -q . || \
     picard CollectInsertSizeMetrics \
        I=$dedup_dir/{}.dedup.sort.bam \
        O=$out_insert_size/{}.dedup.sort.insertion_metrics.txt \
        H=$out_insert_size/{}.dedup.sort.insertion_metrics_histogram.pdf"
log "done: picard insert size metrics"

for i in "$out_insert_size"/*.pdf; do
    [ -e "$i" ] || continue
    magick -density 300 "$i" -quality 100 "${i%.pdf}.png"
done
log "done: insert size histogram plots"

# Extract mean depth - find handles both old and new naming conventions
(for i in $samples; do
    f=$(find "$out_depth_mapping" -maxdepth 1 -name "$i*.mosdepth.summary.txt" | head -1)
    tail -n1 "$f" | cut -f4
done > "$out_temp/exact_mean_coverage_nodedup")

(for i in $samples; do
    f=$(find "$out_depth_dedup" -maxdepth 1 -name "$i*.mosdepth.summary.txt" | head -1)
    tail -n1 "$f" | cut -f4
done > "$out_temp/exact_mean_coverage_dedup")

(for i in $samples; do
    f=$(find "$out_depth_opt" -maxdepth 1 -name "$i*.mosdepth.summary.txt" | head -1)
    tail -n1 "$f" | cut -f4
done > "$out_temp/exact_mean_coverage_dedup_alt")

# Extract mean insert size from picard - find handles both naming conventions
# MEAN_INSERT_SIZE is field 6 in picard InsertSizeMetrics output
(for i in $samples; do
    f=$(find "$out_insert_size" -maxdepth 1 -name "$i*.insertion_metrics.txt" | head -1)
    grep -A1 "^MEDIAN_INSERT" "$f" | tail -1 | cut -f6 | awk '{printf "%.0f\n", $1}'
done > "$out_temp/mean_insert_size")

# Get a list of how many reads are mapping
(for i in $samples; do
    awk '{printf "%.1f\n", $1/1e6}' "$out_samtools/$i.samtools.before_dedup"
done > "$out_temp/mapped_reads_before_dedup")

# Get a list of how many reads remain mapped after deduplication
(for i in $samples; do
    mapped=$(grep "primary mapped" "$out_samtools/$i.samtools.flagstat" | awk '{print $1}')
    dup=$(grep "primary duplicates" "$out_samtools/$i.samtools.flagstat" | awk '{print $1}')
    usable=$((mapped - dup))
    awk -v m="$mapped" -v u="$usable" 'BEGIN{printf "%.1f\t%.1f\n", m/1e6, u/1e6}'
done > "$out_temp/mapped_reads_after_dedup")
log "done: depth and read count summaries"

# Prepare the table with data from the previous report (Before, After, GC, Length)
awk 'OFS="\t" {print $1, $8, $2, $6, $5}' "$proc_dir/all_poolseq_report.tsv" | \
    grep -Ff "$prefix_file" > "$out_temp/half_previous_table"

# Write column names for summary table
out_table="$out_dir/summary_table.tsv"
[ -f "$out_table" ] && out_table="$out_dir/summary_table.tsv"

echo -e "Idn\tRaw_reads\tQC_reads\tGC\tLength\tMapped_reads\tPct_mapped\tDedup_optical\tDedup_all\tDepth_dedup\tDepth_dedup_alt\tDepth_nodedup\tMean_insert_size" > "$out_table"

# Output a summary table with all the data needed for the html report
paste "$out_temp/half_previous_table" \
    "$out_temp/mapped_reads_before_dedup" \
    "$out_temp/pct_mapped" \
    "$out_temp/mapped_reads_after_dedup" \
    "$out_temp/exact_mean_coverage_dedup" \
    "$out_temp/exact_mean_coverage_dedup_alt" \
    "$out_temp/exact_mean_coverage_nodedup" \
    "$out_temp/mean_insert_size" \
    >> "$out_table"
log "done: $out_table"

log "=== Mapping and deduplication metrics complete ==="