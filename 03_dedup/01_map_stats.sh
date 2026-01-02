#!/usr/bin/env bash

#----------------------------------------------------------------------
# Generate mapping and deduplication summary metrics for PoolSeq samples
# using  samtools, mosdepth, and Picard 
#
# Input:
#   - ../data/03_dedup/*.dedup.sort.bam
#   - ../data/02_mappings/*.sort.bam
#   - ../results/01_proc_reads/all_poolseq_report.tsv
# Output:
#   - summary_table.tsv 
#----------------------------------------------------------------------

dedup_dir="../data/03_dedup"
mappings_dir="../data/02_mappings"
proc_dir="../results/01_proc_reads"

out_dir="../results/03_dedup"
out_samtools="../results/03_dedup/samtools_metrics"
out_insert_size="../results/03_dedup/insertion_size_metrics"
out_depth_mapping="../results/03_dedup/depth_mappings" # Depth of mappings
out_depth_dedup="../results/03_dedup/depth_metrics" # Depth after deduplication
out_depth_opt="../results/03_dedup/depth_metrics_alt"   # Depth after optical deduplication only
out_temp="../results/03_dedup/temp"

jobs=10
threads=6

mkdir -p "$out_samtools" "$out_insert_size" "$out_depth_mapping" "$out_depth_dedup" "$out_depth_opt" "$out_temp"

# Run multiple samtools alignment metrics tools
eval "$(conda shell.bash hook)"
conda activate extra

for i in "$dedup_dir"/*.dedup.sort.bam ; do 
    prefix=$(basename "$i" .dedup.sort.bam );    
    samtools flagstat "$i" > "$out_samtools/$prefix.samtools.flagstat"
    samtools view -c -F 4 "$mappings_dir/$prefix.sort.bam" > "$out_samtools/$prefix.samtools.before_dedup"
done


# Calculate depth metrics with mosdepth for all mappings and deduplicated bams
ls "$dedup_dir"/*.sort.bam | parallel -j "$jobs" \
  "mosdepth --by 5000000 -t $threads --no-per-base '"$out_depth_dedup"'/{/.} {}"
    
ls "$dedup_dir"/*.sort.bam | parallel -j "$jobs" \
  "mosdepth --by 5000000 -F 0 -t $threads --no-per-base '"$out_depth_opt"'/{/.} {}"

ls "$mappings_dir"/*.sort.bam | parallel -j "$jobs" \
  "mosdepth --by 5000000 -t $threads --no-per-base '"$out_depth_mapping"'/{/.} {}"


# Run picard CollectInsertSizeMetrics to obtain the histogram plot
eval "$(conda shell.bash hook)"
conda activate bio

find "$dedup_dir" -maxdepth 1 -name "*.dedup.sort.bam" | \
    parallel -j $jobs 'picard CollectInsertSizeMetrics \
        I={} \
        O='"$out_insert_size"'/{/.}.insertion_metrics.txt \
        H='"$out_insert_size"'/{/.}.insertion_metrics_histogram.pdf'

 for i in "$out_insert_size"/*.pdf;  do
    magick convert -density 300 "$i" -quality 100 "${i%.pdf}.png" ; 
done

# List of samples
samples=$(awk 'NR>1 {print $1}' "$proc_dir/all_poolseq_report.tsv")
# Using () subshell to avoid bash code in output
# Extract the mean depth from mosdepth (first mapping, no deduplication)
(for i in $samples; do
    tail -n1 $out_depth_mapping/$i.sort.mosdepth.summary.txt | cut -f4
done > "$out_temp/exact_mean_coverage_nodedup")

# Extract the mean depth from mosdepth (all deduplication)
(for i in $samples; do
    tail -n1 $out_depth_dedup/$i.dedup.sort.mosdepth.summary.txt | cut -f4
done > "$out_temp/exact_mean_coverage_dedup")

# Extract the mean depth from mosdepth (optical only deduplication)
(for i in $samples; do
    tail -n1 $out_depth_opt/$i.dedup.sort.mosdepth.summary.txt | cut -f4
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


# Prepare the table with data from the previous report (Before, After, GC, Length)
awk 'OFS="\t" {print $1, $8, $2, $6, $5}'  "$proc_dir/all_poolseq_report.tsv" | \
    grep -v Idn > "$out_temp/half_previous_table"

# Write column names for summary table
echo -e "Idn\tRaw_reads\tQC_reads\tGC\tLength\tMapped_reads\tDedup_optical\tDedup_all\tMean_depth_dedup\
    \tMean_depth_dedup_alt\tMean_depth_nodedup" > "$out_dir/summary_table.tsv"

# Output a summary table with all the data needed for the html report
paste "$out_temp/half_previous_table" \
    "$out_temp/mapped_reads_before_dedup" \
    "$out_temp/mapped_reads_after_dedup"  \
    "$out_temp/exact_mean_coverage_dedup" \
    "$out_temp/exact_mean_coverage_dedup_alt" \
    "$out_temp/exact_mean_coverage_nodedup" \
    >> "$out_dir/summary_table.tsv"

