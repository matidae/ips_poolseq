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
out_depth="../results/03_dedup/depth_metrics"
out_temp="../results/03_dedup/temp"

jobs=10
threads=6

mkdir -p "$out_samtools" "$out_insert_size" "$out_depth" "$out_temp"

# Run multiple samtools alignment metrics tools
eval "$(conda shell.bash hook)"
conda activate extra
for i in "$dedup_dir"/*.dedup.sort.bam ; do 
    prefix=$(basename "$i" .dedup.sort.bam );    
    samtools flagstat "$i" > "$out_samtools/$prefix.samtools.flagstat"
    samtools view -c -F 4 "$mappings_dir/$prefix.sort.bam" > "$out_samtools/$prefix.samtools.before_dedup"
    samtools view -q 20 -c "$i" > "$out_samtools/$prefix.samtools.q20"
    mosdepth --by 500000 -t $threads  "$out_depth/$prefix.mosdepth" "$i"
done

# Run picard CollectInsertSizeMetrics to obtain the histogram plot
eval "$(conda shell.bash hook)"
conda activate bio
find "$dedup_dir" -name "*.dedup.sort.bam" | \
    parallel -j $jobs 'picard CollectInsertSizeMetrics \
        I={} \
        O='"$out_insert_size"'/{/.}.insertion_metrics.txt \
        H='"$out_insert_size"'/{/.}.insertion_metrics_histogram.pdf'

 for i in "$out_insert_size"/*.pdf;  do
    magick convert -density 300 "$i" -quality 100 "${i%.pdf}.png" ; 
done

# Calculate the exact mean from mosdepth per base output
find "$out_depth" -name "*.mosdepth.per-base.bed.gz" | \
    parallel -j $jobs 'prefix=$(basename {} .mosdepth.per-base.bed.gz); \
    depth=$(zcat {} | awk "{sum+=\$4; count++} END {printf \"%.1f\", sum/count}"); \
    echo -e "$prefix\t$depth"' > "$out_temp"/exact_mean_coverage

#List of samples
samples=$(awk 'NR>1 {print $1}' "$proc_dir/all_poolseq_report.tsv")

# Sort the data by prefix to keep order of the rows
#Using () subshell to avoid bash code in output 
(for i in $samples; do 
    grep -w "$i" "$out_temp/exact_mean_coverage" | awk '{print $2}' ; 
done > "$out_temp/exact_mean_coverage.sort")

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


# Prepare the table with data from the previous report (Before, After, GC)
awk 'OFS="\t" {print $1, $8, $2, $6}'  "$proc_dir/all_poolseq_report.tsv" | \
    grep -v Idn > "$out_temp/half_previous_table"

# Write column names for summary table
echo -e "Idn\tRaw_reads\tQC_reads\tGC\tMapped_reads\tDedup_optical\tDedup_all\tMean_depth" \
    > "$out_dir/summary_table.tsv"

# Output a summary table with all the data needed for the html report
paste "$out_temp/half_previous_table" \
    "$out_temp/mapped_reads_before_dedup" \
    "$out_temp/mapped_reads_after_dedup"  \
    "$out_temp/exact_mean_coverage.sort" \
    >> "$out_dir/summary_table.tsv"

