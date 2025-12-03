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

parallel=10
threads=6

mkdir -p "$out_samtools" "$out_insert_size" "$out_depth"

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
    parallel --halt soon,fail=1 -j 5 'picard CollectInsertSizeMetrics \
        I={} \        
        O='"$out_insert_size"'/{/.}.insertion_metrics.txt \
        H='"$out_insert_size"'/{/.}.insertion_metrics_histogram.pdf'

 for i in "$out_insert_size"/*.pdf;  do 
    convert -density 300 "$i" -quality 100 "${i%.pdf}.png" ; 
done

# Calculate the exact mean from mosdepth per base output
find "$out_depth" -name "*.mosdepth.per-base.bed.gz" | \
    parallel -j $parallel 'prefix=$(basename {} .mosdepth.per-base.bed.gz); \
    depth=$(zcat {} | awk "{sum+=\$4; count++} END {printf \"%.1f\", sum/count}"); \
    echo -e "$prefix\t$depth"' > "$out_depth"/exact_mean_coverage

# Sort the data by prefix to keep order of the rows
for i in $(cat "$proc_dir/all_poolseq_report.tsv" | grep -v Idn | cut -f1); do 
    grep -w "$i" "$out_depth/exact_mean_coverage" | awk '{print $2}' ; 
done > "$out_depth/exact_mean_coverage.sort"

# Get a list of how many reads are mapping
for i in $(cat "$proc_dir/all_poolseq_report.tsv" | grep -v Idn | cut -f1); do  
    awk '{printf "%.1f\n" , $1/1e6}'  "$out_samtools/$i.samtools.before_dedup";
done > "$out_samtools/mapped_reads_before_dedup"

# Get a list of how many reads remain mapped after deduplication
for i in $(cat "$proc_dir/all_poolseq_report.tsv" | grep -v Idn | cut -f1); do 
    cat "$out_samtools/$i.samtools.flagstat" | grep "mapped.*%" | \
    awk 'OFS="\t" {printf "%.1f\n", $1/1e6 }'; 
done > "$out_samtools/mapped_reads_after_dedup"

# Prepare the table with data from the previous report
awk 'OFS="\t" {print $1,$2, $4, $5, $3, $6, $15, $7, $10}'  "$proc_dir/all_poolseq_report.tsv" | \
    grep -v Idn > "$out_dir/half_previous_table"

# Write column names for summary table
echo -e "Idn\tYear\tCountry\tRegion\tTime\tRep\tRaw_reads\tQC_reads\tLength\tMapped_reads\tDedup_reads\tMean_depth" \
    > "$out_dir/summary_table.tsv"

# Output a summary table with all the data needed for the html report
paste "$out_dir/half_previous_table" \
    "$out_samtools/mapped_reads_before_dedup" \
    "$out_samtools/mapped_reads_after_dedup"  \
    "$out_depth/exact_mean_coverage.sort" \
    >> "$out_dir/summary_table.tsv"

# Delete all temporal files
rm  "$out_dir/half_previous_table" \
    "$out_samtools/mapped_reads_before_dedup" \
    "$out_samtools/mapped_reads_after_dedup" \
    "$out_depth/exact_mean_coverage" \
    "$out_depth/exact_mean_coverage.sort" 

