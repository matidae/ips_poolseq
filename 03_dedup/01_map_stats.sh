#!/usr/bin/env bash

#----------------------------------------------------------------------
# Run multiple programs to get alignment stats
# The final ouptut will be used by 02_map_dedup_report.ipynb for an html report and plots
# Run: ./01_map_stats.sh
# Out: summary_table with stats of mapping and deduplication
#----------------------------------------------------------------------
out="../../results/03_dedup"

# Run multiple samtools alignment metrics tools
conda activate extra
for i in $(ls ../../data/03_dedup/*.sort.bam | tail -n +34); do 
    prefix=$(basename $i .dedup.sort.bam );    
    samtools flagstat $i > $out/$prefix.samtools.flagstat    
    samtools view -c -F 4 ../../data/02_mappings/$prefix.sort.bam > $out/$prefix.samtools.before_dedup    
    samtools view -q 20 -c $i > $out/$prefix.samtools.q20
    mosdepth --by 500000 -t 30  $out"/mosdepth/"$prefix.mosdepth $i
done

# Run picard CollectInsertSizeMetrics to obtain the histogram plot 
conda activate bio
find ../../data/03_dedup -name "*.sort.bam" | \
    parallel -j 5 'picard CollectInsertSizeMetrics \
        I={} \
        O=\$out/insertion_metrics/{/.}.insertion_metrics.txt \
        H=\$out/insertion_metrics/{/.}.insertion_metrics_histogram.pdf'

 ln -s $out/insertion_size_metrics/ ../../results/02_mappings/insertion_size_metrics

for i in ../../results/03_dedup/insertion_size_metrics/*.pdf;  do 
    convert -density 300 "$i" -quality 100 "${i%.pdf}.png" ; 
done


# Calculate the exact mean from mosdepth per base output
ls ../../results/03_dedup/mosdepth/*.mosdepth.per-base.bed.gz | \
    parallel -j 40 'prefix=$(basename {} .mosdepth.per-base.bed.gz); \
    depth=$(zcat {} | awk "{sum+=\$4; count++} END {printf \"%.1f\", sum/count}"); \
    echo -e "$prefix\t$depth"' > ../../results/03_dedup/mosdepth/exact_mean_coverage

# Sort the data by prefix to keep order of the rows
for i in $(cat ../../results/01_proc_reads/all_poolseq_report.tsv | grep -v Idn | cut -f1); do 
    grep -w $i ../../results/03_dedup/mosdepth/exact_mean_coverage | awk '{print $2}' ; 
done > ../../results/03_dedup/mosdepth/exact_mean_coverage.sort

# Get a list of how many reads are mapping
for i in $(cat ../../results/01_proc_reads/all_poolseq_report.tsv | grep -v Idn | cut -f1); do  
    awk '{printf "%.1f\n" , $1/1e6}'  ../../results/03_dedup/samtools_metrics/$i.samtools.before_dedup;  
done > ../../results/03_dedup/samtools_metrics/mapped_reads_before_dedup

# Get a list of how many reads remain mapped after deduplication
for i in $(cat ../../results/01_proc_reads/all_poolseq_report.tsv | grep -v Idn | cut -f1); do 
    cat ../../results/03_dedup/samtools_metrics/$i.samtools.flagstat | grep "mapped.*%" | \
    awk 'OFS="\t" {printf "%.1f\n", $1/1e6 }'; 
done > ../../results/03_dedup/samtools_metrics/mapped_reads_after_dedup

# Prepare the table with data from the previous report
awk 'OFS="\t" {print $1,$2, $4, $5, $3, $6, $15, $7}'  ../../results/01_proc_reads/all_poolseq_report.tsv | \
    grep -v Idn > ../../results/03_dedup/half_previous_table

# Write column names for summary table
echo -e "Idn\tYear\tCountry\tRegion\tTime\tRep\tRaw_reads\tQC_reads\tMapped_reads\tDedup_reads\tMean_depth" \
    > ../../results/03_dedup/summary_table

# Output a summary table with all the data needed for the html report
paste ../../results/03_dedup/half_previous_table \
    ../../results/03_dedup/samtools_metrics/mapped_reads_before_dedup \
    ../../results/03_dedup/samtools_metrics/mapped_reads_after_dedup  \
    ../../results/03_dedup/mosdepth/exact_mean_coverage.sort \
    >> ../../results/03_dedup/summary_table
