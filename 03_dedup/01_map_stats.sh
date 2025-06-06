#!/usr/bin/env bash

#----------------------------------------------------------------------
# Run multiple program to get alignment stats
# Run: ./01_map_stats.sh
#----------------------------------------------------------------------
out="../../results/03_dedup"

#Run multiple samtools alignment metrics tools
conda activate extra
for i in $(ls ../../data/03_dedup/*.sort.bam | tail -n +34); do 
    prefix=$(basename $i .dedup.sort.bam );    
    samtools flagstat $i > $out/$prefix.samtools.flagstat    
    samtools view -c -F 4 ../../data/02_mappings/$prefix.sort.bam > $out/$prefix.samtools.before_dedup    
    samtools view -q 20 -c $i > $out/$prefix.samtools.q20
    mosdepth --by 500000 -t 30  $out"/mosdepth/"$prefix.mosdepth $i
done

#Run picard CollectInsertSizeMetrics to obtain the histogram plot 
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


#Calculate the exact mean from mosdepth per base output
ls ../../results/03_dedup/mosdepth/*.mosdepth.per-base.bed.gz | \
    parallel -j 40 'prefix=$(basename {} .mosdepth.per-base.bed.gz); \
    depth=$(zcat {} | awk "{sum+=\$4; count++} END {printf \"%.1f\", sum/count}"); \
    echo -e "$prefix\t$depth"' > ../../results/03_dedup/mosdepth/exact_mean_coverage

#Prepare the tables
awk 'OFS="\t" {print $1,$2, $4, $5, $3, $6, $15, $7}'  ../../results/01_proc_reads/all_poolseq_report.tsv | \
    grep -v Idn > ../../results/03_dedup/half_previous_table

for i in $(cat ../../results/01_proc_reads/all_poolseq_report.tsv | grep -v Idn | cut -f1); do 
    grep -w $i ../../results/03_dedup/mosdepth/exact_mean_coverage | awk '{print $2}' ; 
done > ../../results/03_dedup/mosdepth/exact_mean_coverage.sort

for i in $(cat ../../results/01_proc_reads/all_poolseq_report.tsv | grep -v Idn | cut -f1); do  
    awk '{printf "%.1f\n" , $1/1e6}'  ../../results/03_dedup/samtools_metrics/$i.samtools.before_dedup;  
done > ../../results/03_dedup/samtools_metrics/nreads_before_dedup_all

for i in $(cat ../../results/01_proc_reads/all_poolseq_report.tsv | grep -v Idn | cut -f1); do 
    cat ../../results/03_dedup/samtools_metrics/$i.samtools.flagstat | grep "mapped.*%" | \
    sed 's/% : N\/A)//; s/ + 0 mapped (/ /' | awk 'OFS="\t" {printf "%.1f\t%.1f\n", $1/1e6, $2}'; 
    done > ../../results/03_dedup/samtools_metrics/reads_mapped_all