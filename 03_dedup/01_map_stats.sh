#!/usr/bin/env bash

#----------------------------------------------------------------------
# Run multiple program to get alignment stats
# Run: ./01_map_stats.sh
#----------------------------------------------------------------------
out="../../results/03_dedup"

for i in $(ls ../../data/03_dedup/*.sort.bam | tail -n +34); do 
    prefix=$(basename $i .dedup.sort.bam );    
    samtools flagstat $i > $out/$prefix.samtools.flagstat    
    samtools view -c -F 4 ../../data/02_mappings/$prefix.sort.bam > $out/$prefix.samtools.before_dedup    
    samtools view -q 20 -c $i > $out/$prefix.samtools.q20
    mosdepth --by 500000 -t 30  $out/$prefix.mosdepth $i"/mosdepth"
     
    
    
done



#picard CollectWgsMetrics I=$i O=$out/$prefix.CollectWgsMetrics \
#    R=../../data/reference/Ips_typograpgus_LG16corrected.final.fasta
    
#picard CollectInsertSizeMetrics I=$i O=$out/$prefix.CollectInsertSizeMetrics \
#    H=$out/$prefix.CollectInsertSizeMetrics.histogram.pdf

#find ../../data/03_dedup -name "*.sort.bam" | parallel -j 5 'picard CollectWgsMetrics I={} O=../results/wgs_metrics/{/.}.wgs_metrics.txt R=../../data/reference/Ips_typograpgus_LG16corrected.final.fasta'
 