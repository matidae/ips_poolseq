#!/usr/bin/env bash

#----------------------------------------------------------------------
# Basic Picard MarkDuplicates 
# Run: ./00_picard_dedup.sh
# In: ./ips_project/data/02_mappings/$prefix.sort.bam
# Out: ./ips_project/data/03_dedup/$prefix*.dup_metrics.txt
#----------------------------------------------------------------------


for i in $(ls ../../data/02_mappings/*.sort.bam); do 
    prefix=$(basename $i .sort.bam );
    picard -Xmx16g MarkDuplicates \
        I=$i \
        O=../../data/03_dedup/$prefix.dedup.sort.bam \
        ASSUME_SORT_ORDER=coordinate \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500\
	VALIDATION_STRINGENCY=LENIENT\
	REMOVE_DUPLICATES=true \
	CREATE_INDEX=true \
        M=../../results/03_dedup/markDup_metrics/$prefix.dedup_metrics.txt;   
done
