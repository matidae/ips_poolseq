#!/usr/bin/env bash

#--------------
# Mapping reads with BWA, sorting, merging and indexing
#--------------

while read -r prefix; do
    files=($(grep "$prefix" ./ips_reads/filelist))
    numfiles=$(grep -c "$prefix" ./ips_reads/filelist)    
    
    read1=$(echo ${files[0]} | cut -f5 -d'/' | sed 's/.fq/.qc.fq/')
    read2=$(echo ${files[1]} | cut -f5 -d'/' | sed 's/.fq/.qc.fq/')

    if [ "$numfiles" -eq 2 ]; then
       bwa-mem2 mem \
        -t 40 \
        -R "@RG\tID:$prefix\tSM:$prefix" \
        -o ./ips_mapping/"$prefix".sam \
        ./ips_mapping/index/Ips_typograpgus_LG16corrected.final.bwa_index \
        ./ips_reads/"$prefix"/"$read1" ./ips_reads/"$prefix"/"$read2" \
        2>./ips_mapping/"$prefix".bwa.log
        sambamba view -q -t 40 -S -f bam -o ./ips_mapping/"$prefix".bam ./ips_mapping/"$prefix".sam
        sambamba sort -q -t 40 -o ./ips_mapping/"$prefix".sort.bam ./ips_mapping/"$prefix".bam
        sambamba index -q -t 40 ./ips_mapping/"$prefix".sort.bam
        rm ./ips_mapping/"$prefix".sam
        rm ./ips_mapping/"$prefix".bam
    
    elif [ "$numfiles" -eq 4 ]; then
        read3=$(echo ${files[2]} | cut -f5 -d'/' | sed 's/.fq/.qc.fq/')
        read4=$(echo ${files[3]} | cut -f5 -d'/' | sed 's/.fq/.qc.fq/')
        
        bwa-mem2 mem \
        -t 40 \
        -R "@RG\tID:$prefix\tSM:$prefix" \
        -o ./ips_mapping/"$prefix".1.sam \
        ./ips_mapping/index/Ips_typograpgus_LG16corrected.final.bwa_index \
        ./ips_reads/"$prefix"/"$read1" ./ips_reads/"$prefix"/"$read2" \
        2>./ips_mapping/"$prefix".bwa.log

        bwa-mem2 mem \
        -t 40 \
        -R "@RG\tID:$prefix\tSM:$prefix" \
        -o ./ips_mapping/"$prefix".2.sam \
        ./ips_mapping/index/Ips_typograpgus_LG16corrected.final.bwa_index \
        ./ips_reads/"$prefix"/"$read3" ./ips_reads/"$prefix"/"$read4" \
        2>>./ips_mapping/"$prefix".bwa.log

        sambamba view -q -t 40 -S -f bam -o ./ips_mapping/"$prefix".1.bam ./ips_mapping/"$prefix".1.sam
        sambamba view -q -t 40 -S -f bam -o ./ips_mapping/"$prefix".2.bam ./ips_mapping/"$prefix".2.sam
        sambamba sort -q -t 40 -o ./ips_mapping/"$prefix".sort.1.bam ./ips_mapping/"$prefix".1.bam
        sambamba sort -q -t 40 -o ./ips_mapping/"$prefix".sort.2.bam ./ips_mapping/"$prefix".2.bam
        sambamba merge -q -t 40 ./ips_mapping/"$prefix".sort.bam ./ips_mapping/"$prefix".sort.1.bam ./ips_mapping/"$prefix".sort.2.bam
        rm ./ips_mapping/"$prefix".1.sam ./ips_mapping/"$prefix".2.sam
        rm ./ips_mapping/"$prefix".1.bam ./ips_mapping/"$prefix".2.bam 
        rm ./ips_mapping/"$prefix".sort.1.bam ./ips_mapping/"$prefix".sort.2.bam 
        sambamba index -q -t 40 ./ips_mapping/"$prefix".sort.bam
        rm ./ips_mapping/"$prefix".sort.1.bam ./ips_mapping/"$prefix".sort.1.bam.bai 
        rm ./ips_mapping/"$prefix".sort.1.bam ./ips_mapping/"$prefix".sort.2.bam.bai 

    fi 

done < ips_reads/prefixes