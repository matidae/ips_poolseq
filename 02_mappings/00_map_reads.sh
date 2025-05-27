#!/usr/bin/env bash

#----------------------------------------------------------------------
# Mapping reads with BWA, sorting, merging and indexing with sambamba
# Run: ./00_map_reads.sh
# In: qc reads and genome
# Out: alignment sort.bam files
#----------------------------------------------------------------------


while read -r prefix; do
    files=($(grep "$prefix" ./data/01_proc_reads/filelist))
    numfiles=$(grep -c "$prefix" ./data/01_proc_reads/filelist)    
    
    read1=$(echo ${files[0]} | cut -f6 -d'/' | sed 's/.fq/.qc.fq/')
    read2=$(echo ${files[1]} | cut -f6 -d'/' | sed 's/.fq/.qc.fq/')

    if [ "$numfiles" -eq 2 ]; then
        #Run bwa for the first set of PE reads
        bwa-mem2 mem \
        -t 40 \
        -R "@RG\tID:$prefix\tSM:$prefix" \
        -o ./data/02_mappings/"$prefix".sam \
        ./data/02_mappings/index/Ips_typograpgus_LG16corrected.final.bwa_index \
        ./data/01_proc_reads/"$prefix"/"$read1" ./data/01_proc_reads/"$prefix"/"$read2" \
        2>./data/02_mappings/"$prefix".bwa.log
        
        # Sam files to bam format
        sambamba view -q -t 40 -S -f bam -o ./data/02_mappings/"$prefix".bam ./data/02_mappings/"$prefix".sam
        
        # Sort bam
        sambamba sort -q -t 40 -o ./data/02_mappings/"$prefix".sort.bam ./data/02_mappings/"$prefix".bam
        
        # Index bam file
        sambamba index -q -t 40 ./data/02_mappings/"$prefix".sort.bam

        # Delete sam and bam files not sorted
        rm ./data/02_mappings/"$prefix".sam
        rm ./data/02_mappings/"$prefix".bam
    
    elif [ "$numfiles" -eq 4 ]; then
        read3=$(echo ${files[2]} | cut -f6 -d'/' | sed 's/.fq/.qc.fq/')
        read4=$(echo ${files[3]} | cut -f6 -d'/' | sed 's/.fq/.qc.fq/')
        
        #Run bwa for the first set of PE reads
        bwa-mem2 mem \
        -t 40 \
        -R "@RG\tID:$prefix\tSM:$prefix" \
        -o ./data/02_mappings/"$prefix".1.sam \
        ./data/02_mappings/index/Ips_typograpgus_LG16corrected.final.bwa_index \
        ./data/01_proc_reads/"$prefix"/"$read1" ./data/01_proc_reads/"$prefix"/"$read2" \
        2>./data/02_mappings/"$prefix".bwa.log

        #Run bwa for the second set of PE reads
        bwa-mem2 mem \
        -t 40 \
        -R "@RG\tID:$prefix\tSM:$prefix" \
        -o ./data/02_mappings/"$prefix".2.sam \
        ./data/02_mappings/index/Ips_typograpgus_LG16corrected.final.bwa_index \
        ./data/01_proc_reads/"$prefix"/"$read3" ./data/01_proc_reads/"$prefix"/"$read4" \
        2>>./data/02_mappings/"$prefix".bwa.log

        # Sam files to bam format
        sambamba view -q -t 40 -S -f bam -o ./data/02_mappings/"$prefix".1.bam ./data/02_mappings/"$prefix".1.sam
        sambamba view -q -t 40 -S -f bam -o ./data/02_mappings/"$prefix".2.bam ./data/02_mappings/"$prefix".2.sam
        
        # Sort bam
        sambamba sort -q -t 40 -o ./data/02_mappings/"$prefix".sort.1.bam ./data/02_mappings/"$prefix".1.bam
        sambamba sort -q -t 40 -o ./data/02_mappings/"$prefix".sort.2.bam ./data/02_mappings/"$prefix".2.bam
        
        # Merge both sorted bam files
        sambamba merge -q -t 40 ./data/02_mappings/"$prefix".sort.bam ./data/02_mappings/"$prefix".sort.1.bam ./data/02_mappings/"$prefix".sort.2.bam
        
         # Index bam
        sambamba index -q -t 40 ./data/02_mappings/"$prefix".sort.bam

        # Delete sam, bam, sort and index files before merging
        rm ./data/02_mappings/"$prefix".1.sam ./data/02_mappings/"$prefix".2.sam
        rm ./data/02_mappings/"$prefix".1.bam ./data/02_mappings/"$prefix".2.bam         
        rm ./data/02_mappings/"$prefix".sort.1.bam ./data/02_mappings/"$prefix".sort.1.bam.bai 
        rm ./data/02_mappings/"$prefix".sort.2.bam ./data/02_mappings/"$prefix".sort.2.bam.bai 
        

    fi 

done < data/01_proc_reads/prefixes
