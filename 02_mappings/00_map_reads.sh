#!/usr/bin/env bash

#------------------------------------------------------------------------------
# Mapping reads with BWA, sorting, merging and indexing with sambamba
#
# Input: 
#   - list of QC reads and genome for mapping
# Output: 
#   - merged and sorted alignment as {prefix}.sort.bam files
#------------------------------------------------------------------------------


while read -r prefix; do
    files=($(grep "$prefix" ./data/01_proc_reads/filelist))
    numfiles=${#files[@]}


    for ((i=0; i<numfiles; i+=2)); do    
        r1_base=$(echo ${files[$i]} | cut -f6 -d'/' | sed 's/.fq/.qc.fq/')
        r2_base=$(echo ${files[$i+1]} | cut -f6 -d'/' | sed 's/.fq/.qc.fq/')


        r1="./data/01_proc_reads/$prefix/$r1_base"
        r2="./data/01_proc_reads/$prefix/$r2_base"

        pair_id=$((i/2 + 1))

        sam="./data/02_mappings/${prefix}.${pair_id}.sam"
        bam="./data/02_mappings/${prefix}.${pair_id}.bam"
        sorted="./data/02_mappings/${prefix}.sort.${pair_id}.bam"
   
        #Run bwa for the first set of PE reads
        bwa-mem2 mem \
        -t 40 \
        -R "@RG\tID:${prefix}\tSM:${prefix}" \
        -o "$sam" \
        ./data/02_mappings/index/Ips_typograpgus_LG16corrected.final.bwa_index \
        "$r1" "$r2" \
        2>>./data/02_mappings/"$prefix".bwa.log
        
        # Sam files to bam format
        sambamba view -q -t 40 -S -f bam -o "$bam" "$sam"        
        # Sort bam
        sambamba sort -q -t 40 -o "$sorted" "$bam"
        #Array of sorted bams
        sorted_bams+=("$sorted")
                
        # Delete non sorted sam and bam files
        rm "$sam" "$bam"
    done

    merged_bams="./data/02_mappings/${prefix}.sort.bam"

    # Merge if more than one sorted BAM
    if (( ${#sorted_bams[@]} > 1 )); then
        sambamba merge -q -t 40 "$merged_bams" "${sorted_bams[@]}"
    else
        cp "${sorted_bams[0]}" "$merged_bams"
    fi

    # Index bam file
    sambamba index -q -t 40 "$merged_bam"

    # Delete sam, bam, sort and index files before merging
    for s in "${sorted_bams[@]}"; do
        rm "$s" "$s.bai"
    done      

done < data/01_proc_reads/prefixes
