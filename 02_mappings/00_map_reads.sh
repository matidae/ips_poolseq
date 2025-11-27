#!/usr/bin/env bash

#------------------------------------------------------------------------------
# Mapping reads with BWA, sorting, merging and indexing with sambamba
#
# Input: 
#   - list of processed FASTQ files 
#   - reference genome
# Output: 
#   - alignments as sorted and merged BAM file: $prefix.sort.bam
#------------------------------------------------------------------------------

in_dir="../data/01_proc_reads"
out_dir="../data/02_mappings"

$threads=40
# Input
prefixes="$in_dir/prefixes"
filelist="$in_dir/filelist"

mkdir -p "$out_dir/index"

while read -r prefix; do
    sorted_bams=()
    files=($(grep "/$prefix/" "$filelist"))
    numfiles=${#files[@]}


    for ((i=0; i<numfiles; i+=2)); do    
        r1_base=$(basename "${files[$i]}" | sed 's/.fq/.qc.fq/')
        r2_base=$(basename "${files[$i+1]}" | sed 's/.fq/.qc.fq/')


        r1="$in_dir/$prefix/$r1_base"
        r2="$in_dir/$prefix/$r2_base"

        pair_id=$((i/2 + 1))

        sam="$out_dir/$prefix.$pair_id.sam"
        bam="$out_dir/$prefix.$pair_id.bam"
        sorted="$out_dir/$prefix.sort.$pair_id.bam"
   
        #Run bwa for the first set of PE reads
        bwa-mem2 mem \
        -t "$threads" \
        -R "@RG\tID:${prefix}\tSM:${prefix}" \
        -o "$sam" \
        "$out_dir/index/Ips_typograpgus_LG16corrected.final.bwa_index" \
        "$r1" "$r2" \
        2>> "$out_dir/$prefix.bwa.log"
        
        # Sam files to bam format
        sambamba view -q -t "$threads" -S -f bam -o "$bam" "$sam"        
        # Sort bam
        sambamba sort -q -t "$threads" -o "$sorted" "$bam"
        #Array of sorted bams
        sorted_bams+=("$sorted")
                
        # Delete non sorted sam and bam files
        rm "$sam" "$bam"
    done

    merged_bams="$out_dir/$prefix.sort.bam"

    # Merge if more than one sorted BAM
    if (( ${#sorted_bams[@]} > 1 )); then
        sambamba merge -q -t "$threads" "$merged_bams" "${sorted_bams[@]}"
    else
        cp "${sorted_bams[0]}" "$merged_bams"
    fi

    # Index bam file
    sambamba index -q -t "$threads" "$merged_bams"

    # Delete sam, bam, sort and index files after merging
    for s in "${sorted_bams[@]}"; do
        rm "$s" "$s.bai"
    done      

done < "$prefixes"
