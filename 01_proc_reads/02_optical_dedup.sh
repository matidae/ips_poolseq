#!/usr/bin/env bash

#------------------------------------------------------------------------------
# Clumpify optical duplicate removal for Illumina NovaSeq X PE reads
#
# Input:
#   - fastp-processed FASTQ files (*.qc.fq.gz)
#   - prefix list and filelist
#
# Output:
#   - clumpified FASTQ files (*.clumped.fq.gz)#
#------------------------------------------------------------------------------

in_dir="../data/01_proc_reads"
out_dir="../data/01_proc_reads"   

prefixes="$in_dir/prefixes"
filelist="$in_dir/filelist"

threads=40

while read -r prefix; do

    echo "Processing $prefix..."

    # Extract the fastp-cleaned R1/R2 files
    qc_files=($(grep "/$prefix/" "$filelist" | sed 's/.fq.gz/.qc.fq.gz/'))
    numfiles=${#qc_files[@]}

    for ((i=0; i<numfiles; i+=2)); do
        
        r1="${qc_files[$i]}"
        r2="${qc_files[$i+1]}"

        r1_base=$(basename "$r1")
        r2_base=$(basename "$r2")

        # Output names
        out_r1="$in_dir/$prefix/${r1_base/.qc.fq.gz/.clumped.fq.gz}"
        out_r2="$in_dir/$prefix/${r2_base/.qc.fq.gz/.clumped.fq.gz}"

        echo "  Clumpifying:"
        echo "    $r1"
        echo "    $r2"

        clumpify.sh \
            in="$r1" \
            in2="$r2" \
            out="$out_r1" \
            out2="$out_r2" \
            dedupe \
            optical \            
            subs=0 \
            dupedist=25000 \
            overwrite \
            threads="$threads"

    done

done < "$prefixes"
