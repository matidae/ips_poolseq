#!/usr/bin/env bash

#--------------------------------------------------------------------------------
# Runs variant calling in parallel across genomic regions specified regions_5M.txt
# Uses bcftools mpileup and bcftools call to generate VCF files for each region.
# Concatenates the chunks into the file ips_merged.vcf.gz
#
# Input: 
#   - $DEDUP_DIR/overlap_clip/*.dedup.overlapclip.sort.bam : mappings of all samples
# Output: 
#   - $VARCALL_RESULTS/ips_merged.vcf.gz : vcf file for all samples
#-------------------------------------------------------------------------------

set -euo pipefail
source ./utils/paths.sh

work_dir="$VARCALL_RESULTS"
dedup_dir="$DEDUP_DIR/overlap_clip"
jobs=20

#Input
#Genome file for varcalling
genome="$GENOME"
#List of samples to exclude from varcalling
exclude_file="$EXCLUDE_FILE"

log "=== Variant calling start ==="
#Create work directory if needed
mkdir -p "$work_dir"

#Index the reference genome if not already indexed
if [[ ! -f "$genome.fai" ]]; then
    log "indexing genome: $genome"
    samtools faidx "$genome"
fi

#List of dedup sorted bam files for analysis (optional removal of ones listed in exclude_file)
if [[ -s "$exclude_file" ]]; then
    bam=$(find "$dedup_dir/" -maxdepth 1 -name "*.filt.sort.overlapclip.bam" | grep -vFf "$exclude_file" | sort | tr '\n' ' ')
else
    bam=$(find "$dedup_dir/" -maxdepth 1 -name "*.filt.sort.overlapclip.bam" | sort | tr '\n' ' ')
fi

n_bams=$(echo "$bam" | tr ' ' '\n' | grep -c ".bam" || true)
log "BAMs to process: $n_bams"

# Make windows of the genome in regions of 5M for parallel run of varcalling 
bedtools makewindows -g "$genome.fai" -w 5000000 \
    | awk '{print $1":"$2+1"-"$3}' > "$work_dir/regions_5M.txt"

n_regions=$(wc -l < "$work_dir/regions_5M.txt")
log "genome split into $n_regions regions (5M windows)"

#Run bcftools mpileup and call in parallel taking genome chunks defined in regions_5M
log "running bcftools mpileup | call across $n_regions regions (jobs=$jobs)"
parallel -j $jobs --halt soon,fail=1 "
    region={};
    chunk=\$(echo \$region | sed 's/[:-]/_/g');
    bcftools mpileup -Ou -a FORMAT/AD --max-depth 500 -q 30 -Q 30 \
        -r \$region \
        -f \"$genome\" \
        $bam |
    bcftools call -vmO v -o \"$work_dir/region_{#}_\${chunk}.vcf\"
" < "$work_dir/regions_5M.txt"
log "done: variant calling complete"

# Compress and index each VCF file
for i in "$work_dir"/*.vcf; do
    bgzip "$i"
    tabix -p vcf "$i.gz"
done
log "done: bgzip + tabix"

# Sort files numerically based on the region number before concatenation
ls "$work_dir"/region_*.vcf.gz | sort -V > "$work_dir/vcf_list.txt"

# Concatenate the sorted VCF files
log "concatenating $(wc -l < "$work_dir/vcf_list.txt") chunks to ips_merged.vcf.gz"

bcftools concat -f "$work_dir/vcf_list.txt" -Oz -o "$work_dir/ips_merged.vcf.gz"

# Index the merged VCF
tabix -p vcf "$work_dir/ips_merged.vcf.gz"
log "done: ips_merged.vcf.gz indexed"

# Validate before moving  chunks
log "running bcftools stats"
bcftools stats "$work_dir/ips_merged.vcf.gz" > "$work_dir/ips_merged.stats.txt" \
    && mkdir -p "$work_dir/chunks" \
    && mv "$work_dir"/region_* "$work_dir/chunks/"

log "done: chunks moved to $work_dir/chunks/"
log "=== Variant calling complete ==="