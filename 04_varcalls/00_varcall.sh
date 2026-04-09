#!/usr/bin/env bash

#--------------------------------------------------------------------------------
# Runs variant calling in parallel across genomic regions specified regions_1M.txt
# Uses bcftools mpileup and bcftools call to generate VCF files for each region.
# Concatenates the chunks into the file ips_merged.vcf.gz
#
# Input: 
#   - ../data/03_dedup/*.dedup.sort.bam : mappings of all samples
# Output: 
#   - ips_merged.vcf.gz : vcf file for all samples
#-------------------------------------------------------------------------------
set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate extra

work_dir="../results/04_varcalls"
dedup_dir="../data/03_dedup"
jobs=20

#Input
#Genome file for varcalling
genome="../data/reference/Ips_typographus_LG16corrected.final.LG.fasta"
#List of samples to exclude from varcalling
exclude_file="../data/reference/exclude_samples.txt"

#Index the reference genome if not already indexed
if [[ ! -f "$genome.fai" ]]; then
    samtools faidx "$genome"
fi

#List of dedup sorted bam files for analysis (optional removal of ones listed in exclude_file)
if [[ -s "$exclude_file" ]]; then
    bam=$(find "$dedup_dir/" -maxdepth 1 -name "*.dedup.sort.bam" | grep -v -f "$exclude_file" | sort | tr '\n' ' ')
else
    bam=$(find "$dedup_dir/" -maxdepth 1 -name "*.dedup.sort.bam" | sort | tr '\n' ' ')
fi

# Make windows of the genome in regions of 5M for parallel run of varcalling 
bedtools makewindows -g "$genome.fai" -w 5000000 \
    | awk '{print $1":"$2+1"-"$3}' > "$work_dir/regions_5M.txt"

#Run bcftools mpileup and call in parallel taking genome chunks defined in regions_5M
parallel -j $jobs --halt soon,fail=1 "
    region={};
    chunk=\$(echo \$region | sed 's/[:-]/_/g');
    bcftools mpileup -Ou -I -a FORMAT/AD --max-depth 1000 -q 40 -Q 30 \
        -r \$region \
        -f \"$genome\" \
        $bam |
    bcftools call -vmO v -o \"$work_dir/region_{#}_\${chunk}.vcf\"
" < "$work_dir/regions_5M.txt"

# Compress and index each VCF file
for i in "$work_dir"/*.vcf; do
    bgzip "$i"
    tabix -p vcf "$i.gz"
done


# Sort files numerically based on the region number before concatenation
ls "$work_dir"/region_*.vcf.gz | sort -V > "$work_dir/vcf_list.txt"

# Concatenate the sorted VCF files
bcftools concat -f "$work_dir/vcf_list.txt"  -Oz -o "$work_dir/ips_merged.vcf.gz"

# Index the merged VCF
tabix -p vcf "$work_dir/ips_merged.vcf.gz"

# Validate before moving  chunks
bcftools stats "$work_dir/ips_merged.vcf.gz" > "$work_dir/ips_merged.stats.txt" \
    && mkdir -p "$work_dir/chunks" \
    && mv "$work_dir"/region_* "$work_dir/chunks/"