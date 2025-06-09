#!/usr/bin/env bash

#----------------------------------------------------------------------
# Script to gzip, index and merge all vcf files chunks coming from parallelized bcftools
# Run: ./01_merge_vcf.sh
# Out: summary_table with stats of mapping and deduplication
#----------------------------------------------------------------------
conda init
conda activate extra

wd="../../results/04_varcalls"
# Compress and index each VCF file
for i in "$wd"/*.vcf; do
    bgzip "$i"
    tabix -p vcf "$i.gz"
done

# Sort files numerically based on the region number before concatenation
sorted_vcfs=$(ls "$wd"/region_*.vcf.gz | sort -V)

# Concatenate the sorted VCF files
bcftools concat -Oz -o "$wd/ips_merged.vcf.gz" $sorted_vcfs

# Index the merged VCF
tabix -p vcf "$wd/ips_merged.vcf.gz"
