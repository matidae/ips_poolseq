#!/usr/bin/env bash

#----------------------------------------------------------------------
# Runs variant calling in parallel across genomic regions specified in the 'regions_1M' file
# Uses bcftools mpileup and bcftools call to generate VCF files for each region.
# Concatenates the chunks into the file ips_merged.vcf.gz
#
# Run: ./00_varcall.sh
# In: ../../data/03_dedup/*.dedup.sort.nam
# Out: ../../results/04_varcalls/region_{#}_$chunk.vcf and ips_merged.vcf.gz
#----------------------------------------------------------------------

eval "$(conda shell.bash hook)"
conda activate extra

wd="../../results/04_varcalls"
genome="../../data/reference"
bam=$(find ../../data/03_dedup/ -name "*.dedup.sort.bam" | sort | tr '\n' ' ')

#Run bcftools mpileup and call in parallel taking genome chunks defined in regions_1M
cat regions_1M | parallel -j40 '
    region={};
    chunk=$(echo $region | sed "s/[:-]/_/g");
    bcftools mpileup -Ou -I -a FORMAT/AD --max-depth 5000 \
        -r $region \
        -f '"$genome"/Ips_typograpgus_LG16corrected.final.fasta \
        '"$bam"' |
    bcftools call -vmO v -o '"$wd"/region_{#}_$chunk.vcf
'

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

mkdir "$wd/chunks"
mv "$wd/region_*"  "$wd/chunks/"