#!/usr/bin/env bash

#--------------------------------------------------------------------------------
# Runs variant calling in parallel across genomic regions specified regions_1M.txt
# Uses bcftools mpileup and bcftools call to generate VCF files for each region.
# Concatenates the chunks into the file ips_merged.vcf.gz
#
# Input: 
#   - ../../data/03_dedup/*.dedup.sort.bam : mappings of all samples
# Output: 
#   - ips_merged.vcf.gz : vcf file for all samples
#-------------------------------------------------------------------------------

eval "$(conda shell.bash hook)"
conda activate extra

work_dir="../../results/04_varcalls"
genome="../../data/reference"
bam=$(find ../../data/03_dedup/ -name "*.dedup.sort.bam" | egrep -v SFIN_La_2016 | sort | tr '\n' ' ')

bedtools makewindows -g "$genome/Ips_typograpgus_LG16corrected.final.fasta.fai" -w 1000000 \
    | awk '{print $1":"$2+1"-"$3}' > "$work_dir/regions_1M.txt"

#Run bcftools mpileup and call in parallel taking genome chunks defined in regions_1M
cat "$work_dir/regions_1M.txt" | parallel -j40 "
    region={};
    chunk=\$(echo \$region | sed 's/[:-]/_/g');
    bcftools mpileup -Ou -I -a FORMAT/AD --max-depth 5000 \
        -r \$region \
        -f \"$genome/Ips_typograpgus_LG16corrected.final.fasta\" \
        $bam |
    bcftools call -vmO v -o \"$work_dir/region_{#}_\${chunk}.vcf\"
"
# Compress and index each VCF file
for i in "$work_dir"/*.vcf; do
    bgzip "$i"
    tabix -p vcf "$i.gz"
done

# Sort files numerically based on the region number before concatenation
sorted_vcfs=$(ls "$work_dir"/region_*.vcf.gz | sort -V)

# Concatenate the sorted VCF files
bcftools concat -Oz -o "$work_dir/ips_merged.vcf.gz" $sorted_vcfs

# Index the merged VCF
tabix -p vcf "$work_dir/ips_merged.vcf.gz"

mkdir "$work_dir/chunks"
mv "$work_dir/region_*"  "$work_dir/chunks/"