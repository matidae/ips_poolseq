#----------------------------------------------------------------------
# Runs variant calling in parallel across genomic regions specified in the 'regions_1M' file
# Uses bcftools mpileup and bcftools call to generate VCF files for each region.
# Run: ./00_varcall.sh
# In: ../../data/03_dedup/*.dedup.sort.nam
# Out: ../../results/04_varcalls/region_{#}_$chunk.vcf
#----------------------------------------------------------------------

#!/usr/bin/env bash

BAM_FILES=$(find ../../data/03_dedup/ -name "*.dedup.sort.bam" | sort | tr '\n' ' ')

cat regions_1M | parallel -j40 '
    region={};
    chunk=$(echo $region | sed "s/[:-]/_/g");
    bcftools mpileup -Ou -I -a FORMAT/AD --max-depth 5000 \
        -r $region \
        -f ../../data/reference/Ips_typograpgus_LG16corrected.final.fasta \
        '"$BAM_FILES"' |
    bcftools call -vmO v -o ../../results/04_varcalls/region_{#}_$chunk.vcf
'
