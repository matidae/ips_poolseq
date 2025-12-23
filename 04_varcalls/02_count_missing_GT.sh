#!/usr/bin/env bash

#------------------------------------------------------------------------------
# Count number and proportion of SNPs with missing genotype per sample
#
# Input:
#   - ips_merged.vcf: vcf file for all samples
# Output: 
#   - missing_GT_stats.tsv: file with number and proportion of SNPs per sample with missing genotype
#-------------------------------------------------------------------------------

work_dir="../results/04_varcalls"

total_snps=$(bcftools view "$work_dir/ips.biallelic_q20_m20.maf05.vcf.gz" | grep -v '^#' | wc -l)

bcftools view "$work_dir/ips.biallelic_q20_m20.maf05.vcf.gz" | bcftools query -f '[%SAMPLE\t%GT\n]' \
| awk '$2 == "./."' | cut -f1 | sort | uniq -c | awk '{print $1"\t"$2}' | sort -k2 \
> "$work_dir/missing_GT.tsv"

awk -v total="$total_snps" '{printf "%s %s %.4f\n", $2, $1, $1/total}' "$work_dir/missing_GT.tsv" \
> "$work_dir/missing_GT_stats.tsv"

rm "$work_dir/missing_GT.tsv"
