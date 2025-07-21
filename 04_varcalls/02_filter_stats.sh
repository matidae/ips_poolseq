#!/usr/bin/env bash

#------------------------------------------------------------------------------
# Analyze filtering steps, plot filtered SNPs per sample
#
# Input:
#   - ips_merged.vcf: vcf file for all samples
# Output: 
#   - missing_GT_stats.tsv: file with number and proportion of SNPs per sample with missing genotype
#-------------------------------------------------------------------------------

work_dir="../../results/04_varcalls"

total_snps=$(bcftools view "$work_dir/ips_merged.vcf.gz" -m2 -M2 -v snps -e 'INFO/MQ<20 || QUAL<20' \
| bcftools filter -e 'MAF[0] < 0.05' | grep -v '^#' | wc -l)

bcftools view "$work_dir/ips_merged.vcf.gz" -m2 -M2 -v snps -e 'INFO/MQ<20 || QUAL<20' \
| bcftools filter -e 'MAF[0] < 0.05' | bcftools query -f '[%SAMPLE\t%GT\n]' \
| awk '$2 == "./."' | cut -f1 | sort | uniq -c | awk '{print $1"\t"$2}' | sort -k2 \
> "$work_dir/missing_GT.tsv"

awk -v total="$total_snps" '{printf "%s %s %.4f\n", $2, $1, $1/total}' "$work_dir/missing_GT.tsv" \
> "$work_dir/missing_GT_stats.tsv"

rm "$work_dir/missing_GT.tsv"
