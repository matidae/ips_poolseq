#!/usr/bin/env bash

#------------------------------------------------------------------------------
# Analyze filtering steps, plot filtered SNPs per sample
#
# Input:
#   - ips_merged.vcf: vcf file for all samples
#
# Output: 
#   - stats_filter.tsv: file with number of SNPs per sample filtered by quality and genotype
#-------------------------------------------------------------------------------

work_dir="../../results/04_varcalls"

total_snps=$(bcftools view "$work_dir/ips_merged.vcf.gz" -m2 -M2 -v snps | grep -v '^#' | wc -l)

bcftools view "$work_dir/ips_merged.vcf.gz" -m2 -M2 -v snps -i 'INFO/MQ<20 || QUAL<20' \
  | bcftools query -f '[%SAMPLE\t%GT\n]' | cut -f1 | sort | uniq -c \
  | sort -k1,1nr | awk '{print $1"\t"$2}' | sort -k2 \
  > "$work_dir/qual_filter.tsv"

bcftools view "$work_dir/ips_merged.vcf.gz" -m2 -M2 -v snps -e 'INFO/MQ<20 || QUAL<20' \
  | bcftools query -f '[%SAMPLE\t%GT\n]' | awk '$2 == "./."' | cut -f1 | sort | uniq -c \
  | sort -k1,1nr | awk '{print $1"\t"$2}' | sort -k2 \
  > "$work_dir/genotype_filter.tsv"

paste "$work_dir/qual_filter.tsv" "$work_dir/genotype_filter.tsv" \
| awk -v total="$total_snps" '{print $2, $1, $3, $1/total, $3/total}' > "$work_dir/summary_filter.tsv"
