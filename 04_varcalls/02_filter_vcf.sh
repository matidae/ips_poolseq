#!/usr/bin/env bash

#----------------------------------------------------------------------
# Fiter the vcf file keeping selected SNPs
# Run: ./02_filter_vcf.sh
# Out: Filtered vcf.gz file
#----------------------------------------------------------------------
conda init
conda activate extra

wd="../../results/04_varcalls"

# Filter VCF to keep only biallelic SNPs
bcftools view -m2 -M2 -v snps -Oz -o $wd/ips.bi.vcf.gz $wd/ips_merged.vcf.gz

# Apply MQ ≥ 20 and QUAL ≥ 20 filter
bcftools filter -i 'INFO/MQ>=20 && QUAL>=20' -Oz -o $wd/ips.bi.q20.vcf.gz $wd/ips.bi.vcf.gz

# Remove missing genotypes (./.) directly without unnecessary piping
bcftools view -e 'GT="./."' -Oz -o $wd/ips.bi.q20.full.vcf.gz $wd/ips.bi.q20.vcf.gz

