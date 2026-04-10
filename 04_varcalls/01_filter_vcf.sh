#!/usr/bin/env bash

#-------------------------------------------------------------------------------
# Filter the VCF file by : snps, biallelic, MQ>=40 and QUAL>=40 known genotype 
# and MAF>=0.05 && <=0.95
# Intersect VCF with genic regions 
# Extract AD fields from each VCF file
#
# Input:
#   - ips_merged.vcf: vcf file for all samples
# Output:
#   - genic_variants.vcf.gz : filtered VCF file 
#   - genic_readcounts.tsv : read counts of filtered SNPs from genic regions
#-------------------------------------------------------------------------------

set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate extra

work_dir="../results/04_varcalls"
reference="../data/reference"
bed="Ips_typographus_LG16corrected.liftoff.genes.LGs.noTEs.bed"

# Filter VCF: -v snps (only snps), -m2 -M2 (biallelic)
bcftools view "$work_dir/ips_merged.vcf.gz" -m2 -M2 -v snps -e 'INFO/MQ<40 || QUAL<40' \
  -Oz -o "$work_dir/ips.biallelic_q40_m40.vcf.gz"

# Exclude all SNPs where min(p) <= (1 - min_MAF) and max(p) >= min_MAF #min_MAF=0.05
python3 ./04_varcalls/01aux_filter_MAF.py

# Exclude SNPs with missing genotype GT ./. or .
bcftools view "$work_dir/ips.biallelic_q40_m40.maf05.vcf.gz" -e 'GT="./." || GT="."' \
    -Oz -o "$work_dir/ips.biallelic_q40_m40.maf05.gt.vcf.gz"

# Intersect gene coordinates with VCF file to distinguish genic SNPs
bcftools view -h "$work_dir/ips.biallelic_q40_m40.maf05.gt.vcf.gz" > "$work_dir/genic_variants.vcf"
bedtools intersect -a "$work_dir/ips.biallelic_q40_m40.maf05.gt.vcf.gz" \
    -b "$reference/$bed" -u >> "$work_dir/genic_variants.vcf"
bgzip -f "$work_dir/genic_variants.vcf"
tabix -f -p vcf "$work_dir/genic_variants.vcf.gz"

# Create header for read counts table
bcftools view -h "$work_dir/genic_variants.vcf.gz" | tail -n1 | cut -f-2,4,5,10- \
    > "$work_dir/genic_readcounts.tsv"

# Get table of read counts for genic SNPs
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD\t]\n' "$work_dir/genic_variants.vcf.gz" | \
    sed 's/\t$//' >> "$work_dir/genic_readcounts.tsv"