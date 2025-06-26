#!/usr/bin/env bash

#-------------------------------------------------------------------------------
# Filter the VCF file by : snps, biallelic, MQ>=20 and QUAL>=20 known genotype and MAF>=0.5
# Intersect VCF with genic and intergenic regions 
# Extract AD fields from each VCF file
#
# Input:
#   - ips_merged.vcf: vcf file for all samples
# Output: 
#   - genic_readcounts.tsv : read counts of filtered SNPs from genic regions
#   - intergenic_readcounts.tsv : read counts of filtered SNPs from intergenic regions
#-------------------------------------------------------------------------------

eval "$(conda shell.bash hook)"
conda activate extra

work_dir="../../results/04_varcalls"
reference="../../data/reference"
gff="Ips_typograpgus_LG16corrected.liftoff.gff"
bed="Ips_typograpgus_LG16corrected.liftoff.genes.bed"

# Filter VCF: -v snps (only snps), -m2 -M2 (biallelic) -e 'GT="./."' (exclude missing GT)
bcftools view "$work_dir/ips_merged.vcf.gz" -m2 -M2 -v snps -e 'GT="./." || INFO/MQ<20 || QUAL<20' \
  -Oz -o "$work_dir/ips.filter.vcf.gz" 

# Exclude all SNPs where MAF < 0.05
bcftools filter "$work_dir/ips.filter.vcf.gz" -e 'MAF[0] < 0.05' -Oz -o "$work_dir/ips.filter.m05.vcf.gz"

# Make a bed file of genes to extract genic and intergenic SNPs
awk 'OFS="\t"{ if ($3=="gene") print $1, $4, $5 }' "$reference/$gff" > "$reference/$bed"

# Intersect gene coordinates with VCF file to distinguish genic and intergenic SNPs
# Create header of vcf file since bedtools omits it
bcftools view -h "$work_dir/ips.filter.m05.vcf.gz" | tee "$work_dir/genic_variants.vcf" "$work_dir/intergenic_variants.vcf" > /dev/null

#Extract genic SNPs
bedtools intersect -a "$work_dir/ips.filter.m05.vcf.gz" -b "$reference/$bed" -u >> "$work_dir/genic_variants.vcf"
bgzip "$work_dir/genic_variants.vcf"

# Extract intergenic SNPs
bedtools intersect -a "$work_dir/ips.filter.m05.vcf.gz" -b "$reference/$bed" -v >> "$work_dir/intergenic_variants.vcf"
bgzip  "$work_dir/intergenic_variants.vcf"

# Create header for read counts table
bcftools view -h "$work_dir/ips.filter.m05.vcf.gz" |tail -n1 | cut -f-2,4,5,10- | \
  tee "$work_dir/genic_readcounts.tsv" "$work_dir/intergenic_readcounts.tsv" > /dev/null

# Get table of read counts for genic SNPs
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD\t]\n' "$work_dir/genic_variants.vcf.gz" | \
  sed 's/\t$//' >> "$work_dir/genic_readcounts.tsv"

# Get table of read counts for intergenic SNPs
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD\t]\n' "$work_dir/intergenic_variants.vcf.gz" | \
  sed 's/\t$//' >> "$work_dir/intergenic_readcounts.tsv"
