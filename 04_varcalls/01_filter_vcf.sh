#!/usr/bin/env bash

#----------------------------------------------------------------------
# Filter the vcf file by : snps, biallelic, MQ>=20 and QUAL>=20 known genotype and MAF>=0.5
# Intersects using bedtools with genic and intergenic regions creating two vcf files
# Get readcount of each of this vcf files to obtain depth per SNP
#
# Run: ./02_filter_vcf.sh
# In: ips_merged.vcf
# Out: Filtered all SNPs vcf.gz file, and alos genic and intergenic SNPs .vcf.gz and readcounts
#----------------------------------------------------------------------

eval "$(conda shell.bash hook)"
conda activate extra

wd="../../results/04_varcalls"
reference="../../data/reference"
gff="Ips_typograpgus_LG16corrected.liftoff.gff"
genes="Ips_typograpgus_LG16corrected.liftoff.genes.bed"

# Filter VCF: -v snps (only snps), -m2 -M2 (biallelic) -e 'GT="./."' (exclude missing GT)
bcftools view "$wd/ips_merged.vcf.gz" -m2 -M2 -v snps -e 'GT="./." || INFO/MQ<20 || QUAL<20' \
  -Oz -o "$wd/ips.filter.vcf.gz" 

# Exclude all SNPs where MAF < 0.05
bcftools filter "$wd/ips.filter.vcf.gz" -e 'MAF[0] < 0.05' -Oz -o "$wd/ips.filter.m05.vcf.gz"

# Extract indels and save their positions
bcftools view -v indels $wd/ips_merged.vcf.gz | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' > "$wd/ips_indels.tsv"

# Make a bed file of genes to extract genic SNPs
awk 'OFS="\t"{ if ($3=="gene") print $1, $4, $5 }' "$reference/$gff" > "$reference/$bed"

# Intersect gene coordinates with VCF file to distinguish genic and intergenic SNPs
# Create header of vcf file since bedtools omits it
bcftools view -h "$wd/ips.filter.m05.vcf.gz" | tee "$wd/genic_variants.vcf" "$wd/intergenic_variants.vcf" > /dev/null

#Extract genic SNPs
bedtools intersect -a "$wd/ips.filter.m05.vcf.gz" -b "$reference/$bed" -u >> "$wd/genic_variants.vcf"
bgzip "$wd/genic_variants.vcf"

# Extract intergenic SNPs
bedtools intersect -a "$wd/ips.filter.m05.vcf.gz" -b "$reference/$bed" -v >> "$wd/intergenic_variants.vcf"
bgzip  "$wd/intergenic_variants.vcf"

# Create header of readcounts table
bcftools view -h "$wd/ips.filter.m05.vcf.gz" |tail -n1 | cut -f-2,4,5,10- | \
  tee "$wd/genic_readcounts.tsv" "$wd/intergenic_readcounts.tsv" > /dev/null

# Get table of read counts of genic SNPs
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD\t]\n' "$wd/genic_variants.vcf.gz" | \
  sed 's/\t$//' >> "$wd/genic_readcounts.tsv"

# Get table of read counts of intergenic SNPs
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD\t]\n' "$wd/intergenic_variants.vcf.gz" | \
  sed 's/\t$//' >> "$wd/intergenic_readcounts.tsv"
