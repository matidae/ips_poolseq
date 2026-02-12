#!/usr/bin/env bash

#-------------------------------------------------------------------------------
# Filter the VCF file by : snps, biallelic, MQ>=20 and QUAL>=20 known genotype 
# and MAF>=0.5 && <=0.95
# Intersect VCF with genic and intergenic regions 
# Extract AD fields from each VCF file
#
# Input:
#   - ips_merged.vcf: vcf file for all samples
# Output:
#   - genic_variants.vcf.gz : filtered VCF file 
#   - genic_readcounts.tsv : read counts of filtered SNPs from genic regions
#   - intergenic_readcounts.tsv : read counts of filtered SNPs from intergenic regions
#-------------------------------------------------------------------------------

eval "$(conda shell.bash hook)"
conda activate extra
set -euo pipefail

work_dir="../results/04_varcalls"
reference="../data/reference"
gff="Ips_typograpgus_LG16corrected.liftoff.gff"
bed="Ips_typograpgus_LG16corrected.liftoff.genes.bed"

# Filter VCF: -v snps (only snps), -m2 -M2 (biallelic)
bcftools view "$work_dir/ips_merged.vcf.gz" -m2 -M2 -v snps -e 'INFO/MQ<20 || QUAL<20' \
  -Oz -o "$work_dir/ips.biallelic_q20_m20.vcf.gz" 

# Exclude all SNPs where min(p) <= (1 - min_MAF) and max(p) >= min_MAF #min_MAF=0.05
python3 ./04_varcalls/01aux_filter_MAF.py

#Exclude SNPs with missing genotype GT ./.
bcftools view "$work_dir/ips.biallelic_q20_m20.maf05.vcf.gz" -e 'GT="./."' \
-Oz -o "$work_dir/ips.biallelic_q20_m20.maf05.gt.vcf.gz"

# Make a bed file of genes to extract genic and intergenic SNPs
if [[ ! -f "$reference/$bed" ]]; then
  awk '$3=="gene"' "$reference/$gff" | cut -f1 -d';' | cut -f1,4,5,9 | sed 's/ID=//' | \
   awk 'OFS="\t" {print $1, $2-1, $3, $4}' > "$reference/$bed"
fi

# Intersect gene coordinates with VCF file to distinguish genic and intergenic SNPs
# Create header of vcf file since bedtools omits it
bcftools view -h "$work_dir/ips.biallelic_q20_m20.maf05.gt.vcf.gz" > "$work_dir/genic_variants.vcf"
bcftools view -h "$work_dir/ips.biallelic_q20_m20.maf05.gt.vcf.gz" > "$work_dir/intergenic_variants.vcf"

#Extract genic SNPs
bedtools intersect -a "$work_dir/ips.biallelic_q20_m20.maf05.gt.vcf.gz" -b "$reference/$bed" -u >> "$work_dir/genic_variants.vcf"
bgzip -f "$work_dir/genic_variants.vcf"

# Extract intergenic SNPs
bedtools intersect -a "$work_dir/ips.biallelic_q20_m20.maf05.gt.vcf.gz" -b "$reference/$bed" -v >> "$work_dir/intergenic_variants.vcf"
bgzip -f "$work_dir/intergenic_variants.vcf"

# Create header for read counts table
bcftools view -h "$work_dir/ips.biallelic_q20_m20.maf05.gt.vcf.gz" |tail -n1 | cut -f-2,4,5,10- | \
  tee "$work_dir/genic_readcounts.tsv" "$work_dir/intergenic_readcounts.tsv" > /dev/null

# Get table of read counts for genic SNPs
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD\t]\n' "$work_dir/genic_variants.vcf.gz" | \
  sed 's/\t$//' >> "$work_dir/genic_readcounts.tsv"

# Get table of read counts for intergenic SNPs
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD\t]\n' "$work_dir/intergenic_variants.vcf.gz" | \
  sed 's/\t$//' >> "$work_dir/intergenic_readcounts.tsv"
