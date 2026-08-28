#!/usr/bin/env bash

#-------------------------------------------------------------------------------
# Filter the VCF file by : snps, biallelic, MQ>=30 and QUAL>=30 known genotype 
# and MAF>=0.05 && <=0.95
# Intersect VCF with genic regions 
# Extract AD fields from each VCF file
#
# Input:
#   - $VARCALL_RESULTS/ips_merged.vcf.gz: vcf file for all samples
# Output:
#   - $VARCALL_RESULTS/genic_variants.vcf.gz : filtered VCF file 
#   - $VARCALL_RESULTS/genic_readcounts.tsv : read counts of filtered SNPs from genic regions
#-------------------------------------------------------------------------------

set -euo pipefail
source ./utils/paths.sh

work_dir="$VARCALL_RESULTS"

log "=== VCF filtering start ==="

# Identify indels and their flanking positions. Local misalignment around indels
# is a recurrent source of false SNP calls, so these regions are excluded later.
bcftools view -v indels "$work_dir/ips_merged.vcf.gz" \
  | bcftools query -f '%CHROM\t%POS\t%REF\n' \
  | awk -v p=5 'BEGIN{OFS="\t"} {s=$2-1-p; if(s<0)s=0; print $1, s, $2-1+length($3)+p}' \
  | sort -k1,1 -k2,2n | bedtools merge > "$work_dir/indel_regions.bed"
log "done: $work_dir/indel_regions.bed ($(wc -l < "$work_dir/indel_regions.bed") regions)"

# Filter VCF: -v snps (only snps), -m2 -M2 (biallelic)
bcftools view "$work_dir/ips_merged.vcf.gz" -m2 -M2 -v snps -e 'INFO/MQ<30 || QUAL<30' \
  -t ^LG16 -Oz -o "$work_dir/ips.biallelic_q30_m30.vcf.gz"
log "done: $work_dir/ips.biallelic_q30_m30.vcf.gz" 

# Keep all SNPs where min(p) <= (1 - min_MAF) and max(p) >= min_MAF #min_MAF=0.05
python3 ./04_varcalls/01aux_filter_MAF.py
log "done: $work_dir/ips.biallelic_q30_m30.maf05.vcf.gz"

# Exclude SNPs with missing genotype GT ./. or .
bcftools view "$work_dir/ips.biallelic_q30_m30.maf05.vcf.gz" -e 'GT="./." || GT="."' \
    -Oz -o "$work_dir/ips.biallelic_q30_m30.maf05.gt.vcf.gz"

# Exclude SNPs within 5 bp of a called indel. Reads spanning indels are
# frequently misaligned locally, producing spurious SNP calls.
bedtools intersect -header -v \
    -a "$work_dir/ips.biallelic_q30_m30.maf05.gt.vcf.gz" \
    -b "$work_dir/indel_regions.bed" \
  | bgzip > "$work_dir/ips.biallelic_q30_m30.maf05.gt.noindel.vcf.gz"
tabix -f -p vcf "$work_dir/ips.biallelic_q30_m30.maf05.gt.noindel.vcf.gz"


# Intersect gene coordinates with VCF file to distinguish genic SNPs (use whole gene and filtered out repeats)
declare -A beds=( [genes]="$GENES_BED" [norepeat]="$GENES_BED_NOREPEAT" )
for v in genes norepeat; do
    bed="${beds[$v]}"
    log "intersecting with genic regions: $bed"
    bcftools view -h "$work_dir/ips.biallelic_q30_m30.maf05.gt.noindel.vcf.gz" > "$work_dir/genic_variants_$v.vcf"
    bedtools intersect -a "$work_dir/ips.biallelic_q30_m30.maf05.gt.noindel.vcf.gz" \
    -b "$bed" -u >> "$work_dir/genic_variants_$v.vcf"
    bgzip -f "$work_dir/genic_variants_$v.vcf"
    tabix -f -p vcf "$work_dir/genic_variants_$v.vcf.gz"
    log "done: $work_dir/genic_variants_$v.vcf.gz"

    bcftools view -h "$work_dir/genic_variants_$v.vcf.gz" | tail -n1 | cut -f-2,4,5,10- \
        > "$work_dir/genic_readcounts_$v.tsv"
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD\t]\n' "$work_dir/genic_variants_$v.vcf.gz" | \
        sed 's/\t$//' >> "$work_dir/genic_readcounts_$v.tsv"
    log "done: $work_dir/genic_readcounts_$v.tsv"
done

log "=== VCF filtering complete ==="