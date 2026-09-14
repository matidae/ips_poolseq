#!/usr/bin/env bash

#------------------------------------------------------------------------------
# Count number and proportion of SNPs with missing genotype per sample
#
# Input:
#   - $VARCALL_RESULTS/ips.biallelic_q40_m40.maf05.vcf.gz
# Output: 
#   - $VARCALL_RESULTS/missing_GT_stats.tsv: number and proportion of SNPs per sample with missing genotype
#-------------------------------------------------------------------------------

set -euo pipefail
source ./utils/paths.sh

work_dir="$VARCALL_RESULTS"

log "=== Missing genotype stats start ==="

total_snps=$(bcftools view "$work_dir/ips.biallelic_q40_m40.maf05.vcf.gz" | grep -v '^#' | wc -l)
log "total SNPs: $total_snps"

bcftools query -f '[%SAMPLE\t%GT\n]' "$work_dir/ips.biallelic_q40_m40.maf05.vcf.gz" \
    | awk '$2 == "./."' | cut -f1 | sort | uniq -c | awk '{print $1"\t"$2}' | sort -k2 \
    > "$work_dir/missing_GT.tsv"

awk -v total="$total_snps" 'BEGIN{OFS="\t"} {printf "%s\t%s\t%.4f\n", $2, $1, $1/total}' "$work_dir/missing_GT.tsv" \
    > "$work_dir/missing_GT_stats.tsv"

rm "$work_dir/missing_GT.tsv"
log "done: $work_dir/missing_GT_stats.tsv"

log "=== Missing genotype stats complete ==="