#!/usr/bin/env bash

#----------------------------------------------------------------------
# Categorizes SNPs based on FDR-corrected p-values into 4 categories: 
# drift (null model), directional, fluctuating (saturated model), and mixed.
# Thins the SNPs to one per gene to avoid linkage effects in downstream analyses.
#
# Input:
#   - ../results/08_models/tests_{prefix}_FDR.tsv
#
# Output:
#   - ../results/08_models/s1_filter/tests.{prefix}.[drift|directional|fluctuating|mixed].tsv 
#   - ../results/08_models/s1_filter/tests_summary.tsv: Summary table with counts and percentages
#   - ../results/08_models/s2_thinning/tests.{prefix}.[drift|directional|fluctuating|mixed].thinned.tsv 
#----------------------------------------------------------------------

work_dir="../results/08_models"
out_dir_filter="../results/08_models/s1_filter"
out_dir_gene="../results/08_models/s2_thinning"

bed_file="../data/reference/Ips_typograpgus_LG16corrected.liftoff.genes.bed"
summary_file="$out_dir_filter/tests_summary.tsv"

mkdir -p "$out_dir_filter"
mkdir -p "$out_dir_gene"

for f in "$work_dir"/tests.*.FDR.tsv; do
    base=$(basename "$f" .tsv)
    # Drift : LRT1 ns and LRT2 ns
    awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} ($16 >= 0.05 && $17 >= 0.05)' "$f" > "$out_dir_filter/${base}.drift.tsv"
    # Directional LRT1 sig and LRT2 ns
    awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} ($16 < 0.05 && $17 >= 0.05)' "$f" > "$out_dir_filter/${base}.directional.tsv"
    # Fluctuating : LRT1 ns and LRT2 sig
    awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} ($16 >= 0.05 && $17 < 0.05)' "$f" > "$out_dir_filter/${base}.fluctuating.tsv"
    # Mixed : LRT1 sig and LRT2 sig
    awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} ($16 < 0.05 && $17 < 0.05)' "$f" > "$out_dir_filter/${base}.mixed.tsv"
done


# Header for summary file
echo -e "prefix\ttotal\tdrift\tdirec\tfluct\tmixed\t%drift\t%direct\t%fluct\t%mixed" > "$summary_file"

# Loop over all *_FDR.tsv files
for f in "$work_dir"/tests*.FDR.tsv; do
    sample=$(basename "$f" .tsv)  # sample name without extension

    # Count total SNPs (subtract 1 for header)
    total=$(($(wc -l < "$f") - 1))

    # Count SNPs in each category from filtered files in out_dir
    drift=$(($(wc -l < "$out_dir_filter/${sample}.drift.tsv") - 1))
    direct=$(($(wc -l < "$out_dir_filter/${sample}.directional.tsv") - 1))
    fluct=$(($(wc -l < "$out_dir_filter/${sample}.fluctuating.tsv") - 1))
    mixed=$(($(wc -l < "$out_dir_filter/${sample}.mixed.tsv") - 1))

    # Calculate percentages    
    if [ "$total" -gt 0 ]; then
        drift_p=$(awk "BEGIN{printf \"%.2f\", $drift/$total*100}")
        direct_p=$(awk "BEGIN{printf \"%.2f\", $direct/$total*100}")
        fluct_p=$(awk "BEGIN{printf \"%.2f\", $fluct/$total*100}")
        mixed_p=$(awk "BEGIN{printf \"%.2f\", $mixed/$total*100}")
    else
        drift_p=0; direct_p=0; fluct_p=0; mixed_p=0
    fi

    # Append row to summary
    echo -e "$sample\t$total\t$drift\t$direct\t$fluct\t$mixed\t$drift_p\t$direct_p\t$fluct_p\t$mixed_p" >> "$summary_file"
done


files_to_analyze=$(ls "$out_dir_filter"/tests.*.tsv | egrep "SFIN|WFIN" )

#Intersect SNPs to gene coordinates and select one SNP per gene to avoid linkage effects.
for f in $files_to_analyze; do
    base=$(basename "$f" .tsv)
    awk 'NR==1{next} {print $1, $2-1, $2, $0}' OFS="\t" "$f" > "$out_dir_gene/${base}.bed"
    bedtools intersect -a "$out_dir_gene/${base}.bed" -b "$bed_file" -wa -wb | \
     awk 'BEGIN{OFS="\t"} {print $24, $4, $5, $6, $7, $19, $20}' | \
     awk 'BEGIN{OFS="\t"} !seen[$2,$3,$4]++ {print $0}' > "$out_dir_gene/${base}.annotated.tsv"
    
    # Thinning of SNps, by selecting the SNP with lowest p-value for each gene
    if [[ "$base" == *fluctuating* ]]; then
        # Lowest fluctuating p-value (column 7)
        sort -k1,1 -k7,7g "$out_dir_gene/${base}.annotated.tsv" | \
            awk 'BEGIN{OFS="\t"} !seen[$1]++' > "$out_dir_gene/${base}.thinned.tsv"
    elif [[ "$base" == *directional* ]]; then
        # Lowest directional p-value (column 6)
        sort -k1,1 -k6,6g "$out_dir_gene/${base}.annotated.tsv" | \
            awk 'BEGIN{OFS="\t"} !seen[$1]++' > "$out_dir_gene/${base}.thinned.tsv"
    elif [[ "$base" == *drift* ]]; then
        # Highest p-values for both tests (most confidently neutral)
        sort -k1,1 -k6,6r -k7,7r "$out_dir_gene/${base}.annotated.tsv" | \
            awk 'BEGIN{OFS="\t"} !seen[$1]++' > "$out_dir_gene/${base}.thinned.tsv"
     elif [[ "$base" == *mixed* ]]; then
        # Lowest directional p-value (col 6) among inconclusive SNPs
        sort -k1,1 -k6,6g -k7,7g "$out_dir_gene/${base}.annotated.tsv" | \
            awk 'BEGIN{OFS="\t"} !seen[$1]++' > "$out_dir_gene/${base}.thinned.tsv"
    fi
done
 