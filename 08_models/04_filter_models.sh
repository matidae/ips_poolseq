#!/usr/bin/env bash

#----------------------------------------------------------------------
# Categorizes SNPs based on FDR-corrected p-values into 4 categories: 
# drift (null model), directional, fluctuating (saturated model), and inconclusive.
#
# Input:
#   - ../results/08_models/tests_{prefix}_FDR.tsv
#
# Output:
#   - ../results/08_models/filter/tests_{prefix}_[drift|directional|fluctuating|inconclusive].tsv 
#   - ../results/08_models/filter/tests_summary.tsv: Summary table with counts and percentages
#----------------------------------------------------------------------

work_dir="../results/08_models"
out_dir_filter="../results/08_models/filter"
out_dir_gene="../results/08_models/snp_to_gene"

bed_file="../data/reference/Ips_typograpgus_LG16corrected.liftoff.genes.bed"
summary_file="$out_dir_filter/tests_summary.tsv"

mkdir -p "$out_dir_filter"
mkdir -p "$out_dir_gene"

for f in "$work_dir"/tests_*_FDR.tsv; do
    base=$(basename "$f" .tsv)
    # Drift 
    awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} ($16 >= 0.05 && $17 >= 0.05)' "$f" > "$out_dir_filter/${base}_drift.tsv"
    # Directional
    awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} ($16 < 0.05 && $17 < 0.05)' "$f" > "$out_dir_filter/${base}_directional.tsv"
    # Fluctuating
    awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} ($16 >= 0.05 && $17 < 0.05)' "$f" > "$out_dir_filter/${base}_fluctuating.tsv"
    # Inconclusive
    awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} ($16 < 0.05 && $17 >= 0.05)' "$f" > "$out_dir_filter/${base}_inconclusive.tsv"
done


# Header
echo -e "prefix\ttotal\tdrift\tdirec\tfluct\tinconc\t%drift\t%direct\t%fluct\t%inconc" > "$summary_file"

# Loop over all *_FDR.tsv files
for f in "$work_dir"/tests*_FDR.tsv; do
    sample=$(basename "$f" .tsv)  # sample name without extension

    # Count total SNPs (subtract 1 for header)
    total=$(($(wc -l < "$f") - 1))

    # Count SNPs in each category from filtered files in out_dir
    drift=$(($(wc -l < "$out_dir_filter/${sample}_drift.tsv") - 1))
    direct=$(($(wc -l < "$out_dir_filter/${sample}_directional.tsv") - 1))
    fluct=$(($(wc -l < "$out_dir_filter/${sample}_fluctuating.tsv") - 1))
    inconc=$(($(wc -l < "$out_dir_filter/${sample}_inconclusive.tsv") - 1))

    
    if [ "$total" -gt 0 ]; then
        drift_p=$(awk "BEGIN{printf \"%.2f\", $drift/$total*100}")
        direct_p=$(awk "BEGIN{printf \"%.2f\", $direct/$total*100}")
        fluct_p=$(awk "BEGIN{printf \"%.2f\", $fluct/$total*100}")
        inconc_p=$(awk "BEGIN{printf \"%.2f\", $inconc/$total*100}")
    else
        drift_p=0; direct_p=0; fluct_p=0; inconc_p=0
    fi

    # Append row to summary
    echo -e "$sample\t$total\t$drift\t$direct\t$fluct\t$inconc\t$drift_p\t$direct_p\t$fluct_p\t$inconc_p" >> "$summary_file"
done


out_dir_filter="../results/08_models/filter"
out_dir_gene="../results/08_models/snp_to_gene"
files_to_analyze=$(ls "$out_dir_filter"/tests_*.tsv | egrep "SFIN|WFIN" | egrep "directional|fluctuating|drift" )

#Intersect SNPs to gene coordinates and select one SNP per gene to avoid linkage effects.
for f in $files_to_analyze; do
    base=$(basename "$f" .tsv)
    awk 'NR==1{next} {print $1, $2-1, $2, $0}' OFS="\t" "$f" > "$out_dir_gene/${base}.bed"
    bedtools intersect -a "$out_dir_gene/${base}.bed" -b "$bed_file" -wa -wb | \
     awk 'BEGIN{OFS="\t"} {print $24, $4, $5, $6, $7, $19, $20}' | \
     awk 'BEGIN{OFS="\t"} !seen[$2,$3,$4]++ {print $0}' > "$out_dir_gene/${base}.annotated.tsv"
    
    # Select one SNP per gene
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
        sort -k1,1 -k6,6gr -k7,7gr "$out_dir_gene/${base}.annotated.tsv" | \
            awk 'BEGIN{OFS="\t"} !seen[$1]++' > "$out_dir_gene/${base}.thinned.tsv"
    fi
done
 