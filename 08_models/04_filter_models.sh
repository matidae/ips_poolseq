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
out_dir="../results/08_models/filter"
summary_file="$out_dir/tests_summary.tsv"

mkdir -p "$out_dir"
for f in "$work_dir"/tests_*_FDR.tsv; do
    base=$(basename "$f" .tsv)
    # Drift 
    awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} ($16 >= 0.05 && $17 >= 0.05)' "$f" > "$out_dir/${base}_drift.tsv"
    # Directional
    awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} ($16 < 0.05 && $17 < 0.05)' "$f" > "$out_dir/${base}_directional.tsv"
    # Fluctuating
    awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} ($16 >= 0.05 && $17 < 0.05)' "$f" > "$out_dir/${base}_fluctuating.tsv"
    # Inconclusive
    awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} ($16 < 0.05 && $17 >= 0.05)' "$f" > "$out_dir/${base}_inconclusive.tsv"
done


# Header
echo -e "prefix\ttotal\tdrift\tdirec\tfluct\tinconc\t%drift\t%direct\t%fluct\t%inconc" > "$summary_file"

# Loop over all *_FDR.tsv files
for f in "$work_dir"/tests*_FDR.tsv; do
    sample=$(basename "$f" .tsv)  # sample name without extension

    # Count total SNPs (subtract 1 for header)
    total=$(($(wc -l < "$f") - 1))

    # Count SNPs in each category from filtered files in out_dir
    drift=$(($(wc -l < "$out_dir/${sample}_drift.tsv") - 1))
    direct=$(($(wc -l < "$out_dir/${sample}_directional.tsv") - 1))
    fluct=$(($(wc -l < "$out_dir/${sample}_fluctuating.tsv") - 1))
    inconc=$(($(wc -l < "$out_dir/${sample}_inconclusive.tsv") - 1))

    
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
