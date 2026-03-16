#!/usr/bin/env bash

#----------------------------------------------------------------------
# Counts SNPs per selection category falling in inversions vs collinear regions of the genome.
#
# Inputs:
#   - {filter_dir}/tests.{prefix}.FDR_{category}.tsv
#   - {thinned_dir}/tests.{prefix}.FDR_{category}.thinned.tsv
#   - {inversion_bed}: BED file with inversion coordinates
#
# Output:
#   - {out_dir}/inversion_counts_summary.tsv  (s1_filter results)
#   - {out_dir}/inversion_counts_thinned.tsv  (s2_thinning results)
#----------------------------------------------------------------------

filter_dir="../results/08_models/s1_filter"
thinned_dir="../results/08_models/s2_thinning"
out_dir="../results/08_models/inversions"
inversion_bed="../data/reference/inversions.bed"
tmp_dir="${out_dir}/tmp"

populations=("SFIN_E" "SFIN_L" "WFIN_E" "WFIN_L")
categories=("drift" "directional" "fluctuating" "mixed")

mkdir -p "$out_dir" "$tmp_dir"

header="prefix\tcategory\tsubset\tn_total\tn_inv\tn_col\tpct_inv\tpct_col"

count_intersect() {
    local bed="$1" prefix="$2" category="$3" label="$4" out="$5"
    local n_inv n_col total pct_inv pct_col

    n_inv=$(bedtools intersect -a "$bed" -b "$inversion_bed" -wa -wb | wc -l)
    n_col=$(bedtools intersect -a "$bed" -b "$inversion_bed" -v | wc -l)
    total=$(( n_inv + n_col ))

    pct_inv=$(awk "BEGIN{printf \"%.1f\", ${n_inv}/${total}*100}")
    pct_col=$(awk "BEGIN{printf \"%.1f\", ${n_col}/${total}*100}")
    echo -e "${prefix}\t${category}\t${label}\t${total}\t${n_inv}\t${n_col}\t${pct_inv}\t${pct_col}" >> "$out"
}

to_bed() {
    local file="$1" chrom_col="$2" pos_col="$3" skip_header="$4"
    awk -v c="$chrom_col" -v p="$pos_col" -v skip="$skip_header" \
        'BEGIN{OFS="\t"} NR==1 && skip{next} {print $c, $p-1, $p, $0}' "$file"
}

top10() {
    local file="$1" sort_col="$2" reverse="$3"
    local n
    n=$(( $(wc -l < "$file") / 10 ))
    if [[ "$reverse" == "1" ]]; then
        sort -k${sort_col},${sort_col}gr "$file" | head -n "$n"
    else
        sort -k${sort_col},${sort_col}g  "$file" | head -n "$n"
    fi
}

# Processing s1_filter results
out_filter="${out_dir}/inversion_counts_summary.tsv"
echo -e "$header" > "$out_filter"

for prefix in "${populations[@]}"; do
    for category in "${categories[@]}"; do
        input="$filter_dir/tests.${prefix}.FDR.${category}.tsv"

        # sort col: directional=18, others=19; drift=descending, others=ascending
        # +3 offset because to_bed merges chrom/start/end columns
        sort_col=22; [[ "$category" == "directional" ]] && sort_col=21
        reverse=0;   [[ "$category" == "drift" ]] && reverse=1

        all="$tmp_dir/${prefix}_${category}_all.bed"
        top="$tmp_dir/${prefix}_${category}_top10.bed"

        to_bed "$input" 1 2 1 > "$all"
        top10  "$all" "$sort_col" "$reverse" > "$top"

        count_intersect "$all" "$prefix" "$category" "all"      "$out_filter"
        count_intersect "$top" "$prefix" "$category" "top10pct" "$out_filter"
    done
done

# Processing s2_thinning results
out_thinned="${out_dir}/inversion_counts_thinned.tsv"
echo -e "$header" > "$out_thinned"

for prefix in "${populations[@]}"; do
    for category in "${categories[@]}"; do
        input="$thinned_dir/tests.${prefix}.FDR.${category}.thinned.tsv"

        reverse=0; [[ "$category" == "drift" ]] && reverse=1

        all="$tmp_dir/${prefix}_${category}_thinned_all.bed"
        top="$tmp_dir/${prefix}_${category}_thinned_top10.bed"

        to_bed "$input" 2 3 0 > "$all"
        top10  "$all" 10 "$reverse" > "$top"

        count_intersect "$all" "$prefix" "$category" "all" "$out_thinned"
        count_intersect "$top" "$prefix" "$category" "top10pct" "$out_thinned"
    done
done

rm -rf "$tmp_dir"