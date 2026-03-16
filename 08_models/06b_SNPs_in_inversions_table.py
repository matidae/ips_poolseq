#!/usr/bin/env python3

#----------------------------------------------------------------------
# Generates an HTML report with 4 tables of SNP counts and percentages in inversions vs. collinear regions.
#
# Inputs:
#   - {out_dir}/inversion_counts_summary.tsv  (s1_filter)
#   - {out_dir}/inversion_counts_thinned.tsv  (s2_thinning)
#
# Output:
#   - {out_dir}/inversion_counts_report.html
#----------------------------------------------------------------------

import os
import pandas as pd

out_dir      = "../results/08_models/inversions"
filter_file  = f"{out_dir}/inversion_counts_summary.tsv"
thinned_file = f"{out_dir}/inversion_counts_thinned.tsv"
output_file  = f"{out_dir}/SNPs_inversion_enrichment.html"

cat_colors = {
    "drift":        "#7A7A7A",
    "directional":  "#3A7DBF",
    "fluctuating":  "#E07B39",
    "mixed":        "#6A0DAD",
}

def build_table(df_subset, title, subtitle):
    prefixes   = df_subset["prefix"].unique()
    categories = ["drift", "directional", "fluctuating", "mixed"]

    header_cells = "<th>Category</th>" + "".join(
        f"<th>{p}</th>" for p in prefixes
    )
    subheader_cells = "<th></th>" + "".join(
        "<th>n &nbsp;|&nbsp; %inv</th>" for _ in prefixes
    )

    rows_html = ""
    for cat in categories:
        color = cat_colors.get(cat, "#333")
        cells = f"<td class='cat-label' style='color:{color};'>{cat}</td>"
        for prefix in prefixes:
            row = df_subset[
                (df_subset["prefix"] == prefix) & (df_subset["category"] == cat)
            ]
            if row.empty:
                cells += "<td>-</td>"
            else:
                r       = row.iloc[0]
                n_total = int(r["n_total"])
                pct_inv = float(r["pct_inv"])
                cells += f"<td>{n_total:,} &nbsp;|&nbsp; {pct_inv}%</td>"
        rows_html += f"<tr>{cells}</tr>"

    return f"""
    <div class='table-section'>
        <h2>{title}</h2>
        <p class='subtitle'>{subtitle}</p>
        <table>
            <thead>
                <tr>{header_cells}</tr>
                <tr class='subheader'>{subheader_cells}</tr>
            </thead>
            <tbody>{rows_html}</tbody>
        </table>
    </div>
    """


def main():
    missing = [f for f in [filter_file, thinned_file] if not os.path.exists(f)]
    if missing:
        print(f"Missing input files: {missing}")
        return

    df_filter  = pd.read_csv(filter_file,  sep="\t")
    df_thinned = pd.read_csv(thinned_file, sep="\t")

    table1 = build_table(
        df_filter[df_filter["subset"] == "all"],
        "All SNPs",
        "All SNPs passing FDR threshold per selection category."
    )
    table2 = build_table(
        df_filter[df_filter["subset"] == "top10pct"],
        "Top 10% SNPs",
        "Top 10% most significant SNPs per category sorted by p-value."
    )
    table3 = build_table(
        df_thinned[df_thinned["subset"] == "all"],
        "All thinned SNPs",
        "All gene-thinned SNPs (one SNP per gene) per selection category."
    )
    table4 = build_table(
        df_thinned[df_thinned["subset"] == "top10pct"],
        "Top 10% thinned SNPs",
        "Top 10% most significant gene-thinned SNPs per category sorted by p-value."
    )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Enrichment of SNPs categories in inversions and collinear regions of the genome</title>
    <link rel="icon" type="image/png" href="../favicon_sbb.png">
    <style>
        body {{
            font-family: Arial, sans-serif;
            font-size: 13px; color: #222;
            background: #fff; padding: 40px;
            max-width: 1200px; margin: 0 auto;
        }}

        h1 {{
            font-size: 1.4rem;font-weight: 600; margin-bottom: 4px;
        }}

        .page-subtitle {{
            color: #666; font-size: 0.95rem;
            margin-bottom: 10px;line-height:1.4;
        }}

        .section-header {{
            font-size: 0.75rem; font-weight: 700;
            letter-spacing: 0.1em;text-transform: uppercase;
            color: #222; border-bottom: 2px solid #eee;
            padding-bottom: 6px;margin-bottom: 24px; margin-top: 28px;
        }}

        .table-section {{
            margin-bottom: 40px;
        }}

        .table-section h2 {{
            font-size: 0.95rem;
            font-weight: 600;
            margin-bottom: 4px;
        }}

        .subtitle {{
            color: #666;
            font-size: 0.92rem;
            margin-bottom: 12px;
        }}

        table {{
            border-collapse: collapse;
            width: 100%;
        }}

        thead tr:first-child th {{
            background: #f2f2f2;
            font-weight: 600;
            font-size: 12px;
            padding: 8px 14px;
            text-align: center;
            border: 1px solid #ddd;
        }}

        thead tr:first-child th:first-child {{
            text-align: left;
            width: 120px;
        }}

        thead tr.subheader th {{
            background: #fafafa;
            color: #888;
            font-size: 11px;
            font-weight: 400;
            padding: 4px 14px;
            text-align: center;
            border: 1px solid #ddd;
        }}

        tbody tr {{
            border-bottom: 1px solid #eee;
        }}

        tbody tr:hover {{
            background: #fafafa;
        }}

        tbody td {{
            padding: 8px 14px; text-align: center;
            border: 1px solid #eee; font-variant-numeric: tabular-nums;
        }}

        td.cat-label {{
            text-align: left; font-weight: 600;
            font-size: 12px;border: 1px solid #eee;
        }}

        .note {{
            color: #999; font-size: 11px; margin-top: 6px;
        }}

        hr {{
            border: none; border-top: 1px solid #eee; margin: 8px 0 40px;
        }}
    </style>
</head>
<body>
    <h1>Enrichment of SNPs categories in inversions and collinear regions of the genome</h1>
    <p class="page-subtitle">
    Proportion of SNPs in each category located within annotated genomic inversions 
    vs. collinear regions <br>  
    Tables show the number of SNPs (n) and the percentage of SNPs in each selection 
    category that fall within inversions (%inv). FDR &lt; 0.05 for selecting top SNPs.
    </p>

    <div class='section-header'>Full SNP set</div>
    {table1}
    <hr>
    {table2}

    <div class='section-header'> Gene-thinned SNP set</div>
    {table3}
    <hr>
    {table4}
</body>
</html>"""

    with open(output_file, "w") as f:
        f.write(html)
    print(f"Saved: {output_file}")


if __name__ == "__main__":
    main()