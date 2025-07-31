#!/usr/bin/env python3

#----------------------------------------------------------------------
# Builds HTML table of null variance and divergence calculations.
# In:
#   - null_variance_summary.tsv: summary statistics for dz2 and depth variance
# Outputs:
#   - null_variance_summary.html: HTML table with null variance
#----------------------------------------------------------------------

import csv

work_dir = "../../results/06_SNPs_stats"

# Input file
null_var_in = f"{work_dir}/null_variance_summary.tsv"
# Output file
null_var_table_out = f"{work_dir}/null_variance_summary.html"

def read_tsv(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        return list(reader)

def tsv_to_html_table(rows):
    html = ['<table>']

    # Header row
    html.append('<tr>')
    for cell in rows[0]:
        html.append(f'<th>{cell}</th>')
    html.append('</tr>')

    # Data rows
    for row in rows[1:]:
        html.append('<tr>')
        for cell in row:
            html.append(f'<td>{cell}</td>')
        html.append('</tr>')

    html.append('</table>')
    return '\n'.join(html)

def write_html(html_content, output_path):
    html = f"""<html><head><meta charset="UTF-8">
  <style>
    table {{ border-collapse: collapse; width: 100%; font-family: sans-serif; }}
    th, td {{ border: 1px solid #999; padding: 4px; text-align: right; }}
    th {{ background-color: #f2f2f2; }}
  </style>
</head>
<body>
{html_content}
</body>
</html>"""
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)

def main():
    rows = read_tsv(null_var_in)
    html_table = tsv_to_html_table(rows)
    write_html(html_table, null_var_table_out)

if __name__ == "__main__":
    main()
