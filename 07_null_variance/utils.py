import csv

depths_dir = "../../results/05_SNPs_depths"
stats_dir = "../../results/06_SNPs_stats"

def parse_counts(ad_field):    
    ref, alt = map(int, ad_field.split(','))
    return ref, alt

def load_paired_samples():
    samples_reps = {}
    with open(f"{stats_dir}/samples_paired.tsv") as fh:
        for line in fh:
            sample, idx1, idx2 = line.strip().split('\t')
            samples_reps[sample] = [int(idx1), int(idx2)]
    return samples_reps

def load_depth_threshold():
    depth_threshold = f"{depths_dir}/genic_depth_stats.tsv"
    with open(depth_threshold ,"r") as f:
        for line in f:
            pass
        min_depth, max_depth = map(float, line.strip().split('\t')[-2:])
        return min_depth, max_depth

def load_null_variance_recalc():
    null_var = {}
    null_var_file = f"{stats_dir}/null_variance_summary.recalc.tsv"    
    with open(null_var_file, newline='') as null_var_fh:
        next(null_var_fh)
        null_var_read = csv.reader(null_var_fh, delimiter='\t')
        for row in null_var_read:
            null_var[row[0]] = float(row[1])
    return null_var

def tsv_to_html_table(tsv_in, title):
    # Open TSV to build HTML table
    with open(tsv_in, 'r', encoding='utf-8') as tsv_fh:
        rows = list(csv.reader(tsv_fh, delimiter='\t'))

    table = ['<table>']
    # Header row
    table.append('<tr>')
    for cell in rows[0]:
        table.append(f'<th>{cell}</th>')
    table.append('</tr>')
    # Data rows
    for row in rows[1:]:
        table.append('<tr>')
        for cell in row:
            table.append(f'<td>{cell}</td>')
        table.append('</tr>')
    table.append('</table>')    
    table_content =  '\n'.join(table)
    
    html_header = f"""<html><head><meta charset="UTF-8">
  <style>
    table {{ border-collapse: collapse; width: 100%; font-family: sans-serif; }}
    th, td {{ border: 1px solid #999; padding: 4px; text-align: right; }}
    th {{ background-color: #f2f2f2; }}
    h2 {{ font-family: sans-serif; }}
  </style>
</head>
<body>
<h2>{title}</h2>
{table_content}
</body>
</html>"""
    html_out =  tsv_in.replace("tsv", "html")
    with open(html_out, 'w', encoding='utf-8') as f:
        f.write(html_header)