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
    table.append('<thead><tr>')
    for cell in rows[0]:
        table.append(f'<th>{cell}</th>')
    table.append('</tr></thead><tbody>')
    # Data rows
    for row in rows[1:]:
        table.append('<tr>')
        for cell in row:
            table.append(f'<td>{cell}</td>')
        table.append('</tr>')
    table.append('</tbody></table>')    
    table_content =  '\n'.join(table)
    
    html_header = f"""<html><head><meta charset="UTF-8">
    <style>
    table {{ border-collapse: collapse; width: 100%; font-family: sans-serif; }}
    th, td {{ border: 1px solid #999; padding: 4px; text-align: right; }}
    th {{ background-color: #f2f2f2; cursor: pointer; }}
    th:hover {{ background-color: #e0e0e0; }}
    h2 {{ font-family: sans-serif; }}
  </style>
    <script>
    // Simple sortable table script
    document.addEventListener('DOMContentLoaded', () => {{
      const getCellValue = (tr, idx) => tr.children[idx].innerText || tr.children[idx].textContent;

      const comparer = (idx, asc) => (a, b) => {{
        const v1 = getCellValue(asc ? a : b, idx);
        const v2 = getCellValue(asc ? b : a, idx);
        const n1 = parseFloat(v1.replace(/,/g, ''));
        const n2 = parseFloat(v2.replace(/,/g, ''));
        const bothNumeric = !isNaN(n1) && !isNaN(n2);
        return bothNumeric ? n1 - n2 : v1.localeCompare(v2);
      }};

      document.querySelectorAll('th').forEach(th => th.addEventListener('click', (() => {{
        const table = th.closest('table');
        const tbody = table.querySelector('tbody');
        Array.from(tbody.querySelectorAll('tr'))
          .sort(comparer(Array.from(th.parentNode.children).indexOf(th), this.asc = !this.asc))
          .forEach(tr => tbody.appendChild(tr));
      }})));
    }});
  </script>
</head>
<body>
<h2>{title}</h2>
{table_content}
</body>
</html>"""
    html_out =  tsv_in.replace("tsv", "html")
    with open(html_out, 'w', encoding='utf-8') as f:
        f.write(html_header)