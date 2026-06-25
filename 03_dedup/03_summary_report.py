#!/usr/bin/env python3

#----------------------------------------------------------------------
# Generate HTML report for mapping and deduplication metrics.
#
# Input:
#    - $DEDUP_RESULTS/summary_table_all.tsv
#    - $DEDUP_RESULTS/depth_metrics/$prefix.depth.500K_bins.png
#    - $DEDUP_RESULTS/insertion_size_metrics/$prefix.dedup.insertion_metrics_histogram.png
#
# Output:
#    - $DEDUP_RESULTS/map_and_dedup.html
#----------------------------------------------------------------------

import os, sys
import csv
import html
sys.path.append("./utils")
from utils import prefix_order, load_config, log


# Working directory
cfg = load_config()
out_dir = cfg["DEDUP_RESULTS"]

html_top = """
<html>
  <head>
  <title>Ips PoolSeq: Mapping & Deduplication QC</title>
  <link rel="icon" type="image/png" href="../favicon_sbb.png">
    <style>
      body { font-family: system-ui, -apple-system, "Segoe UI", sans-serif; font-size:13px; color:#2b2b2b; margin:24px 32px }
      h2 { font-weight:600; letter-spacing:0.2px; margin-bottom:4px }
      table { border-collapse:separate; border-spacing:0; width:95%; border:1px solid #e3e3e3; border-radius:8px; overflow:hidden }
      th, td { text-align:center; vertical-align:middle; font-size:13px; padding:1px 0px; border-bottom:1px solid #ededed }      
      thead th { position:sticky; top:0; background-color:#f4f4f5; font-weight:600; border-bottom:2px solid #d8d8d8; z-index:2; padding:10px 10px; }
      .top_header { background-color:#ececee; font-weight:600; color:#555; letter-spacing:0.3px }
      tbody tr:hover { background-color:#f0f6fb; }
      td:first-child { font-weight:600; text-align:left; padding-left:14px }
      td:not(:first-child) { font-variant-numeric: tabular-nums }
      img { max-width:120px; height:auto; border-radius:4px; transition:transform 0.15s ease }
      img:hover { transform:scale(5); z-index:10; position:relative; transform-origin:right center }
      .square_image { max-width:50px }
      .legend-box { padding:2px 10px; border:1px solid #ddd; border-radius:3px; display:inline-block }
      thead th:hover { background-color: #e8e8ea; }
      
    </style>
  </head>
  <body>
<br>
<div style="align-items:center; margin-bottom:10px;">
  <h2 style="margin:0 0 16px 0;">Mapping output, before and after deduplication</h2>
<div style="font-size:12px;">
  <span style="font-weight:bold;">Mean depth:</span>&nbsp;&nbsp;
  <span class="legend-box" style="background:#f5c8c8;"></span>&nbsp;&lt;100&nbsp;&nbsp;
  <span class="legend-box" style="background:#f5e0c8;"></span>&nbsp;100-200&nbsp;&nbsp;
  <span class="legend-box" style="background:#f5f5c8;"></span>&nbsp;200-300&nbsp;&nbsp;
  <span class="legend-box" style="background:#c8f0c8;"></span>&nbsp;&gt;300
  <br><br>
  <span style="font-weight:bold;">Read length:</span>&nbsp;&nbsp;
  <span class="legend-box" style="background:#f5c8c8;"></span>&nbsp;&lt;120&nbsp;&nbsp;
  <span class="legend-box" style="background:#f5e0c8;"></span>&nbsp;120-135&nbsp;&nbsp;
  <span class="legend-box" style="background:#f5f5c8;"></span>&nbsp;135-145&nbsp;&nbsp;
  <span class="legend-box" style="background:#c8f0c8;"></span>&nbsp;&gt;145
</div>
</div>
<br>
"""
js_coloring = """
  <script>
  //Color depth cells based on thresholds
document.querySelectorAll('td.cov_cell').forEach(function(cov) {
    let val = parseFloat(cov.textContent);
    if (!isNaN(val)) {
        if (val < 100)       cov.style.backgroundColor = '#f5c8c8';
        else if (val < 200)  cov.style.backgroundColor = '#f5e0c8';
        else if (val < 300)  cov.style.backgroundColor = '#f5f5c8';
        else                 cov.style.backgroundColor = '#c8f0c8';
    }
});
    // Color read length column based on thresholds
document.querySelectorAll('td.length_cell').forEach(function(cell) {
    let val = parseFloat(cell.textContent);
    if (!isNaN(val)) {
        if (val < 120)       cell.style.backgroundColor = '#f5c8c8';
        else if (val < 135)  cell.style.backgroundColor = '#f5e0c8';
        else if (val < 145)  cell.style.backgroundColor = '#f5f5c8';
        else                 cell.style.backgroundColor = '#c8f0c8';
    }
});
    // Group rows by prefix without replicate (e.g. LAAU_E_2023)
    const rows = document.querySelectorAll('tbody tr');
    let groupColor = false;
    let lastGroup = null;

    rows.forEach(row => {
      const id = row.cells[0].textContent.trim();
      // Strip trailing a/b to get group key
      const group = id.replace(/([EL])[ab]_/, '$1_');
      if (group !== lastGroup) {
        groupColor = !groupColor;
        lastGroup = group;
      }
      row.style.backgroundColor = groupColor ? '#ffffff' : '#fafafa';
    });

    // Add sorting functionality to table headers
    document.querySelectorAll('thead th').forEach(function(th, idx) {
    th.style.cursor = 'pointer';
    th.addEventListener('click', function() {
        const tbody = document.querySelector('tbody');
        const rows = Array.from(tbody.querySelectorAll('tr'));
        const asc = th.dataset.sort !== 'asc';
        th.dataset.sort = asc ? 'asc' : 'desc';
        rows.sort(function(a, b) {
            const va = parseFloat(a.cells[idx].textContent);
            const vb = parseFloat(b.cells[idx].textContent);
            if (isNaN(va) || isNaN(vb)) return 0;
            return asc ? va - vb : vb - va;
        });
        rows.forEach(function(row) { tbody.appendChild(row); });
    });
});
  </script>
"""

#Columns in HTML table
#Idn Year	Time	Country	Region	Reps	Before After	Mapped Dedup Cov	Cov_mapped	Cov_dedup

#Write table header
def write_table_header(out, header):
    out.write("<table>\n")
    out.write("<thead><tr>\n")
    for col in header:
        out.write(f"<th>{html.escape(col)}</th>\n")
    out.write("<th>Depth plot</th>\n")
    out.write("<th>Insertion size</th>\n")
    out.write("</tr></thead>\n")

def write_table_row(out, row):
    #Write table body
    #  for row in data_rows:
    out.write("<tr>\n")
    for col_idx, col in enumerate(row):
        if col_idx in [8]:  # column depth_cell class for coloring
            out.write(f'<td class="cov_cell">{html.escape(col)}</td>\n')
        elif col_idx in [4]: # column read length for coloring
          out.write(f'<td class="length_cell">{html.escape(col)}</td>\n')    
        else:            
            out.write(f"<td>{html.escape(col)}</td>\n")
    # Add image cell
    sample_id = row[0]  #get prefix for plot files

    depth_img = f"./depth_metrics/{sample_id}.depth.500K_bins.png"
    out.write(f' <td><img src="{depth_img}" alt="{depth_img}"></td>\n')

    insertion_size_img = f"./insertion_size_metrics/{sample_id}.dedup.sort.insertion_metrics_histogram.png"
    out.write(f'<td><img class="square_image" src="{insertion_size_img}" style="background-color: white;" ></td>\n')
    out.write("</tr>\n")

def generate_html_report(table_path, output_path):  
  # Read input table
  with open(table_path, 'r') as tsvfile:
      reader = csv.reader(tsvfile, delimiter='\t')
      rows = list(reader)
  
  prefix_sort = prefix_order()
  header = rows[0]
  data_rows = rows[1:]
  data_rows.sort(key=lambda row: prefix_sort[row[0]]) 
  
  # Write HTML report
  with open(output_path, 'w') as out:
      # Write HTML top
      out.write(html_top)
      
      # Write table header
      write_table_header(out, header)
      
      # Write table body
      out.write("<tbody>\n")
      for row in data_rows:
          write_table_row(out, row)#, MAPPED_READS_COL, DEDUP_READS_COL)
      out.write("</tbody>\n</table>\n")
      
      # Add JavaScript for coloring
      out.write(js_coloring)
      out.write("</body>\n</html>\n")


def main():
    # Define paths
    table_path = os.path.join(out_dir, "summary_table_all.tsv")
    output_path = os.path.join(out_dir, "map_and_dedup.html")
    log(f"=== Generating HTML report for mapping and deduplication metrics ===")
    # Generate report
    log(f"done: {output_path}")
    generate_html_report(table_path, output_path)
    log("=== HTML report generation complete ===")

if __name__ == "__main__":
    main()