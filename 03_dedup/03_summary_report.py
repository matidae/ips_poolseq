import os, sys
import csv
import html
sys.path.append("./utils")
from utils import prefix_order


"""
Generate HTML report for mapping and deduplication metrics.

Input:
    - ../results/03_dedup/summary_table.tsv
    - ../results/03_dedup/depth_metrics/*.depth.500K_bins.png
    - ../results/03_dedup/insertion_size_metrics/*.insertion_metrics_histogram.png

Output:
    - ../results/03_dedup/map_and_dedup.html
"""

# Working directory 
out_dir = "../results/03_dedup"

html_top = """
<html>
  <head>
  <title>Ips PoolSeq QC report</title>
  <link rel="icon" type="image/png" href="../favicon_sbb.png">
    <style>
      tr:nth-child(odd) { background-color: #fff; }
      tr:nth-child(even) { background-color: #f2f2f2; }

      table { border-collapse: collapse; width: auto; }

      th, td {
        border: 1px solid #ccc;
        text-align: center;
        vertical-align: middle;
        font-size: 12px;
        padding: 4px 5px;
      }
      th { background-color: #ddd; border-color: #bbb; }
      img { max-width: 150px; height: auto; }
      body { font-family: sans-serif; font-size: 12px; }
      .top_header { background-color: #d3d3d3; padding: 6px; border-color: #bbb; }
      .title { font-family: sans-serif; font-size: 18px; font-weight: bold; }
      .legend {border:1px solid #ccc; margin-right:6px; display:inline-block; width:12px; height:12px;align-items: center;}
      img:hover { transform: scale(3);
                    z-index: 10; 
                    position: relative;}
      .square_image { max-width: 75px }
    </style>   
  </head>
  <body>
<br><span class="title">Mapping output, before and after deduplication.<br><br></span>  
  <div style="display: flex; gap: 16px; align-items: center; font-family: sans-serif;">
    <div><span><strong>Mean depth:</strong></span></div>
    <div><span class="legend"; style="background:#ff9999;"></span> 0 - 100</div>
    <div><span class="legend"; style="background:#ffcc99;"></span> 100 - 200  </div>
    <div><span class="legend"; style="background:#ffffcc;"></span> 200 - 300 </div>
    <div><span class="legend"; style="background:#99ff99;"></span> > 300 </div>
  </div> 
  <br>
"""
js_coloring = """
    <script>
      document.querySelectorAll('td.cov_cell').forEach(function(cov) {
        let val = parseFloat(cov.textContent);
        if (!isNaN(val)) {
          if (val < 100) {
            cov.style.backgroundColor = '#ff9999';
          } else if (val < 200) {
            cov.style.backgroundColor = '#ffcc99';
          } else if (val < 300) {
            cov.style.backgroundColor = '#ffffcc';
          } else {
            cov.style.backgroundColor = '#99ff99';
          }
        }
      });
    </script>
"""

#Columns in HTML table
#Idn Year	Time	Country	Region	Reps	Before After	Mapped Dedup Cov	Cov_mapped	Cov_dedup

#Write table header
def write_table_header(out, header):
    out.write("<table border='1'>\n")
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
        if col_idx in [8]:  # columns 11 and 12 with cov_cell class for color
            out.write(f'<td class="cov_cell">{html.escape(col)}</td>\n')
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
    table_path = os.path.join(out_dir, "summary_table.tsv")
    output_path = os.path.join(out_dir, "map_and_dedup.html")
    
    # Generate report
    generate_html_report(table_path, output_path)

if __name__ == "__main__":
    main()