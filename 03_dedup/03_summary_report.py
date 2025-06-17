import csv
import html

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
        font-size: 14px;
        padding: 5px 6px;
      }
      th { background-color: #ddd; border-color: #bbb; }
      img { max-width: 300px; height: auto; }
      body { font-family: sans-serif; font-size: 14px; }
      .top_header { background-color: #d3d3d3; padding: 6px; border-color: #bbb; }
      .title { font-family: sans-serif; font-size: 18px; font-weight: bold; }
      .legend {border:1px solid #ccc; margin-right:6px; display:inline-block; width:12px; height:12px;align-items: center;}
    </style>   
  </head>
  <body>
<br><br><span class="title">Mapping output, before and after deduplication.<br><br></span>


  
  <div style="display: flex; gap: 16px; align-items: center; font-family: sans-serif;">
    <div><span><strong>Mean depth:</strong></span></div>
    <div><span class="legend"; style="background:#ff9999;"></span> 0 - 100</div>
    <div><span class="legend"; style="background:#ffcc99;"></span> 100 - 200  </div>
    <div><span class="legend"; style="background:#ffffcc;"></span> 200 - 300 </div>
    <div><span class="legend"; style="background:#99ff99;"></span> > 300 </div>
  </div> 
  <br>
"""
table = "../../results/03_dedup/summary_table"
#Columns in HTML table
#Idn Year	Time	Country	Region	Reps	Before After	Mapped Dedup Cov	Cov_mapped	Cov_dedup
out_html_report = "../../results/03_dedup/map_and_dedup.html"
#Open file for reading and for writing
with open(table) as tsvfile, open(out_html_report, "w") as out:
    reader = csv.reader(tsvfile, delimiter='\t')
    rows = list(reader)
    header = rows[0]
    data_rows = rows[1:]    
    
    #Write HTML top and css
    out.write(html_top)

    #Write table header
    out.write("<table border='1'>\n")
    out.write("  <thead><tr>\n")
    for col in header:
        out.write(f"<th>{html.escape(col)}</th>\n")
    out.write("<th>Depth</th>\n")
    out.write("<th>Insertion size</th>\n")
    out.write("</tr></thead>\n")

    #Write table body
    out.write("<tbody>\n")
    for row in data_rows:
        out.write("<tr>\n")
        for col_idx, col in enumerate(row):
            if col_idx in [10, 11]:  # columns 11 and 12 with cov_cell class for color
                out.write(f'<td class="cov_cell">{html.escape(col)}</td>\n')
            else:
                out.write(f"<td>{html.escape(col)}</td>\n")
        # Add image cell
        sample_id = row[0]  #get prefix for plot files

        depth_img = f"./depth/{sample_id}.depth.500K_bins.png"
        out.write(f' <td><img src="{depth_img}" alt="{depth_img}"></td>\n')

        insertion_size_img = f"./insertion_size_metrics/{sample_id}.dedup.sort.insertion_metrics_histogram.png"
        out.write(f'<td><img src="{insertion_size_img}" width="160px"></td>\n')
        out.write("</tr>\n")
    out.write("</tbody>\n</table>\n")

#Add JS script for coloring:
    out.write("""
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
""")
