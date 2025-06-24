"""
Generates an html file with a table summarizing the Fastp report, one data sample per row
In: 
Out: ips_poolseq_report.html
"""
import subprocess
import sys
work_dir = "../../results/01_proc_reads/"
html = """<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Ips PoolSeq QC report</title>
    <link rel="icon" type="image/png" href="../favicon_sbb.png">
    <style>
        tr:nth-child(odd) {
            background-color: #ffffff; /* Light gray */
        }
        tr:nth-child(even) {
            background-color: #f2f2f2; /* White */
        }
        table { border-collapse: collapse; width: 100%; }
        th, td { border: 1px solid #ccc;  text-align: center; vertical-align: middle; padding: 0 8px; }
        th { background-color: #ddd; padding: 6px; border-color:#bbb}
        img { max-width: 250px; height: auto; }
        body {font-family: sans-serif; font-size:12px}
        .top_header { background-color: #d3d3d3; padding: 6px; border-color:#bbb}
        .reads{font-family:mono; font-size:10px; text-align:left}
        .red_samples {
        color: red;
        font-weight: bold;
    }
    </style>  
  <script>
  window.onload = function() {    
    const table = document.getElementById('Tab');  
      for (const row of table.tBodies[0].rows) {
      const cov = row.cells[21];
      const val = parseFloat(cov.textContent);
        const gc = row.cells[13];
        const val_gc = parseFloat(gc.textContent);

        if (!isNaN(val_gc)) {
          if (val_gc > 40) {
            gc.style.backgroundColor = '#FFD9E6';
        } else if(val_gc > 50){
        gc.style.backgroundColor = '#F4A6A6';
        }}

        const gcb = row.cells[21];
        const val_gcb = parseFloat(gcb.textContent);

        if (!isNaN(val_gcb)) {
          if (val_gcb > 40 && val_gcb <50) {
            gcb.style.backgroundColor = '#FFD9E6';
        } else if(val_gcb > 50){
        gcb.style.backgroundColor = '#F4A6A6';
        }}

      if (!isNaN(val)) {
        if (val < 100) {
          cov.style.backgroundColor = '#ff9999';
        } else if (val < 200) {
          cov.style.backgroundColor = '#ffcc99';
        } else if (val < 300) {
          cov.style.backgroundColor = '#ffffcc';
        } else {
          cov.style.backgroundColor = '#99ff99';  // val >= 300
        }
      }
    }
  };
</script>
</head>
<body>    
<br>
    <table id="Tab">
    <h3>Mapping output, before and after deduplication (plots are after)</h3>
    <thead>
        <tr>
        <th colspan=6; class="top_header">Dataset</th>
        <th colspan=13; class="top_header">Reads before filtering</th>
        <th colspan=12; class="top_header">Reads after filtering</th>
        </tr>
        <tr>
        <th>Prefix</th>
        <th>Year</th>
        <th>Time</th>
        <th>Country</th>
        <th>Region</th>
        <th>Rep</th>
        <th>Qubit</th>        
               
        <th>Reads before (M)</th>
        <th>Reads discarded (%)</th>
        <th>Reads discarded (M)</th>
        <th>Low qual (M)</th>
        <th>Many Ns (M)</th>
        <th>Short reads (M)</th>
        <th>GC before(%)</th>
        <th>Q20 before</th>     
        <th>Insert bp</th>
        <th>Duplics (%)</th>   
        <th>Quality plot before QC</th>
        <th>Base ratio plot before QC</th>

        <th>Reads after (M)</th>
        <th>Gbp total</th>
        <th>Coverage</th>
        <th>Read length</th>
        <th>GC (%)</th>
        <th>Q20 after</th>
        <th>Quality plot after QC</th>
        <th>Base ratio plot after QC</th>

        <th>Files</th>       
        <th>Overrepresented reads (in 1e6 sample) </th>   
        </tr> 
        </thead>
  <tbody>"""

with open(f"{work_dir}/all_poolseq_report.tsv", "r") as fh:
    next(fh)
    for line in fh:
        idn, year, time, country, region, rep, after, gbp, cov, length, gc, insert , \
        dups, q20a, before, lost, lqual, ns, short, q20b, qubit, pct_lost, gc_before, file = line.strip().split("\t")
    
        top5_reads_out = subprocess.run(f"grep -A12 {idn} {work_dir}top5_reads | egrep -v '^0|_' ", shell=True, capture_output=True, text=True)
        top5_reads = top5_reads_out.stdout.replace("1.fq", "").replace("2.fq", "")
    
        plotQ = "./before/" + file.replace(".fq.gz",".qual.png").split(',')[0]
        plotB = "./before/" + file.replace(".fq.gz",".base.png").split(',')[0]
        plotQa = "./after/" + file.replace(".fq.gz",".qual.png").split(',')[0]
        plotBa = "./after/" + file.replace(".fq.gz",".base.png").split(',')[0]

        html += f"""
            <tr>
                <td>{idn}</td>
                <td>{year}</td>
                <td>{time}</td>
                <td>{country}</td>
                <td>{region}</td>
                <td>{rep}</td>
                <td>{qubit}</td>   
                                                
                <td>{before}</td>
                <td>{pct_lost}</td>
                <td>{lost}</td>
                <td>{lqual}</td>
                <td>{ns}</td>
                <td>{short}</td>
                <td>{gc_before}</td>
                <td>{q20b}</td>
                <td>{insert}</td>
                <td>{dups}</td>
                <td><img src="{plotQ}" /></td>
                <td><img src="{plotB}" /></td>

                        <td>{after}</td>
                <td>{gbp}</td>
                <td>{cov}</td>
                <td>{length}</td>
                <td>{gc}</td>
                <td>{q20a}</td>
                <td><img src="{plotQa}" /></td>
                <td><img src="{plotBa}" /></td>

                <td>{file.replace(',','<br>')}</td>
                <td class="reads"><pre>{top5_reads}</pre></td>            
            </tr>
        """
html += """
   </tbody> </table>
   </body>
   <script>
    const keywords = ["SFIN_Ea_2015", "SFIN_Eb_2015", "WFIN_Ea_2015", "WFIN_Eb_2015", 
    "SFIN_La_2016", "SFIN_Lb_2016r", "SFIN_Ea_2017", "SFIN_Eb_2017", "SFIN_La_2017", "SFIN_Lb_2017"];

    document.querySelectorAll("td").forEach(cell => {
        keywords.forEach(keyword => {
            if (cell.textContent.includes(keyword)) {
                cell.classList.add("red_samples");
            }
        });
    });
</script>
</html>
"""
# Write to file
with open(f"{work_dir}ips_poolseq_report.html", "w", encoding="utf-8") as f:
    f.write(html)
