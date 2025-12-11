"""
Reads fastp JSON output files and generate summary tables and plots.

Input:
    - fastp JSON report files (e.g., *.json)

Output:
    - all_poolseq_report.tsv : summary table with quality metrics for all samples
    - $prefix_quality.png : quality score plot for each sample
    - $prefix_base_content.png : Base composition plot for each sample
    - ips_poolseq_report.html : html report with quality metrics for all samples
"""

import os
import glob
import json
import matplotlib.pyplot as plt

in_dir = "../data/01_proc_reads"
out_dir = "../results/01_proc_reads"

#Ips typographus genome size
ips_genome_size = 224977219

# Input
qubit_in = "../data/reference/sample_qubit"
json_list = f"{in_dir}/*/*.json"

# Output
report_out = f"{out_dir}/all_poolseq_report.tsv"


def plot_quality(quality_curves, filename, filter):   
   pos=list(range(1,151))
   plt.rcParams['font.size'] = 6
   plt.figure(figsize=(4, 2))
   plt.ylim(0,40)
   plt.ylabel("quality")
   plt.xlabel("position")   
   plt.title(filename)
   plt.plot(pos, quality_curves, linestyle="-", color="darkblue", lw=1)
   plt.style.use("ggplot")
   plt.tight_layout()
   plot_filename = f"{out_dir}/{filter}/" + filename.replace('.fq.gz','.qual.png')
   plt.savefig(plot_filename, dpi=300)
   plt.close()

def plot_base_content(content_curves, filename, filter):   
   pos=list(range(1,151))   
   plt.rcParams['font.size'] = 6
   plt.figure(figsize=(4, 2))
   plt.ylim(0,1)   
   plt.ylabel("base ratio")
   plt.title(filename)
   plt.xlabel("position")   
   plt.plot(pos, content_curves["A"], linestyle="-", color="darkblue", lw=1)
   plt.plot(pos, content_curves["T"], linestyle="-", color="darkred", lw=1)
   plt.plot(pos, content_curves["G"], linestyle="-", color="darkgreen", lw=1)
   plt.plot(pos, content_curves["C"], linestyle="-", color="orange", lw=1)
   plt.legend(["A", "T", "G", "C"], loc="upper right")
   plt.style.use("ggplot")   
   plt.tight_layout()
   plot_filename = f"{out_dir}/{filter}/" + filename.replace('.fq.gz','.base.png')
   plt.savefig(plot_filename, dpi=300)
   plt.close()

def parse_json(json_files):    
    data = {}
    sample_qubit = {i.rstrip().split()[1] : i.split()[0] for i in open(qubit_in)}

    #Parse JSON files from Fastp to retrieve QC info.
    for file in json_files:
        with open(file, 'r') as fh:
            fastp_out = json.load(fh)
            prefix = file.split("/")[-2]
            prefix_short = prefix[:12]
            idn = "_".join(prefix.split("_")[:3])
            
            reads_after = fastp_out["summary"]["after_filtering"]["total_reads"]
            total_bp = fastp_out["summary"]["after_filtering"]["total_bases"]
            coverage = total_bp/ips_genome_size
            length = (fastp_out["summary"]["after_filtering"]["read1_mean_length"] + \
                       fastp_out["summary"]["after_filtering"]["read2_mean_length"])*0.5
            gc_content = fastp_out["summary"]["after_filtering"]["gc_content"]
            
            dups = fastp_out["duplication"]["rate"]

            reads_before = fastp_out["summary"]["before_filtering"]["total_reads"]
            reads_lost = reads_before - reads_after
            low_qual = fastp_out["filtering_result"]["low_quality_reads"]
            too_short = fastp_out["filtering_result"]["too_short_reads"]
            
            filename = file.split('/')[-1].replace("qc_report.json", "fq.gz")
            quality_curves_after = fastp_out["read1_after_filtering"]["quality_curves"]["mean"]
            content_curves_after = fastp_out["read1_after_filtering"]["content_curves"]
            # Plots quality and base content per read after Fastp QC
            plot_quality(quality_curves_after, filename, "after")
            plot_base_content(content_curves_after, filename, "after")

            quality_curves_before = fastp_out["read1_before_filtering"]["quality_curves"]["mean"]
            content_curves_before = fastp_out["read1_before_filtering"]["content_curves"]
            # Plots quality and base content per read before Fastp QC
            plot_quality(quality_curves_before, filename, "before")
            plot_base_content(content_curves_before, filename, "before")

            qubit = sample_qubit[prefix_short]
            gc_content_before = fastp_out["summary"]["before_filtering"]["gc_content"]
            pct_reads_lost = round(100 - (reads_after/reads_before)*100, 1)
            
            row = [idn, reads_after, total_bp, coverage, length, gc_content, dups, \
                reads_before, reads_lost, low_qual, too_short, qubit, pct_reads_lost, gc_content_before, filename]
                
            if idn not in data.keys():
                data[idn] = row
            else:
                #sum:  "reads_after", "total_bp", "reads_before", "reads_lost", "low_qual",  "too_short"
                index_to_sum = [1, 2, 7, 8, 9, 10]

                for k in index_to_sum:
                    data[idn][k] += row[k]
                #Calculate coverage
                coverage = data[idn][2]/ips_genome_size
                data[idn][3] = coverage
                #Append filename
                data[idn][-1] += "," + row[-1]
                #Reads lost
                pct_reads_lost = round(100 - (data[idn][1]/data[idn][7])*100, 1)
                data[idn][12] = pct_reads_lost
    return data

def write_tsv_report(data):
    col_names = ["Idn", "After", "Gbp", "Cov", "Length", "GC", "Insert", "Dups", \
                "Before",  "Lost", "LQual", "Short", "Qubit", "%Reads_lost", "GCbefore", "File"]
    header = "\t".join(col_names)
    data_sort = sorted(data.items(), key=lambda x: x[1][0])
    data_sort_new=[]

    report_fh = open(report_out, 'w')
    report_fh.write(header + '\n')
    for k, v in data_sort:
        elems = [
            v[0],                              # idn
            round(v[1]/1e6, 1),                # reads_after (M)
            round(v[2]/1e9, 1),                # total_bp (Gb)
            round(v[3], 1),                    # coverage
            int(v[4]),                         # length
            round(v[5] * 100, 1),              # gc_content %
            round(v[6] * 100, 1),              # dups %
            round(v[7]/1e6, 1),                # reads_before (M)
            round(v[8]/1e6, 1),                # reads_lost (M)
            round(v[9]/1e6, 1),                # low_qual (M)
            round(v[10]/1e6, 1),               # too_short (M)
            v[11],                             # qubit
            v[12],                             # pct_reads_lost
            round(v[13] * 100, 1),             # gc_content_before %
            v[14]
        ]            
        elems_str = list(map(str, elems))
        data_sort_new.append(elems)
        report_fh.write("\t".join(elems_str) + "\n")
    report_fh.close()

def write_html_report():
    html = """<!DOCTYPE html>
    <html>
    <head>
      <meta charset="utf-8">
      <title>Fastp QC summary: before vs after filtering </title>
      <link rel="icon" type="image/png" href="../favicon_sbb.png">
      <style>
          tr:nth-child(odd) { background-color: #ffffff; }
          tr:nth-child(even) { background-color: #f2f2f2; }
          table { border-collapse: collapse; width: 100%; }
          th, td { border: 1px solid #ccc;  text-align: center; vertical-align: middle; padding: 0 8px; }
          th { background-color: #ddd; padding: 6px; border-color:#bbb}
          img { max-width: 200px; height: auto; }
          body {font-family: sans-serif; font-size:12px}
          .top_header { background-color: #d3d3d3; padding: 6px; border-color:#bbb}
          .reads{font-family:mono; font-size:10px; text-align:left}
          .legend-box { padding: 2px 8px; border: 1px solid #ccc; display: inline-block; }

        }
       </style>  
    <script>
    window.onload = function () {
        const table = document.getElementById('Tab');
        for (const row of table.tBodies[0].rows) {        
            const col_gc_before = 7;            
            const col_cov = 13;
            const col_gc_after = 15;

            // Color by  Depth
            const cov_cell = row.cells[col_cov];
            const cov_val = parseFloat(cov_cell.textContent);

            if (!isNaN(cov_val)) {
                if (cov_val < 100) {
                    cov_cell.style.backgroundColor = '#ff9999';
                } else if (cov_val < 200) {
                    cov_cell.style.backgroundColor = '#ffcc99';
                } else if (cov_val < 300) {
                    cov_cell.style.backgroundColor = '#ffffcc';
                } else {
                    cov_cell.style.backgroundColor = '#99ff99';
                }
            }
            // Color by GC before
            const gcb_cell = row.cells[col_gc_before];
            const gcb_val = parseFloat(gcb_cell.textContent);
            if (!isNaN(gcb_val)) {
                if (gcb_val > 50) {
                    gcb_cell.style.backgroundColor = '#F4A6A6';
                } else if (gcb_val > 40) {
                    gcb_cell.style.backgroundColor = '#FFD9E6';
                }
            }
            // Color by GC after
            const gca_cell = row.cells[col_gc_after];
            const gca_val = parseFloat(gca_cell.textContent);
            if (!isNaN(gca_val)) {
                if (gca_val > 50) {
                    gca_cell.style.backgroundColor = '#F4A6A6';
                } else if (gca_val > 40) {
                    gca_cell.style.backgroundColor = '#FFD9E6';
                }
            }
        }
    };
    </script>
    </head>
    <body>    
      <br>
      <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:10px;">
        <h2 style="margin:0;">Fastp QC summary: before vs after filtering</h2>
        <div style="font-size:11px; white-space:nowrap;">
          <span style="font-weight:bold;">Depth:</span>&nbsp;&nbsp;
          <span class="legend-box" style="background:#ff9999;"></span>&nbsp;&lt;100&nbsp;&nbsp;
          <span class="legend-box" style="background:#ffcc99"></span>&nbsp;100-199&nbsp;&nbsp;
          <span class="legend-box" style="background:#ffffcc"></span>&nbsp;200-299&nbsp;&nbsp;
          <span class="legend-box" style="background:#99ff99"></span>&nbsp;≥300&nbsp;&nbsp;
          <span style="font-weight:bold;">&nbsp;&nbsp;&nbsp;&nbsp;GC:</span>&nbsp;&nbsp;
          <span class="legend-box" style="background:#FFD9E6"></span>&nbsp;40-50%&nbsp;&nbsp;
          <span class="legend-box" style="background:#F4A6A6"></span>&nbsp;&gt;50%
        </div>
    </div>
    <table id="Tab">
      <thead>
          <tr>
            <th colspan=2; class="top_header"></th>
            <th colspan=9; class="top_header">Reads before filtering</th>
            <th colspan=9; class="top_header">Reads after filtering</th>
          </tr>
          <tr>
            <th>Prefix</th> 
            <th>Qubit</th>
            <th>Reads before (M)</th>
            <th>Discarded (M)</th>
            <th>Discarded (%)</th>
            <th>Low qual (M)</th> 
            <th>Short (M)</th>
            <th>GC (%)</th>
            <th>Duplics (%)</th>
            <th>Quality plot before QC</th>
            <th>Base ratio plot before QC</th>
            <th>Reads after (M)</th>
            <th>Gbp total</th>
            <th>Depth</th>
            <th>Length</th>
            <th>GC (%)</th>
            <th>Quality plot after QC</th>
            <th>Base ratio plot after QC</th>
          </tr> 
        </thead>
    <tbody>"""
    with open(f"{out_dir}/all_poolseq_report.tsv", "r") as fh:
        next(fh)
        for line in fh:        
            idn, after, gbp, cov, length, gc, dups, before, lost, \
            lqual, short, qubit, pct_lost, gc_before, file = line.strip().split("\t")
        
            plotQ = "./before/" + file.replace(".fq.gz",".qual.png").split(',')[0]
            plotB = "./before/" + file.replace(".fq.gz",".base.png").split(',')[0]
            plotQa = "./after/" + file.replace(".fq.gz",".qual.png").split(',')[0]
            plotBa = "./after/" + file.replace(".fq.gz",".base.png").split(',')[0]

            html += f"""
                <tr>
                  <td>{idn}</td>
                  <td>{qubit}</td>
                  <td>{before}</td>
                  <td>{lost}</td>
                  <td>{pct_lost}</td>
                  <td>{lqual}</td>
                  <td>{short}</td>
                  <td>{gc_before}</td>
                  <td>{dups}</td>
                  <td><img src="{plotQ}" /></td>
                  <td><img src="{plotB}" /></td>
                  <td>{after}</td>
                  <td>{gbp}</td>
                  <td>{cov}</td>
                  <td>{length}</td>
                  <td>{gc}</td>
                  <td><img src="{plotQa}" /></td>
                  <td><img src="{plotBa}" /></td>
                </tr>
            """
    html += "</tbody> </table> </body> </html>"

    # Write to file
    with open(f"{out_dir}/ips_poolseq_report.html", "w", encoding="utf-8") as f:
        f.write(html)
    
def main ():
    #Create output directories if they doesn't exist
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(f"{out_dir}/after", exist_ok=True)
    os.makedirs(f"{out_dir}/before", exist_ok=True)

    #Parse Fastp JSON report
    json_files = glob.glob(json_list)
    data = parse_json(json_files)
    write_tsv_report(data)
    write_html_report()

if __name__ == "__main__":
    main()

