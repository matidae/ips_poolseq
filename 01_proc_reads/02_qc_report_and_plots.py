"""
Reads fastp JSON output files and generate summary tables and plots.

Input:
    - fastp JSON report files (e.g., *.json)

Output:
    - all_poolseq_report.tsv : summary table with quality metrics for all samples
    - $prefix_quality.png : quality score plot for each sample
    - $prefix_base_content.png : Base composition plot for each sample
"""

import os
import glob
import json
import matplotlib.pyplot as plt

in_dir = "../data/01_proc_reads"
out_dir = "../results/01_proc_reads/"

ips_genome_size = 224977219

# Input files
qubit_in = f"{in_dir}/sample_qubit"
json_list = f"{in_dir}/*/*.json"

# Output files
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
   plot_filename = out_dir + filter + "/" + filename.replace('.fq.gz','.qual.png')
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
   plot_filename = out_dir + filter + "/" + filename.replace('.fq.gz','.base.png')
   plt.savefig(plot_filename, dpi=300)
   plt.close()


def main ():
    #Creat output directories if they doesn't exist
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(f"{out_dir}/after", exist_ok=True)
    os.makedirs(f"{out_dir}/before", exist_ok=True)

    json_files = glob.glob(json_list)

    data = {}    
    col_names = ["Idn     ", \
                "After", "Gbp", "Cov", "Length", "GC",\
                "Insert", "Dups", \
                "Before",  "Lost", "LQual", "Short", "Qubit", "%Reads_lost", "GCbefore", "File"]
    header = "\t".join(col_names)

    sample_qubit = {i.rstrip().split()[1] : i.split()[0] for i in open(qubit_in) }
    
    for file in json_files:
        with open(file, 'r') as f:
            fastp_out = json.load(f)
            
            prefix = file.split("/")[-2]
            print(prefix)
            prefix_short = prefix[:12]
            idn = "_".join(prefix.split("_")[:3])      
            
            reads_after = fastp_out["summary"]["after_filtering"]["total_reads"]
            total_bp = fastp_out["summary"]["after_filtering"]["total_bases"]
            coverage = total_bp/ips_genome_size
            length = (fastp_out["summary"]["after_filtering"]["read1_mean_length"] + fastp_out["summary"]["after_filtering"]["read2_mean_length"])*0.5
            gc_content = fastp_out["summary"]["after_filtering"]["gc_content"]
            
            dups = fastp_out["duplication"]["rate"]

            reads_before = fastp_out["summary"]["before_filtering"]["total_reads"]
            reads_lost = reads_before - reads_after
            low_qual = fastp_out["filtering_result"]["low_quality_reads"]
            too_short = fastp_out["filtering_result"]["too_short_reads"]
            
            filename = file.split('/')[-1].replace("qc_report.json", "fq.gz")
            quality_curves_after = fastp_out["read1_after_filtering"]["quality_curves"]["mean"]
            content_curves_after = fastp_out["read1_after_filtering"]["content_curves"]
           # plot_quality(quality_curves_after, filename, "after")
            #plot_base_content(content_curves_after, filename, "after")

            quality_curves_before = fastp_out["read1_before_filtering"]["quality_curves"]["mean"]
            content_curves_before = fastp_out["read1_before_filtering"]["content_curves"]
            #plot_quality(quality_curves_before, filename, "before")
            #plot_base_content(content_curves_before, filename, "before")

            qubit = sample_qubit[prefix_short]
            gc_content_before = fastp_out["summary"]["before_filtering"]["gc_content"]            
            pct_reads_lost = round(100 - (reads_after/reads_before)*100, 1)
            
            row = [idn, \
                reads_after, total_bp, coverage, length, gc_content,\
                dups, \
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
    
    data_sort = sorted(data.items(), key=lambda x: x[1][0])
    data_sort_new=[]

    report_fh = open(report_out, 'w')
    report_fh.write(header + '\n')
    for k, v in data_sort:

        elems = [
            v[0],                                # idn
            round(v[1]/1e6, 1),                # reads_after (M)
            round(v[2]/1e9, 1),                # total_bp (Gb)
            round(v[3], 1),                      # coverage
            int(v[4]),                           # length
            round(v[5] * 100, 1),                # gc_content %            
            round(v[6] * 100, 1),                # dups %            
            round(v[7]/1e6, 1),                # reads_before (M)
            round(v[8]/1e6, 1),               # reads_lost (M)
            round(v[9]/1e6, 1),               # low_qual (M)
            round(v[10]/1e6, 1),               # too_short (M)
            v[11],                               # qubit
            v[12],                               # pct_reads_lost
            round(v[13] * 100, 1),               # gc_content_before %
            v[14]
        ]            
        elems_str = list(map(str, elems))
        data_sort_new.append(elems)
        report_fh.write("\t".join(elems_str) + "\n")

    report_fh.close()

if __name__ == "__main__":
    main()

