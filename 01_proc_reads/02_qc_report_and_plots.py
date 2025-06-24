"""
Using the JSON output of Fastp produces a table and plots for each sample
In: json fastp outputs
Out: all_poolseq_report.tsv table with summary of all samples
     png plots for quality and base content for each sample
"""
import glob
import json
import matplotlib.pyplot as plt

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

out_dir = "../../results/01_proc_reads/"
json_files = glob.glob(out_dir + "json_reports/*.json")

data = {}
ips_genome_size = 224977219
col_names = ["Idn     ", "Year", "Time", "Country", "Region", "Reps", \
             "After", "Gbp", "Cov", "Length", "GC",\
             "Insert", "Dups", "q20A", \
              "Before",  "Lost", "LQual", "Ns", "Short", "q20B", "Qubit", "%Reads_lost", "GCbefore", "File"]
header = "\t".join(col_names)

sample_qubit = {i.rstrip().split()[1] : i.split()[0] for i in open(out_dir + "sample_qubit") }

for file in json_files:
    with open(file, 'r') as f:
        fastp_out = json.load(f)        
        prefix = file.split("/")[5].replace(".qc_report.json", "")        
        prefix_short = prefix[:12]
        idn = "_".join(prefix.split("_")[:3])        
               
        year = prefix.split("_")[2]
        time = prefix.split("_")[1][0]
        region = prefix.split("_")[0][0]
        country = prefix.split("_")[0][1:]
        replicate = prefix.split("_")[1][1]
        
        reads_after = fastp_out["summary"]["after_filtering"]["total_reads"]
        total_bp = fastp_out["summary"]["after_filtering"]["total_bases"]
        coverage = total_bp/ips_genome_size
        length = (fastp_out["summary"]["after_filtering"]["read1_mean_length"] + fastp_out["summary"]["after_filtering"]["read2_mean_length"])*0.5
        gc_content = fastp_out["summary"]["after_filtering"]["gc_content"]
        
        insert = fastp_out["insert_size"]["peak"]
        dups = fastp_out["duplication"]["rate"]
        over_q20 = fastp_out["summary"]["after_filtering"]["q20_rate"]
        before_q20 = fastp_out["summary"]["before_filtering"]["q20_rate"]

        reads_before = fastp_out["summary"]["before_filtering"]["total_reads"]
        reads_lost = reads_before - reads_after
        low_qual = fastp_out["filtering_result"]["low_quality_reads"]
        many_Ns = fastp_out["filtering_result"]["too_many_N_reads"]
        too_short = fastp_out["filtering_result"]["too_short_reads"]
        
        filename = file.split('/')[5].replace("qc_report.json", "fq.gz")
        quality_curves_after = fastp_out["read1_after_filtering"]["quality_curves"]["mean"]
        content_curves_after = fastp_out["read1_after_filtering"]["content_curves"]
        plot_quality(quality_curves_after, filename, "after")
        plot_base_content(content_curves_after, filename, "after")

        quality_curves_before = fastp_out["read1_before_filtering"]["quality_curves"]["mean"]
        content_curves_before = fastp_out["read1_before_filtering"]["content_curves"]
        plot_quality(quality_curves_before, filename, "before")
        plot_base_content(content_curves_before, filename, "before")

        qubit = sample_qubit[prefix_short]
        gc_content_before = fastp_out["summary"]["before_filtering"]["gc_content"]
        total_bp_before = fastp_out["summary"]["before_filtering"]["total_bases"]
        pct_reads_lost = round(100 - (reads_after/reads_before)*100, 1)
        
        row = [idn, year, time, country, region, replicate, \
               reads_after, total_bp, coverage, length, gc_content,\
               insert, dups, over_q20,\
               reads_before, reads_lost, low_qual, many_Ns, too_short, before_q20, qubit, pct_reads_lost, gc_content_before, filename]
               
        if idn not in data.keys():
            data[idn] = row            
        else:
            #sum:  "reads_after", "total_bp", "reads_before", "reads_lost", "low_qual", "many_Ns", "too_short"
            index_to_sum = [6, 7, 14, 15, 16, 17, 18]
            for k in index_to_sum:
                data[idn][k] += row[k]
            #Calculate coverage            
            coverage = data[idn][7]/ips_genome_size
            data[idn][8] = coverage            
            #Append filename
            data[idn][-1] += "," + row[-1]
            #Reads lost
            pct_reads_lost = round(100 - (data[idn][6]/data[idn][14])*100, 1)
            data[idn][21] = pct_reads_lost

data_sort = sorted(data.items(), key=lambda x: x[1][1:6])
data_sort_new=[]
out = open(out_dir + 'all_poolseq_report.tsv', 'w')
out.write(header + '\n')
for k, v in data_sort:    
    after =  round(v[6]/1e6, 1)
    elems=v[:6]+ [round(v[6]/1e6, 1) ,round(v[7]/1e9, 1), round(v[8], 1), int(v[9]), round(v[10]*100,1), \
          int(v[11]), round(v[12]*100,1), round(v[13]*100,1), \
            round(v[14]/1e6, 1), round(v[15]/1e6, 1),round(v[16]/1e6, 1),round(v[17]/1e6, 1), round(v[18]/1e6, 1), \
                round(v[19]*100,1), v[20], v[21], round(v[22]*100, 1), v[23]]
    elems_str = list(map(str, elems))
    data_sort_new.append(elems)
    out.write("\t".join(elems_str) + "\n")

out.close()
