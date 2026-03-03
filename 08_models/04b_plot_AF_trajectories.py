#!/usr/bin/env python3

#----------------------------------------------------------------------
# Plots allele frequency trajectories for the top SNPs from each category 
# (fluctuating, directional, mixed, drift).
#
# Inputs:
#   - {thinned_dir}/tests.{prefix}.FDR_{category}.thinned.tsv
#   - {work_dir}/z_year.{prefix}.tsv
#
# Output:
#   - {out_dir}/{prefix}_trajectories_AFs.png
#   - {out_dir}/{prefix}_trajectories_delta_AFs.png
#----------------------------------------------------------------------

import os
import math
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

work_dir     = "../results/08_models"
thinned_dir  = "../results/08_models/s2_thinning"
out_dir      = "../results/08_models/plots/trajectories"

top_n_snps = 10
category_colors = {
    "fluctuating": "#E07B39", "directional": "#3A7DBF",
    "mixed":       "#6A0DAD", "drift": "#4A4A4A",
}

categories = ["fluctuating", "directional", "mixed", "drift"]

def get_prefixes(thinned_dir):
    files = glob.glob(f"{thinned_dir}/tests.*.FDR.*.thinned.tsv")
    prefixes = sorted(set([os.path.basename(i).split(".")[1] for i in files]))    
    return prefixes

#Reads z files and returns raw allele frequencies for requested SNPs.
#Back-transforms Z to AF: p = sin(z/2)^2
#Applies minor allele polarization (mean Z > pi/2 -> flip)
def read_z_trajectories(z_file):
    trajectories = {}
    years        = []
    with open(z_file) as f:
        header     = f.readline().strip().split("\t")
        year_cols  = header[4:]
        # Extract year from prefix
        years = [int(col.split("_")[-1]) for col in year_cols]        

        for line in f:
            cols = line.strip().split("\t")
            key  = cols[0] + "_" + cols[1]

            zlist = []
            for j in range(4, len(cols)):
                val = cols[j].split(",")[0]
                if val == "NA":
                    zlist.append(None)
                else:
                    zlist.append(float(val))

            if None in zlist:
                continue

            # Minor allele polarization
            zbar = sum(zlist) / len(zlist)
            if zbar > math.pi / 2:
                zlist = [math.pi - z for z in zlist]

            # Back-transform to allele frequency
            af = [(math.sin(z / 2)) ** 2 for z in zlist]
            trajectories[key] = af
    # All years including those without data from min to max for x axis spacing
    all_years = list(range(min(years), max(years) + 1))    

    return trajectories, years, all_years


def plot_one(prefix, snps_dict, trajectories_dict, years, all_years,
             plot_type, out_dir, snps_dict_bg, trajectories_dict_bg, categories):
    n_cats = len(categories)
    fig, axes = plt.subplots(n_cats, 1, figsize=(14, 5 * n_cats), sharey=True, sharex=True)
    if n_cats == 1:
        axes = [axes]

    title_suffix = "" if plot_type == "AFs" else "(delta from initial timepoint)"
    fig.suptitle(f"{prefix} - Trajectory of SNPs allele frequencies {title_suffix}",
                 fontsize=13, fontweight="bold")

    ylabel = "Raw allele frequency" if plot_type == "AFs" else "delta allele frequency (from t=0)"

    for ax, category in zip(axes, categories):
        snps         = snps_dict.get(category, [])
        trajectories = trajectories_dict.get(category, {})
        color        = category_colors[category]

        base_rgb = np.array(matplotlib.colors.to_rgb(color))
        white    = np.array([1, 1, 1])

        colors = [tuple(white * (1 - t) + base_rgb * t)
            for t in np.linspace(0.4, 1.0, max(len(snps), 1))]
        
        # Background - next 1000 SNPs in grey
        bg_snps = snps_dict_bg.get(category, [])
        bg_trajectories = trajectories_dict_bg.get(category, {})

        for (key, _) in bg_snps: 
            vals_bg = bg_trajectories[key]
            if plot_type == "delta_AFs":
                vals_bg = [v - vals_bg[0] for v in vals_bg]
            ax.plot(years, vals_bg, color="grey", linewidth=0.4, alpha=0.2, rasterized=True)
            

        for (key, _), c in zip(snps, colors):
            vals = trajectories[key]
            if plot_type == "delta_AFs":
                vals = [v - vals[0] for v in vals]

            # Plot at real year positions
            ax.plot(years, vals, color=c, linewidth=1.5, marker="o", markersize=4,
                    markerfacecolor="white", markeredgewidth=1.2, markeredgecolor=c, alpha=0.85)

        # x_axis - all years as ticks, only data years labeled
        ax.set_xticks(all_years)
        if ax == axes[-1]:
            labels = [str(y) if y in years else "" for y in all_years]
            ax.set_xticklabels(labels, ha="right", fontsize=12)
            ax.tick_params(axis="both", labelsize=12)            
            ax.set_xlabel("Year", fontsize=13)
        else:
            ax.set_xticklabels([])

        ax.set_xlim(all_years[0] - 0.3, all_years[-1] + 0.5)
        ax.set_title(category.capitalize(), fontsize=11, fontweight="bold", color=color, loc="left")
        fig.text(0.04, 0.5, ylabel, va="center", rotation="vertical", fontsize=13)
        sns.despine(ax=ax)

    plt.tight_layout()
    plt.subplots_adjust(left=0.08)  
    out_path = os.path.join(out_dir, f"{prefix}_trajectories_{plot_type}.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()

def generate_html(prefixes, out_dir):
    desc_AF = "Raw allele frequency trajectories of SNPs per selection category."

    desc_DELTA = "Allele frequency change relative to the initial timepoint. \
        Removes the effect of initial allele frequency making trajectories directly comparable."
    
    table_html = f"""
    <table>
        <thead>
            <tr>
                <th>Category</th>
                <th>LRT1</th>
                <th>LRT2</th>
                <th>Description</th>
            </tr>
        </thead>
        <tbody>
            <tr><td style="color:#7A7A7A; font-weight:600;">Drift</td><td>ns</td><td>ns</td><td>No departure from drift</td></tr>
            <tr><td style="color:#3A7DBF; font-weight:600;">Directional</td><td>sig</td><td>ns</td><td>Linear trend, no oscillation</td></tr>
            <tr><td style="color:#E07B39; font-weight:600;">Fluctuating</td><td>ns</td><td>sig</td><td>Oscillation, no net directional trend</td></tr>
            <tr><td style="color:#6A0DAD; font-weight:600;">Mixed</td><td>sig</td><td>sig</td><td>Directional trend + oscillation</td></tr>
        </tbody>
    </table>
    <p class="description2">
        FDR p-value: 0.05 -
        LRT1: LRT directional vs drift -
        LRT2: LRT saturated vs drift        
    </p>
    """
    
    header_row = f"""
    <div class="plot-row">
        <div class="col-header">
            <h3>Delta allele frequency (from t=0)</h3>
            <p class="description">{desc_DELTA.strip()}</p>
        </div>
        <div class="col-header">
            <h3>Allele frequency</h3>
            <p class="description">{desc_AF.strip()}</p>
        </div>
    </div>
    """    

    sections = ""
    for prefix in prefixes:
        af_img    = f"{prefix}_trajectories_AFs.png"
        delta_img = f"{prefix}_trajectories_delta_AFs.png"

        sections += f"""
        <div class="prefix-section">
            <h2>{prefix}</h2>
            <div class="plot-row">
                <div class="img-wrapper">
                    <img src="{delta_img}" alt="{prefix} delta">
                </div>
                <div class="img-wrapper">
                    <img src="{af_img}" alt="{prefix} AF">
                </div>
            </div>
        </div>
        <hr>
        """

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <link rel="icon" type="image/png" href="../favicon_sbb.png">
    <meta charset="UTF-8">
    <title>SNPs allele frequency trajectory plots</title>
    <style>
    body {{
        font-family: sans-serif;
        font-size: 14px;
        color: #222;
        padding: 24px;
    }}

    h1 {{ margin-bottom: 4px; }}
    h2 {{ margin: 24px 0 12px; }}
    h3 {{ margin-bottom: 8px; border-bottom: 1px solid #eee; padding-bottom: 4px; }}

    .description {{
        color: #555;
        font-size: 0.92rem;
        line-height: 1.6;
        margin-bottom: 12px;
    }}
        .description2 {{
        color: #555;
        font-size: 0.8rem;
        line-height: 1.6;
        margin-bottom: 12px;
    }}


    .plot-row {{
        display: grid;
        grid-template-columns: 1fr 1fr;
        gap: 24px;
    }}

    img {{
        width: 90%;
        display: block;
        cursor: zoom-in;
        transition: transform 0.2s ease;
        transform-origin: top left;
        position: relative;
        z-index: 1;
    }}

    img:hover {{
        transform: scale(1.15);
        box-shadow: 0 8px 32px rgba(0,0,0,0.2);
        z-index: 999;
        cursor: zoom-out;
    }}

    hr {{ border: none; border-top: 1px solid #eee; margin: 32px 0; }}
    table {{ border-collapse: collapse; margin-bottom: 12px; }}
    th, td {{ padding: 6px 16px; text-align: left; border: 1px solid #eee; font-size: 12px; }}
    th {{ background: #f5f5f5; font-weight: 600; }}
</style>
</head>
<body>
    <h1>SNPs allele frequency trajectories</h1>
    <br>
    {table_html}
    <br>
    <p class="description">
        Allele frequency trajectories per selection category - top 10 SNPs by p-value are highlighted and the  next 1000 in grey.       
    </p>       
    {header_row}
    {sections}
</body>
</html>"""

    out_file = os.path.join(out_dir, "trajectories_report.html")
    with open(out_file, "w") as f:
        f.write(html)

def main():
    os.makedirs(out_dir, exist_ok=True)
    sns.set_style("whitegrid")
    prefixes = get_prefixes(thinned_dir)

    for prefix in prefixes:  
        z_file = f"{work_dir}/z_year.{prefix}.tsv"

        # Load all valid trajectories once for this prefix
        trajectories_all, years, all_years = read_z_trajectories(z_file)

        snps_dict = {}
        snps_dict_bg = {}
        trajectories_dict = {}
        trajectories_dict_bg = {}
        
        # For each category, select first TOP_N valid SNPs
        for category in categories:
            thinned_file = f"{thinned_dir}/tests.{prefix}.FDR.{category}.thinned.tsv"

            snps = []
            snps_bg = []

            with open(thinned_file) as f:
                for line in f:
                    cols = line.strip().split("\t")
                    key = cols[1] + "_" + cols[2]
                    if key not in trajectories_all:
                        continue
                    if len(snps) < top_n_snps:
                        snps.append((key, key))
                    elif len(snps_bg) < 1000:
                        snps_bg.append((key, key))
                    if len(snps) == top_n_snps and len(snps_bg) == 1000:
                        break
            snps_dict[category] = snps
            snps_dict_bg[category] = snps_bg

            trajectories_dict[category] = {key: trajectories_all[key] for key, _ in snps}
            trajectories_dict_bg[category] = {key: trajectories_all[key] for key, _ in snps_bg}

        # Generate plots
        plot_one(prefix, snps_dict, trajectories_dict, years, all_years,
                 "AFs", out_dir, snps_dict_bg, trajectories_dict_bg, categories)

        plot_one(prefix, snps_dict, trajectories_dict, years, all_years,
                 "delta_AFs", out_dir, snps_dict_bg, trajectories_dict_bg, categories)
    #Generate HTML report    
    generate_html(prefixes, out_dir)

if __name__ == "__main__":
    main()