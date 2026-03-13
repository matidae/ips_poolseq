#!/usr/bin/env python3

#----------------------------------------------------------------------
# Plots allele frequency trajectories for the top SNPs from each category
# (fluctuating, directional, mixed, drift).
#
# Inputs:
#   - {thinned_dir}/tests.{prefix}.FDR_{category}.thinned.tsv
#   - {work_dir}/z_year.{prefix}.tsv
#
# Output (AF and Z mode):
#   - {prefix}_trajectories_Z.png
#   - {prefix}_trajectories_delta_Z.png
#   - SNPs_Z_trajectories.html
#----------------------------------------------------------------------

import os
import math
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

work_dir    = "../results/08_models"
thinned_dir = "../results/08_models/s2_thinning"
out_dir     = "../results/08_models/trajectories"

top_n_snps = 10
category_colors = {"fluctuating": "#E07B39", "directional": "#3A7DBF",
    "mixed": "#6A0DAD", "drift": "#4A4A4A",}

categories = ["fluctuating", "directional", "mixed", "drift"]

# Get list of prefixes from thinned_dir
def get_prefixes(thinned_dir):
    files = glob.glob(f"{thinned_dir}/tests.*.FDR.*.thinned.tsv")
    return sorted(set(os.path.basename(f).split(".")[1] for f in files))

# Reads z_year file and returns trajectories for all valid SNPs (applies minor allele polarization).
# - mode="Z":  returns z-transformed values
# - mode="AF": back-transforms z to p
def read_z_trajectories(z_file, mode):
    trajectories = {}
    with open(z_file) as f:
        header    = f.readline().strip().split("\t")
        year_cols = header[4:]
        years     = [int(col.split("_")[-1]) for col in year_cols]

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

            # Minor allele polarisation
            if sum(zlist) / len(zlist) > math.pi / 2:
                zlist = [math.pi - z for z in zlist]

            if mode == "AF":
                trajectories[key] = [(math.sin(z / 2)) ** 2 for z in zlist]
            else:
                trajectories[key] = zlist

    all_years = list(range(min(years), max(years) + 1))
    return trajectories, years, all_years


# Plotting function for one prefix and one plot type (raw or delta)
def plot_one(prefix, snps_dict, trajectories_dict, years, all_years,
             plot_type, out_dir, snps_dict_bg, trajectories_dict_bg, categories,
             mode):
    is_delta = plot_type.startswith("delta")

    if mode == "AF":
        suptitle_base = "Trajectory of SNPs allele frequencies"
        ylabel = "Raw allele frequency" if not is_delta else "delta  allele frequency (from t=0)"
    else:
        suptitle_base = "Trajectory of z-transformed allele frequencies"
        ylabel = "z-transformed AF" if not is_delta else "delta  z-transformed AF (from t=0)"

    suffix = " (delta from initial timepoint)" if is_delta else ""
    n_cats = len(categories)

    fig, axes = plt.subplots(n_cats, 1, figsize=(14, 5 * n_cats), sharey=True, sharex=True)
    if n_cats == 1:
        axes = [axes]

    fig.suptitle(f"{prefix} - {suptitle_base}{suffix}", fontsize=13, fontweight="bold")

    for ax, category in zip(axes, categories):
        snps = snps_dict.get(category, [])
        trajectories = trajectories_dict.get(category, {})
        bg_snps = snps_dict_bg.get(category, [])
        bg_traj = trajectories_dict_bg.get(category, {})
        color = category_colors[category]

        base_rgb = np.array(matplotlib.colors.to_rgb(color))
        white = np.array([1, 1, 1])
        colors = [tuple(white * (1 - t) + base_rgb * t) for t in np.linspace(0.4, 1.0, max(len(snps), 1))]

        # Background SNPs
        for (key, _) in bg_snps:
            vals = bg_traj[key]
            if is_delta:
                vals = [v - vals[0] for v in vals]
            ax.plot(years, vals, color="grey", linewidth=0.4, alpha=0.2, rasterized=True)

        # Top SNPs
        for (key, _), c in zip(snps, colors):
            vals = trajectories[key]
            if is_delta:
                vals = [v - vals[0] for v in vals]
            ax.plot(years, vals, color=c, linewidth=1.5, marker="o", markersize=4,
                    markerfacecolor="white", markeredgewidth=1.2, markeredgecolor=c, alpha=0.85)

        # X axis
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


# HTML report

TABLE_HTML = """
<table>
    <thead>
        <tr><th>Category</th><th>LRT1</th><th>LRT2</th><th>Description</th></tr>
    </thead>
    <tbody>
        <tr><td style="color:#7A7A7A;font-weight:600;">Drift</td><td>ns</td><td>ns</td><td>No departure from drift</td></tr>
        <tr><td style="color:#3A7DBF;font-weight:600;">Directional</td><td>sig</td><td>ns</td><td>Linear trend, no oscillation</td></tr>
        <tr><td style="color:#E07B39;font-weight:600;">Fluctuating</td><td>ns</td><td>sig</td><td>Oscillation, no net directional trend</td></tr>
        <tr><td style="color:#6A0DAD;font-weight:600;">Mixed</td><td>sig</td><td>sig</td><td>Directional trend + oscillation</td></tr>
    </tbody>
</table>
<p class="description2">FDR p-value: 0.05 — LRT1: directional vs drift — LRT2: saturated vs drift</p>
"""

CSS = """
body { font-family: sans-serif; font-size: 14px; color: #222; padding: 24px; }
h1   { margin-bottom: 4px; }
h2   { margin: 24px 0 12px; }
h3   { margin-bottom: 8px; border-bottom: 1px solid #eee; padding-bottom: 4px; }
.description  { color: #555; font-size: 0.92rem; line-height: 1.6; margin-bottom: 12px; }
.description2 { color: #555; font-size: 0.80rem; line-height: 1.6; margin-bottom: 12px; }
.plot-row { display: grid; grid-template-columns: 1fr 1fr; gap: 24px; }
img { width: 90%; display: block; transition: transform 0.2s ease;
      transform-origin: top left; position: relative; z-index: 1; }
hr  { border: none; border-top: 1px solid #eee; margin: 32px 0; }
table { border-collapse: collapse; margin-bottom: 12px; }
th, td { padding: 6px 16px; text-align: left; border: 1px solid #eee; font-size: 12px; }
th { background: #f5f5f5; font-weight: 600; }
"""


def generate_html(prefixes, out_dir, mode):
    if mode == "AF":
        page_title = "SNPs allele frequency trajectories"
        h1 = "SNPs allele frequency trajectories"
        intro = ("Allele frequency trajectories per selection category — "
                       "top 10 SNPs by p-value are highlighted and the next 1000 are in grey.")
        col_left = ("Delta allele frequency (from t=0)",
                       "delta p = p<sub>t</sub> − p<sub>0</sub>. Removes the effect of initial "
                       "allele frequency, making trajectories directly comparable "
                       "(ranked by FDR-adjusted p-value).")
        col_right  = ("Allele frequency",
                       "Raw allele frequency trajectories (p) across all years per selection "
                       "category (ranked by FDR-adjusted p-value).")
        img_suffixes = ("delta_AFs", "AFs")
        out_file    = os.path.join(out_dir, "SNPs_trajectories.html")
    else:
        page_title = "SNPs z-transformed allele frequency trajectories"
        h1 = "SNPs z-transformed allele frequency trajectories"
        intro = ("z-transformed allele frequency trajectories per selection category — "
                       "top 10 SNPs by p-value are highlighted and the next 1000 are in grey.")
        col_left = ("Delta of z-transformed allele frequency (from t=0)",
                       "delta z = z<sub>t</sub> − z<sub>0</sub>. Removes the effect of initial "
                       "allele frequency, making trajectories directly comparable "
                       "(ranked by FDR-adjusted p-value).")
        col_right = ("z-transformed allele frequency",
                       "Raw z-transformed allele frequency trajectories across all years per "
                       "selection category (ranked by FDR-adjusted p-value).")
        img_suffixes = ("delta_Z", "Z")
        out_file = os.path.join(out_dir, "SNPs_Z_trajectories.html")

    header_row = f"""
    <div class="plot-row">
        <div class="col-header">
            <h3>{col_left[0]}</h3>
            <p class="description">{col_left[1]}</p>
        </div>
        <div class="col-header">
            <h3>{col_right[0]}</h3>
            <p class="description">{col_right[1]}</p>
        </div>
    </div>
    """

    sections = ""
    for prefix in prefixes:
        delta_img = f"{prefix}_trajectories_{img_suffixes[0]}.png"
        raw_img = f"{prefix}_trajectories_{img_suffixes[1]}.png"
        sections += f"""
        <div class="prefix-section">
            <h2>{prefix}</h2>
            <div class="plot-row">
                <div class="img-wrapper"><img src="{delta_img}" alt="{prefix} delta"></div>
                <div class="img-wrapper"><img src="{raw_img}" alt="{prefix} raw"></div>
            </div>
        </div>
        <hr>
        """

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <link rel="icon" type="image/png" href="../favicon_sbb.png">
    <meta charset="UTF-8">
    <title>{page_title}</title>
    <style>{CSS}</style>
</head>
<body>
    <h1>{h1}</h1>
    <br>
    {TABLE_HTML}
    <br>
    <p class="description">{intro}</p>
    {header_row}
    {sections}
</body>
</html>"""

    with open(out_file, "w") as f:
        f.write(html)


# Main pipeline
def run_mode(mode, prefixes):
    # Map mode  plot_type pair (raw, delta)
    if mode == "AF":
        plot_types = ("AFs", "delta_AFs")
    else:
        plot_types = ("Z", "delta_Z")

    for prefix in prefixes:
        z_file = f"{work_dir}/z_year.{prefix}.tsv"
        trajectories_all, years, all_years = read_z_trajectories(z_file, mode)

        snps_dict, snps_dict_bg = {}, {}
        trajectories_dict, traj_dict_bg = {}, {}

        for category in categories:
            thinned_file = f"{thinned_dir}/tests.{prefix}.FDR.{category}.thinned.tsv"
            snps, snps_bg = [], []

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
            trajectories_dict[category] = {k: trajectories_all[k] for k, _ in snps}
            traj_dict_bg[category] = {k: trajectories_all[k] for k, _ in snps_bg}

        for plot_type in plot_types:
            plot_one(prefix, snps_dict, trajectories_dict, years, all_years,
                     plot_type, out_dir, snps_dict_bg, traj_dict_bg, categories, mode)

    generate_html(prefixes, out_dir, mode)


def main():
    os.makedirs(out_dir, exist_ok=True)
    sns.set_style("whitegrid")
    prefixes = get_prefixes(thinned_dir)
    #Modes are AFs types: p or z
    modes = ["AF", "Z"]
    for mode in modes:
        run_mode(mode, prefixes)

if __name__ == "__main__":
    main()