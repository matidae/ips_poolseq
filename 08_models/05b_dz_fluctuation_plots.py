#!/usr/bin/env python3

#----------------------------------------------------------------------
# Produces multiple plots from the dz analysis of fluctuating SNPs 
#
# Inputs:
#   - {work_dir}/{prefix}_interval_types.tsv
#   - {work_dir}/{prefix}_mean_deltaZ_ttests.tsv
#   - {work_dir}/{prefix}_correlations.tsv
#   - {work_dir}/{prefix}_deltaZ_per_snp.tsv
#   - {work_dir}/{prefix}_sync_scores.tsv
# Outputs:
#   - {prefix}_01_mean_dz_per_interval.png
#   - {prefix}_02_correlation_heatmap.png
#   - {prefix}_03_dz_distribution_per_interval.png
#   - {prefix}_04_snp_heatmap.png
#   - {prefix}_07_dz_scatterplots.png
#   - {prefix}.html
#----------------------------------------------------------------------

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import seaborn as sns
from scipy.stats import pearsonr, linregress

work_dir    = "../results/08_models/analysis"
out_dir     = "../results/08_models/plots/deltaZ"
DPI         = 150
PVAL_THRESH = 0.05
COLOR_G     = "#E07B39"
COLOR_Y     = "#3A7DBF"
UNIQUE_BBOX = dict(facecolor="#FFCC80", edgecolor="#E07B39", pad=3)
TOP_N    = 10
ALPHA_BG = 0.08

# Get list of population prefixes
def get_prefixes(work_dir):
    files = glob.glob(f"{work_dir}/*.fluctuating_interval_types.tsv")
    return sorted(os.path.basename(f).split(".")[0] for f in files)

# Load files needed for plotting
def load_files(work_dir, prefix):
    base           = f"{work_dir}/{prefix}.fluctuating"
    interval_types = pd.read_csv(f"{base}_interval_types.tsv", sep="\t")
    ttests         = pd.read_csv(f"{base}_mean_deltaZ_ttests.tsv", sep="\t")
    correlations   = pd.read_csv(f"{base}_correlations.tsv", sep="\t", index_col=0)
    dz_per_snp     = pd.read_csv(f"{base}_deltaZ_per_snp.tsv", sep="\t", index_col=0)
    sync_scores    = pd.read_csv(f"{base}_sync_scores.tsv", sep="\t")
    return interval_types, ttests, correlations, dz_per_snp, sync_scores

# Highlights x tick labels not in shared_intervals with yellow background.
def highlight_unique_xticklabels(ax, shared_intervals):    
    ax.figure.canvas.draw()
    for label in ax.get_xticklabels():
        if label.get_text() and label.get_text() not in shared_intervals:
            label.set_bbox(UNIQUE_BBOX)

#Highlights y tick labels not in shared_intervals with yellow background.
def highlight_unique_yticklabels(ax, shared_intervals):    
    ax.figure.canvas.draw()
    for label in ax.get_yticklabels():
        if label.get_text() and label.get_text() not in shared_intervals:
            label.set_bbox(UNIQUE_BBOX)

#Barplot of mean_dz for each interval, colored by up/down type, with * for significant.
def plot_mean_dz_per_interval(ttests, interval_types, prefix, out_dir, shared_intervals):
    fig, ax = plt.subplots(figsize=(12, 6))
    colors = [COLOR_G if t == "up" else COLOR_Y for t in interval_types["type"]]
    ax.bar(ttests["interval"], ttests["mean_dz"],
           color=colors, edgecolor="white", linewidth=0.8, alpha=0.85)
    
    for i, (_, row) in enumerate(ttests.iterrows()):
        y_offset = 0.002 if row["mean_dz"] >= 0 else -0.004
        sig = "*" if row["p_value"] < PVAL_THRESH else ""
        if sig:
            ax.text(i, row["mean_dz"] + y_offset, sig, ha="center", va="bottom", fontsize=12, fontweight="bold")
    
    ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
    ax.set_xlabel("Interval", fontsize=14)
    ax.set_ylabel("Mean dz", fontsize=14)
    ax.tick_params(axis="both", labelsize=12)    
    sns.despine(ax=ax)
    plt.tight_layout()
    highlight_unique_xticklabels(ax, shared_intervals)
    out_path = os.path.join(out_dir, f"{prefix}_01_mean_dz_per_interval.png")
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close()

# Heatmap of correlations between per-SNP dz values across intervals, and masking of adjacent intervals
def plot_correlation_heatmap(correlations, dz_per_snp, prefix, out_dir, shared_intervals):
    corr = correlations.astype(float)
    n = corr.shape[0]

    # Build mask: upper triangle + diagonal + adjacent
    mask = np.triu(np.ones_like(corr, dtype=bool), k=0)
    for i in range(n - 1):
        mask[i+1, i] = True

    # Drop rows and columns that are entirely masked (no visible data)
    visible_rows = [i for i in range(n) if i >= 2]
    visible_cols = [j for j in range(n) if j <= n - 3]
    corr_trimmed = corr.iloc[visible_rows, visible_cols]
    mask_trimmed = mask[np.ix_(visible_rows, visible_cols)]


    # Set adjacent cells (orig_i == orig_j + 2) to NaN before plotting
    row_intervals = corr_trimmed.index.tolist()
    col_intervals = corr_trimmed.columns.tolist()
    orig_rows = corr.index.tolist()
    orig_cols = corr.columns.tolist() 

   # Add adjacent cells to mask_trimmed
    for i, row_iv in enumerate(row_intervals):
        for j, col_iv in enumerate(col_intervals):
            orig_i = orig_rows.index(row_iv)
            orig_j = orig_cols.index(col_iv)
            if orig_i == orig_j + 2:
                mask_trimmed[i, j] = True

    # Now compute ns_cells on truly visible cells only
    ns_cells = []
    for i, row_iv in enumerate(row_intervals):
        for j, col_iv in enumerate(col_intervals):
            if mask_trimmed[i, j]:
                continue
            r, p = pearsonr(dz_per_snp[row_iv], dz_per_snp[col_iv])
            if p >= 0.05:
                corr_trimmed.iloc[i, j] = np.nan
                ns_cells.append((i, j))

    fig, ax = plt.subplots(figsize=(6, 4))
    sns.heatmap(corr_trimmed, mask=mask_trimmed, annot=True, fmt=".2f", center=0,
                cmap="RdBu_r", linewidths=0.5, linecolor="white",
                vmin=-1, vmax=1, ax=ax, annot_kws={"size": 8}, 
                cbar_kws={"label": "", "shrink": 0.9})
    ax.collections[0].colorbar.ax.tick_params(labelsize=7)

    # Annotate adjacent cells with NA
    for i, row_iv in enumerate(row_intervals):
        for j, col_iv in enumerate(col_intervals):
            orig_i = orig_rows.index(row_iv)
            orig_j = orig_cols.index(col_iv)
            if orig_i == orig_j + 2:
                ax.add_patch(Rectangle((j, i), 1, 1, facecolor="#ffffff", edgecolor="none", zorder=0))
                ax.text(j + 0.5, i + 0.5, "NA", ha="center", va="center", fontsize=8, color="#333333")
    # Annotate ns cells
    for i, j in ns_cells:
        ax.add_patch(Rectangle((j, i), 1, 1, facecolor="#e0e0e0", edgecolor="none", zorder=1))
        ax.text(j + 0.5, i + 0.5, "ns", ha="center", va="center", fontsize=8, color="#333333")

    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="both", labelsize=7)    
    ax.xaxis.set_label_position("top")
    ax.xaxis.tick_top()    
    plt.tight_layout()
    highlight_unique_xticklabels(ax, shared_intervals)
    highlight_unique_yticklabels(ax, shared_intervals)
    out_path = os.path.join(out_dir, f"{prefix}_02_correlation_heatmap.png")
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close()

# Violin plots of per-SNP dz distributions for each interval, colored by up/down type
def plot_dz_distribution(dz_per_snp, interval_types, prefix, out_dir, shared_intervals):
    df_long  = dz_per_snp.reset_index().melt(id_vars="SNP", var_name="interval", value_name="dz")
    type_map = dict(zip(interval_types["interval"], interval_types["type"]))
    intervals = dz_per_snp.columns.tolist()
    palette  = {iv: (COLOR_G if type_map.get(iv) == "up" else COLOR_Y) for iv in intervals}
    fig, ax = plt.subplots(figsize=(14, 6))
    sns.violinplot(data=df_long, x="interval", y="dz", hue="interval", legend=False, order=intervals, 
                   palette=palette, inner="box", linewidth=1.2, alpha=0.75, ax=ax)
    ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
    ax.set_xlabel("Interval", fontsize=15)
    ax.set_ylabel("dz", fontsize=15)
    ax.tick_params(axis="both", labelsize=14)
    sns.despine(ax=ax)
    plt.tight_layout()
    highlight_unique_xticklabels(ax, shared_intervals)
    out_path = os.path.join(out_dir, f"{prefix}_03_dz_distribution_per_interval.png")
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close()

# Produces dz scatter plots between interval pairs with strongest pos. and neg. correlations (excluding adjacent ones)
def plot_dz_scatterplots(dz_per_snp, interval_types, prefix, out_dir, shared_intervals):
    intervals = dz_per_snp.columns.tolist()
    type_map  = dict(zip(interval_types["interval"], interval_types["type"]))
    n         = len(intervals)
    # Find strongest positive and negative correlation among non-adjacent pairs
    best_pos = (None, None, -np.inf)
    best_neg = (None, None,  np.inf)
    # Only consider non-adjacent pairs (i, j) where j >= i + 2 to avoid shared timepoint dependency
    for i in range(n):
        for j in range(i + 2, n):  # non-adjacent pairs only
            iv_i = intervals[i]
            iv_j = intervals[j]
            r, p = pearsonr(dz_per_snp[iv_i].values, dz_per_snp[iv_j].values)
            if r > best_pos[2]:
                best_pos = (iv_i, iv_j, r)
            if r < best_neg[2]:
                best_neg = (iv_i, iv_j, r)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    # Plot positive correlation pair in left panel, negative in right panel
    for ax, (iv_x, iv_y, r) in zip(axes, [best_pos, best_neg]):
        x = dz_per_snp[iv_x].values
        y = dz_per_snp[iv_y].values

        ax.scatter(x, y, alpha=0.4, s=3, color="black", rasterized=True)

        slope, intercept, r_val, p_val, _ = linregress(x, y)
        x_line = np.linspace(x.min(), x.max(), 100)
        # Draw regression line
        ax.plot(x_line, slope * x_line + intercept, color="red", linewidth=1.2)
        # Annotate with Pearson r value
        ax.text(0.05, 0.92, f"r = {r_val:.2f}",
                transform=ax.transAxes, fontsize=12, fontweight="bold",
                color="red" if r_val < 0 else "black")
        # Add dashed lines at x=0 and y=0
        ax.axhline(0, color="grey", linewidth=0.5, linestyle="--")
        ax.axvline(0, color="grey", linewidth=0.5, linestyle="--")
        # Color x and y labels by up/down type
        x_color = COLOR_G if type_map.get(iv_x) == "up" else COLOR_Y
        y_color = COLOR_G if type_map.get(iv_y) == "up" else COLOR_Y        
        xlabel = ax.set_xlabel(iv_x, fontsize=14, color=x_color, fontweight="bold")
        ylabel = ax.set_ylabel(iv_y, fontsize=14, color=y_color, fontweight="bold")

        # Highlight not in shared_intervals labels 
        if iv_x not in shared_intervals:
            xlabel.set_bbox(UNIQUE_BBOX)
        if iv_y not in shared_intervals:
            ylabel.set_bbox(UNIQUE_BBOX)

        ax.tick_params(labelsize=11)
        sns.despine(ax=ax)
    plt.tight_layout()
    out_path = os.path.join(out_dir, f"{prefix}_07_dz_scatterplots.png")
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close()

# Spaghetti plot of per-SNP dz trajectories across intervals, highlighting top and bottom SNPs by sync score
def plot_spaghetti(dz_per_snp, sync_scores, prefix, out_dir, shared_intervals):
    gs = sync_scores.set_index("SNP") if "SNP" in sync_scores.columns else sync_scores
    # Top N SNPs by sync score 
    top_snps_ids = gs.nsmallest(TOP_N, "rank").index.tolist()
    top_gs = gs.loc[top_snps_ids, "sync_score"].values
    gs_min, gs_max = top_gs.min(), top_gs.max()
    gs_norm = (top_gs - gs_min) / (gs_max - gs_min + 1e-9)
    # Orange gradient for top SNPs in phase with genome-wide trend
    cmap  = plt.cm.Oranges
    intervals = dz_per_snp.columns.tolist()
    x = np.arange(len(intervals))
    fig, ax = plt.subplots(figsize=(11, 6))

    for _, row in dz_per_snp.iterrows():
        ax.plot(x, row.values, color="gray",
                alpha=ALPHA_BG, linewidth=0.4, rasterized=True)

    for i, snp in enumerate(top_snps_ids):
        if snp not in dz_per_snp.index:
            continue
        color = cmap(0.35 + gs_norm[i] * 0.65)
        ax.plot(x, dz_per_snp.loc[snp].values, color=color, linewidth=1, alpha=0.9, rasterized=True)

    # Bottom N SNPs by sync score
    bottom_snps_ids = gs.nlargest(TOP_N, "rank").index.tolist()
    bottom_gs       = gs.loc[bottom_snps_ids, "sync_score"].values
    bs_min, bs_max  = bottom_gs.min(), bottom_gs.max()
    bs_norm         = (bottom_gs - bs_min) / (bs_max - bs_min + 1e-9)
    # Blue gradient for bottom SNPs out of phase with genome-wide trend
    cmap_blue = plt.cm.Blues

    for i, snp in enumerate(bottom_snps_ids):
        if snp not in dz_per_snp.index:
            continue
        color = cmap_blue(0.35 + bs_norm[i] * 0.65)
        ax.plot(x, dz_per_snp.loc[snp].values, color=color,
                linewidth=0.5, alpha=0.9, rasterized=True)
    
    ax.axhline(0, color="black", linewidth=0.5, linestyle="--", alpha=0.4)
    ax.set_xticks(x)
    ax.set_xticklabels(intervals, fontsize=10)
    ax.set_xlabel("Interval", fontsize=10)
    ax.set_ylabel("dz", fontsize=10)
    ax.tick_params(axis="y", labelsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    highlight_unique_xticklabels(ax, shared_intervals)
    out_path = os.path.join(out_dir, f"{prefix}_08_spaghetti.png")
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close()

# -----------------------------------------------------------------------
# HTML generation
# -----------------------------------------------------------------------
#

DESCRIPTIONS = {

    "01": (
        "- Mean change in z per temporal interval (dz) across all fluctuating SNPs "
        "(~13,000 after thinning to one SNP per gene at FDR &lt; 0.05). (Kelly figure 2C). <br>"
        "- Allele frequencies are polarized, so that the dz refers always to the minor allele (Kelly, Methods D).<br>"
        "- Orange bars (up intervals) indicate intervals where the mean minor allele frequency "
        "increased; blue bars (down intervals) where it decreased. <br>"
        "- The up/down classification is derived from the sign of the mean dz across all fluctuating SNPs in that interval. <br>"        
        "- * indicate that the mean is significantly different from zero by one-sample t-test (p-value &lt; 0.05). "
    ),

    "03": (
        "- Distribution of per-SNP dz values across all thinned fluctuating SNPs for each temporal interval. (not in paper)<br>"
        "- The white line shows the median and the thick bar the interquartile range.<br>"
    ),

    "02": (
        "- Pearson correlations between per-SNP dz values across all pairs of temporal intervals (Kelly figure 2B).<br> "
        "- Correlations between non-adjacent intervals, adjacent intervals are excluded (NA) because the shared timepoint parameter "
        "induces a statistical dependency between their estimates. (Kelly p. 509).<br>"   
        "- Same-type interval pairs (up-up or down-down) are expected to show positive correlations, and opposite-type pairs "
        "(up-down) negative correlations. Non-adjacent intervals should exhibit no correlation under drift.<br>" 
        "- ns cells indicate non-significant correlations (p &gt;= 0.05)." 
    ),

    "07": (
        "- Pairwise scatter plots of per-SNP dz values between an up-up and up-down intervals with the strongest correlations (Kelly figure 2A).<br>"
        "- Each point is one fluctuating SNP. The regression line and Pearson r are shown.<br>"        
        "- Under drift all pairwise correlations should be near zero.<br>" 
    ),

    "08": (
        "- dz trajectories across all temporal intervals for all thinned fluctuating SNPs (grey lines). <br>"
        "- The top 10 SNPs by sync score (Cg score in Kelly's) are highlighted in orange "
        " and the bottom 10 in blue (more out of phase). <br>"
        "- The sync score is computed per SNP as the sum of signed dz contributions: if a SNP dz "
        "in a time interval has the same sign as the genome-wide mean dz for that time interval, it contributes "
        "positively; if opposite, negatively  (Kelly p. 510). <br>"
        "- A high sync score indicates consistent in-phase oscillation with the genome-wide trend; a low score indicates out-of-phase oscillation. "
    ),
}

PLOT_ORDER = [
    ("01", "Mean dz per interval"),
    ("03", "dz distribution per interval"),
    ("02", "Correlation heatmap between intervals"),    
    ("07", "dz scatter plots between interval pairs"),    
    ("08", "dz trajectory plot - highlight top SNPs in phase and out of phase"),
]

PLOT_FILES = {
    "01": "_01_mean_dz_per_interval.png",
    "02": "_02_correlation_heatmap.png",
    "03": "_03_dz_distribution_per_interval.png",    
    "07": "_07_dz_scatterplots.png",
    "08": "_08_spaghetti.png",
}


def build_html(population, prefixes_in_pop, out_dir):
    rows_html = ""
    for plot_id, plot_title in PLOT_ORDER:
        desc  = DESCRIPTIONS[plot_id]
        cells = ""
        for prefix in prefixes_in_pop:
            img_file = f"{prefix}{PLOT_FILES[plot_id]}"
            cells += f"""
            <div>
                <h3>{prefix}</h3>
                <img src="{img_file}" alt="{prefix} {plot_title}">
            </div>"""
        rows_html += f"""
        <div class="plot-section">
            <h2>{plot_title}</h2>
            <p class="description">{desc}</p>
            <div class="plot-row">
                {cells}
            </div>
        </div>
        <hr>"""
    notice = "Highlighted labels indicate intervals not present in both E and L samples." if "SFIN" in population else ""
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <link rel="icon" type="image/png" href="../favicon_sbb.png">
    <title>Analysis of fluctuating SNPs for {population}</title>
    <style>
        body {{ font-family: -apple-system, Arial, sans-serif;
            font-size: 13px; color: #222; background: #fff; padding: 32px; }}
        h1 {{ font-size: 1.3rem; font-weight: 600; margin-bottom: 4px; }}
        h2 {{ font-size: 1rem; font-weight: 600; margin-bottom: 8px; }}
        h3 {{ font-size: 0.9rem; font-weight: 600; color: #555; margin-bottom: 8px; }}
        .page-subtitle {{ color: #333; font-size: 0.95rem; margin-bottom: 32px; line-height:1.3}}
        .plot-section {{ margin-bottom: 24px; }}
        .description {{ color: #555; font-size: 0.92rem; line-height: 1.6; 
            margin-bottom: 16px;max-width: 100%;}}
        .plot-row {{ display: grid; grid-template-columns: 1fr 1fr; gap: 16px;}}
        img {{ width: 80%; display: block; }}
        hr {{ border: none; border-top: 1px solid #eee; margin: 32px 0; }}
    </style>
</head>
<body>
    <h1>Analysis of fluctuating SNPs for {population} samples</h1>
    <p class="page-subtitle"> Allele frequency change at fluctuating SNPs (FDR &lt; 0.05) quantified as z (Fisher-Ford transformation)<br>
      <em>z = 2·arcsin(√p) <br> dz = z<sub>t</sub> − z<sub>t−1</sub></em><br>
      {notice}
    <hr>
    {rows_html}
</body>
</html>"""

    out_file = os.path.join(out_dir, f"{population}_fluctuating_SNPs.html")
    with open(out_file, "w") as f:
        f.write(html)
    print(f"  Saved HTML: {out_file}")

def generate_html(prefixes, out_dir):
    sfin = [p for p in prefixes if p.startswith("SFIN")]
    wfin = [p for p in prefixes if p.startswith("WFIN")]
    if sfin:
        build_html("SFIN", sfin, out_dir)
    if wfin:
        build_html("WFIN", wfin, out_dir)

def main():
    sns.set_style("whitegrid")
    os.makedirs(out_dir, exist_ok=True)
    prefixes = get_prefixes(work_dir)

    # Load all data
    all_data = {}
    for prefix in prefixes:
        all_data[prefix] = load_files(work_dir, prefix)
    # Compute shared intervals per pair (E <-> L)
    shared_map = {}
    for prefix in prefixes:
        if prefix not in all_data:
            continue
        pop, season = prefix.rsplit("_", 1)
        paired_season = "L" if season == "E" else "E"
        paired_prefix = f"{pop}_{paired_season}"
        if paired_prefix in all_data:
            my_intervals     = set(all_data[prefix][0]["interval"].tolist())
            paired_intervals = set(all_data[paired_prefix][0]["interval"].tolist())
            shared_map[prefix] = my_intervals & paired_intervals
            print(f"  {prefix} shared intervals with {paired_prefix}: {len(shared_map[prefix])}")
        else:
            # No pair found, all intervals treated as shared (no highlighting)
            shared_map[prefix] = set(all_data[prefix][0]["interval"].tolist())

    # Generate plots
    for prefix in prefixes:
        interval_types, ttests, correlations, dz_per_snp, sync_scores = all_data[prefix]
        shared_intervals = shared_map[prefix]

        plot_mean_dz_per_interval(ttests, interval_types, prefix, out_dir, shared_intervals)
        plot_dz_distribution(dz_per_snp, interval_types, prefix, out_dir, shared_intervals)
        plot_correlation_heatmap(correlations, dz_per_snp, prefix, out_dir, shared_intervals)
        plot_dz_scatterplots(dz_per_snp, interval_types, prefix, out_dir, shared_intervals)
        plot_spaghetti(dz_per_snp, sync_scores, prefix, out_dir, shared_intervals)

    generate_html(prefixes, out_dir)

if __name__ == "__main__":
    main()