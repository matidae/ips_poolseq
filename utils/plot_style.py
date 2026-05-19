#!/usr/bin/env python3

#----------------------------------------------------------------------
# plot_style.py - unified color palette and style for all figures.
#
# Usage: 
#
#    from plot_style import apply_style, C, CMAP_DIV, style_boxplot
#    apply_style()   # call once at the top of every script
#
#    # Single-series default
#    ax.plot(x, y, color=C["teal"])
#
#    # Binary categorical (e.g. in-phase / anti-phase, up / down)
#    ax.bar(x, y, color=C["rust"])       # category A
#    ax.bar(x, y, color=C["steel"])      # category B
#
#    # Multi-category (up to 5)
#    for i, group in enumerate(groups):
#        ax.plot(x, y, color=C["cat"][i])
#
#    # Diverging colormap (e.g. heatmaps)
#    ax.imshow(data, cmap=CMAP_DIV)
#----------------------------------------------------------------------

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

#  Core palette 
#
#   teal   - single-series default (boxplots, scatter, line)
#   rust   - category A  (warm, replaces default matplotlib orange)
#   steel  - category B  (cool, replaces default matplotlib blue)
#   sage   - category C  (muted green)
#   plum   - category D  (muted purple)
#   sand   - category E  (warm neutral, 5th category or stacked bar fill)
#
C = {
    #  primary single-series
    "teal":         "#2E7D7D",
    "teal_light":   "#7DBFBF",
    "teal_dark":    "#1A4F4F",

    #  binary pair
    "rust":         "#C1440E",
    "rust_light":   "#F0B49A",
    "rust_dark":    "#7A2606",

    "steel":        "#2E5F8A",
    "steel_light":  "#A8C4DE",
    "steel_dark":   "#1A3A57",

    #  extended categorical
    "sage":         "#4A7C59",
    "sage_light":   "#AECFB8",

    "plum":         "#6B4C8A",
    "plum_light":   "#C4ADDA",

    "sand":         "#A8824A",
    "sand_light":   "#DEC99A",

    #  neutrals
    "grey_dark":    "#2C3E50",   # titles, tick labels
    "grey_mid":     "#64748B",   # axis labels, annotations
    "grey_light":   "#CBD5E1",   # grid lines, spines
    "grey_bg":      "#EEF2F5",   # axes background
    "white":        "#FFFFFF",

    #  semantic
    "median":       "#5A4A6F",   # median lines on boxplots (muted plum)
    "ref_line":     "#94A3B8",   # reference / threshold lines
    "band":         "#E4E4E4",   # alternating background bands

    # Pipeline depth sequence (raw -> qc -> mapped -> optical dedup -> full dedup)
    "pipe": ["#9E4A2E", "#D4724E", "#C4A882", "#4E9E9E", "#2E7D7D"],
}

# Ordered list for categorical cycling (up to 5 categories)
C["cat"] = [C["rust"], C["steel"], C["teal"], C["sage"], C["plum"]]

# Light fills matching the categorical order (for violin / bar / box fills)
C["cat_light"] = [
    C["rust_light"], C["steel_light"], C["teal_light"],
    C["sage_light"], C["plum_light"],
]

## Colormaps 
# Diverging: steel -> white -> rust  (replaces green→red)
CMAP_DIV = mcolors.LinearSegmentedColormap.from_list(
    "steel_rust",
    [C["steel_dark"], C["steel"], C["steel_light"],
     "#FFFFFF",
     C["rust_light"], C["rust"], C["rust_dark"]],
)

# Sequential teal (single-variable heatmaps / density)
CMAP_SEQ = mcolors.LinearSegmentedColormap.from_list(
    "seq_teal",
    [C["grey_bg"], C["teal_light"], C["teal"], C["teal_dark"]],
)

# rcParams 
# Call once per script, before any plotting. (font size can be customized)
def apply_style(font_size: int = 11):    
    plt.rcParams.update({
        # Font
        "font.family":        "DejaVu Sans",
        "font.size":          font_size,
        "axes.titlesize":     font_size + 4,
        "axes.titleweight":   "bold",
        "axes.titlecolor":    C["grey_dark"],
        "axes.labelsize":     font_size + 1,
        "axes.labelcolor":    C["grey_mid"],
        "xtick.labelsize":    font_size,
        "ytick.labelsize":    font_size,
        "legend.fontsize":    font_size - 1,

        # Colors
        "axes.facecolor":     C["grey_bg"],
        "figure.facecolor":   C["white"],
        "axes.edgecolor":     C["grey_light"],
        "axes.linewidth":     0.8,
        "xtick.color":        C["grey_mid"],
        "ytick.color":        C["grey_mid"],
        "text.color":         C["grey_dark"],

        # Grid
        "axes.grid":          True,
        "axes.grid.axis":     "y",
        "grid.color":         C["grey_light"],
        "grid.linewidth":     0.7,
        "axes.axisbelow":     True,

        # Spines
        "axes.spines.top":    False,
        "axes.spines.right":  False,
        "axes.spines.left":   False,
        "axes.spines.bottom": True,

        # Ticks
        "xtick.major.pad":    6,
        "ytick.major.pad":    6,

        # Legend
        "legend.frameon":     True,
        "legend.framealpha":  0.92,
        "legend.edgecolor":   C["grey_light"],

        # Lines / markers
        "lines.linewidth":    1.6,
        "lines.markersize":   5,

        # Figure
        "figure.dpi":         100,
        "savefig.dpi":        100,
        "savefig.bbox":       "tight",
        "savefig.facecolor":  C["white"],

        # Default color cycle
        "axes.prop_cycle": plt.cycler(color=C["cat"]),
    })


#  Boxplot helper 
def style_boxplot(bp, face=None, edge=None, median=None, whisker=None):
    face   = face   or C["teal_light"]
    edge   = edge   or C["teal_dark"]
    median = median or C["median"]
    whisker = whisker or C["teal"]

    for patch in bp["boxes"]:
        patch.set(facecolor=face, edgecolor=edge, linewidth=1.3, alpha=0.9)
    for line in bp["medians"]:
        line.set(color=median, linewidth=2.2)
    for line in bp["whiskers"]:
        line.set(color=whisker, linewidth=1.0, linestyle="--")
    for cap in bp["caps"]:
        cap.set(color=whisker, linewidth=1.3)


#  Alternating year-band overlay 
def add_year_bands(ax, years, reps_per_year=2, color=None, alpha=0.55):
    color = color or C["band"]
    for i, _ in enumerate(years):
        if i % 2 == 0:
            x0 = i * reps_per_year + 0.5
            ax.axvspan(x0, x0 + reps_per_year, color=color, alpha=alpha, zorder=0)

#  Reference h-lines 
def add_hlines(ax, values, color=None, lw=0.9, ls="--"):
   
    color = color or C["ref_line"]
    for v in values:
        ax.axhline(v, color=color, linewidth=lw, linestyle=ls, zorder=2)