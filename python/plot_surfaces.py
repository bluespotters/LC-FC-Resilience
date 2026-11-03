"""
Cortical (surfplot) and subcortical (subcortex_visualization) maps of parcel values.

Parcels are identified by ROI_Number (1-32 = Tian S2 subcortex, 33+ = Schaefer cortex,
in atlas order).
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from brainspace.datasets import load_parcellation
from neuromaps.datasets import fetch_fslr
from subcortex_visualization.plotting import plot_subcortical_data
from surfplot import Plot

N_SUBCORTICAL = 32
LABEL_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Tian_Subcortex_S2_3T_label.txt")

# Tian S2 label -> subcortex_visualization (Melbourn S2) region name
S2_TO_REGION = {
    "aHIP": "hippocampus_anterior", "pHIP": "hippocampus_posterior",
    "lAMY": "amygdala_lateral", "mAMY": "amygdala_medial",
    "THA-DP": "thalamus_DP", "THA-VP": "thalamus_VP", "THA-VA": "thalamus_VA", "THA-DA": "thalamus_DA",
    "NAc-shell": "accumbens_shell", "NAc-core": "accumbens_core",
    "pGP": "pallidum_posterior", "aGP": "pallidum_anterior",
    "aPUT": "putamen_anterior", "pPUT": "putamen_posterior",
    "aCAU": "caudate_anterior", "pCAU": "caudate_posterior",
}


def make_colormap(name, custom_colors=None):
    if name == "teal_tomato":                # diverging effect maps
        return mcolors.LinearSegmentedColormap.from_list("teal_tomato", ["#008080", "#FFFFFF", "#FF6347"])
    if name == "significance_flow":          # 1 = remained, 2 = gained, 3 = lost
        return mcolors.ListedColormap(["#696969", "#008080", "#FF6347"])
    if name == "custom":
        return mcolors.LinearSegmentedColormap.from_list("custom", [c.strip() for c in custom_colors.split(",")])
    return plt.get_cmap(name)


def plot_subcortex(results, value_col, fdr_col, alpha, cmap, vmin, vmax, out_dir,
                   mask_nonsignificant=True, transparent=False):
    with open(LABEL_FILE) as f:
        labels = [line.strip() for line in f]

    rows = []
    for _, r in results.iterrows():
        label = labels[int(r["ROI_Number"]) - 1]                 # e.g. "THA-VA-rh"
        value = r[value_col]
        if mask_nonsignificant and fdr_col:
            value = value if r[fdr_col] < alpha else 0           # show pFDR < alpha only
        if pd.isna(value):                                       # leave undrawn (flow maps)
            continue
        rows.append({"region": S2_TO_REGION[label.rsplit("-", 1)[0]], "value": value,
                     "Hemisphere": "L" if label.endswith("-lh") else "R"})
    if not rows:
        print("No subcortical regions to plot")
        return

    fig = plt.figure(figsize=(12, 6))
    plot_subcortical_data(
        subcortex_data=pd.DataFrame(rows), atlas="Melbourne_S2", hemisphere="both",
        line_thickness=1.5, line_color="black", fill_title=f"Subcortical {value_col}",
        cmap=cmap, vmin=vmin, vmax=vmax,
        # centre the scale at 0 for diverging maps, otherwise at the middle of the range
        # (keeps the colour scale linear for e.g. 0-0.1 mean maps or 1-3 flow maps)
        midpoint=0 if vmin < 0 < vmax else (vmin + vmax) / 2.0,
    )
    path = os.path.join(out_dir, f"subcortical_{value_col}_plot.png")
    plt.tight_layout()
    plt.savefig(path, dpi=400, bbox_inches="tight", transparent=transparent)
    plt.close("all")
    print("Saved", path)


def plot_cortex(results, value_col, fdr_col, alpha, cmap, vmin, vmax, out_dir,
                parcellation, add_outline=True, layout="grid"):
    results = results.sort_values("ROI_Number").reset_index(drop=True)
    values = results[value_col].values
    significant = results[fdr_col].values < alpha if (add_outline and fdr_col) else None

    lh, rh = fetch_fslr()["inflated"]
    atlas = load_parcellation("schaefer", scale=parcellation, join=True)
    data = np.full(atlas.shape, np.nan)
    outline = np.full(atlas.shape, np.nan)
    for i, label in enumerate(l for l in np.unique(atlas) if l != 0):   # parcels in atlas order
        if i < len(values):
            data[atlas == label] = values[i]
            if significant is not None and significant[i]:
                outline[atlas == label] = 1

    p = Plot(lh, rh, views=["lateral", "medial"], layout=layout,
             zoom=1.25 if layout == "row" else 1.5, brightness=1)
    p.add_layer(data, cmap=cmap, cbar=True, color_range=(vmin, vmax))
    if significant is not None:
        p.add_layer(outline, cmap="gray", as_outline=True, cbar=False)
    path = os.path.join(out_dir, f"surface_{value_col}_plot.png")
    p.build().savefig(path, dpi=700, bbox_inches="tight")
    print("Saved", path)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results_csv")
    ap.add_argument("value_column")
    ap.add_argument("fdr_column", nargs="?", default=None,
                    help="FDR-corrected p-values: outline (cortex) / show only (subcortex) significant parcels")
    ap.add_argument("--parcellation", type=int, choices=[200, 400], default=200)
    ap.add_argument("--colormap", default="teal_tomato",
                    help="teal_tomato, significance_flow, custom or any matplotlib colormap")
    ap.add_argument("--custom-colors", help='comma-separated hex colours for --colormap custom, e.g. "#FFFFFF,#FF6347"')
    ap.add_argument("--vmin", type=float, required=True)
    ap.add_argument("--vmax", type=float, required=True)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--layout", choices=["grid", "column", "row"], default="grid")
    ap.add_argument("--transparent-subcortical", action="store_true")
    a = ap.parse_args()

    results = pd.read_csv(a.results_csv)
    fdr_col = a.fdr_column if a.fdr_column in results.columns else None
    cmap = make_colormap(a.colormap, a.custom_colors)
    out_dir = os.path.dirname(os.path.abspath(a.results_csv))

    sub = results[results["ROI_Number"] <= N_SUBCORTICAL]
    ctx = results[results["ROI_Number"] > N_SUBCORTICAL]
    if not sub.empty:
        plot_subcortex(sub, a.value_column, fdr_col, a.alpha, cmap, a.vmin, a.vmax, out_dir,
                       mask_nonsignificant=fdr_col is not None, transparent=a.transparent_subcortical)
    if not ctx.empty:
        plot_cortex(ctx, a.value_column, fdr_col, a.alpha, cmap, a.vmin, a.vmax, out_dir,
                    a.parcellation, add_outline=fdr_col is not None, layout=a.layout)


if __name__ == "__main__":
    main()
