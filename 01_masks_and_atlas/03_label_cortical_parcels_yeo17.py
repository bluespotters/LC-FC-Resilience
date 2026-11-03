"""
Sort Yeo7 & Yeo-17 labels together for the Schaefer parcels.
"""

import argparse
import os

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--centroids-dir", required=True)
    ap.add_argument("--parcellation", type=int, default=200)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()

    name = "Schaefer2018_{}Parcels_{}Networks_order_FSLMNI152_2mm.Centroid_RAS.csv"
    yeo7 = pd.read_csv(os.path.join(a.centroids_dir, name.format(a.parcellation, 7)))
    yeo17 = pd.read_csv(os.path.join(a.centroids_dir, name.format(a.parcellation, 17)))
    c17 = yeo17[["R", "A", "S"]].values

    rows = []
    for i, r in yeo7.iterrows():
        j = int(np.argmin(np.sqrt(((c17 - r[["R", "A", "S"]].values.astype(float)) ** 2).sum(axis=1))))
        rows.append({"yeo7_index": i, "yeo7_name": r["ROI Name"], "yeo7_network": r["ROI Name"].split("_")[2],
                     "yeo17_index": j, "yeo17_name": yeo17.loc[j, "ROI Name"],
                     "yeo17_network": yeo17.loc[j, "ROI Name"].split("_")[2]})
    os.makedirs(a.out_dir, exist_ok=True)
    out = os.path.join(a.out_dir, f"cortical_yeo17_yeo7_labels_{a.parcellation}.csv")
    pd.DataFrame(rows).to_csv(out, index=False)
    print("Saved", out)


if __name__ == "__main__":
    main()

