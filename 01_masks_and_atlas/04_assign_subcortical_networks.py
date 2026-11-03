"""
Assign Tian S2 subcortical parcels to a network.
"""

import argparse
import os

import nibabel as nib
import numpy as np
import pandas as pd

FMRI_FILE = os.path.join("func", "denoised_func_data_nonaggr_wm_csf_tfilt.nii.gz")
N_SUB = 32
YEO17_TO_YEO7 = {"VisCent": "Vis", "VisPeri": "Vis", "SomMotA": "SomMot", "SomMotB": "SomMot",
                 "DorsAttnA": "DorsAttn", "DorsAttnB": "DorsAttn", "SalVentAttnA": "SalVentAttn",
                 "SalVentAttnB": "SalVentAttn", "LimbicA": "Limbic", "LimbicB": "Limbic", "ContA": "Cont",
                 "ContB": "Cont", "ContC": "Cont", "DefaultA": "Default", "DefaultB": "Default",
                 "DefaultC": "Default", "TempPar": "Default"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--subjects-dir", required=True)
    ap.add_argument("--subjects", required=True)
    ap.add_argument("--atlas", required=True)
    ap.add_argument("--tian-labels", required=True)
    ap.add_argument("--cortical-labels", required=True)
    ap.add_argument("--parcellation", type=int, default=200)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()

    atlas = nib.load(a.atlas).get_fdata()
    n_total = N_SUB + a.parcellation
    cortex = pd.read_csv(a.cortical_labels).sort_values("yeo7_index")
    tian = [l.strip() for l in open(a.tian_labels)]

    # Subcortical x cortical correlation, averaged over participants
    corr_sum, n = np.zeros((N_SUB, a.parcellation)), 0
    for sub in (s.strip() for s in open(a.subjects) if s.strip()):
        path = os.path.join(a.subjects_dir, sub, FMRI_FILE)
        if not os.path.exists(path):
            continue
        data = nib.load(path).get_fdata()
        ts = np.array([data[atlas == label].mean(axis=0) for label in range(1, n_total + 1)])
        corr_sum += np.corrcoef(ts)[:N_SUB, N_SUB:]
        n += 1
    mean_corr = corr_sum / n
    print(f"Averaged subcortical-cortical correlations over {n} participants")

    yeo17_net = cortex["yeo17_network"].values
    networks = list(YEO17_TO_YEO7)
    rows17, rows7 = [], []
    for i in range(N_SUB):
        by_net = {net: mean_corr[i, yeo17_net == net].mean() for net in networks if (yeo17_net == net).any()}
        best = max(by_net, key=by_net.get)
        rows17.append({"ROI_Number": i + 1, "ROI_Name": tian[i], "Network": best})
        rows7.append({"ROI_Number": i + 1, "ROI_Name": tian[i], "Network": YEO17_TO_YEO7[best]})
    for k, r in enumerate(cortex.itertuples()):
        rows17.append({"ROI_Number": N_SUB + k + 1, "ROI_Name": r.yeo17_name, "Network": r.yeo17_network})
        rows7.append({"ROI_Number": N_SUB + k + 1, "ROI_Name": r.yeo7_name, "Network": r.yeo7_network})

    os.makedirs(a.out_dir, exist_ok=True)
    for rows, tag in [(rows17, "yeo17"), (rows7, "yeo7")]:
        pd.DataFrame(rows).to_csv(os.path.join(a.out_dir, f"roi_network_assignments_{tag}_{a.parcellation}.csv"), index=False)
    print("Subcortical assignments:", {r["ROI_Name"]: r["Network"] for r in rows17[:N_SUB]})


if __name__ == "__main__":
    main()
