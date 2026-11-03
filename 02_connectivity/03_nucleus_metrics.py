"""
Generate FC metrics.
"""

import argparse
import glob
import os

import numpy as np
import pandas as pd


def network_means(fc, assignments):
    """Mean seed connectivity to each network (parcels ordered by ROI_Number)."""
    networks = assignments.sort_values("ROI_Number")["Network"].values
    return {f"{net}_mean": fc[:, networks == net].mean(axis=1) for net in sorted(set(networks))}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix-dir", required=True)
    ap.add_argument("--parcellation", required=True)
    ap.add_argument("--seed", required=True)
    ap.add_argument("--assignments-dir", required=True,
                    help="folder with roi_network_assignments_yeo17_<parc>.csv and _yeo7_<parc>.csv")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    suffix = f"_{a.parcellation}_{a.seed}_matrix.csv"
    files = sorted(glob.glob(os.path.join(a.matrix_dir, f"*{suffix}")))
    subs = [os.path.basename(f)[: -len(suffix)] for f in files]
    fc = np.array([pd.read_csv(f, header=None).values[-1, :-1] for f in files])   # seed row - no self
    n_parcels = fc.shape[1]

    stem = os.path.join(a.out_dir, f"{a.parcellation}_{a.seed}")
    roi = pd.DataFrame(fc, columns=[f"ROI_{i + 1}" for i in range(n_parcels)])
    roi.insert(0, "sub", subs)
    roi.to_csv(f"{stem}_ROI_correlations.csv", index=False)

    for partition, out_suffix in [("yeo17", ""), ("yeo7", "_yeo7")]:
        assignments = pd.read_csv(os.path.join(a.assignments_dir, f"roi_network_assignments_{partition}_{a.parcellation}.csv"))
        if len(assignments) != n_parcels:
            raise ValueError(f"{partition}: {len(assignments)} assignments for {n_parcels} parcels")
        table = pd.DataFrame({"sub": subs, "global_mean": fc.mean(axis=1), **network_means(fc, assignments)})
        table.to_csv(f"{stem}_network_connectivity{out_suffix}.csv", index=False)
    print(f"{len(subs)} participants, {n_parcels} parcels -> {stem}_*.csv")


if __name__ == "__main__":
    main()
