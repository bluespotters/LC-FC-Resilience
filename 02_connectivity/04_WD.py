"""
Weighted degree of each cortical parcel, used as a tissue control.
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd

N_SUBCORTICAL = 32


def weighted_degree(matrix):
    m = np.abs(matrix)
    np.fill_diagonal(m, 0)
    return m.sum(axis=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix-dir", required=True)
    ap.add_argument("--parcellation", required=True)
    ap.add_argument("--seed", default="LC")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()

    suffix = f"_{a.parcellation}_{a.seed}_matrix.csv"
    rows = []
    for f in sorted(glob.glob(os.path.join(a.matrix_dir, f"*{suffix}"))):
        m = pd.read_csv(f, header=None).values[:-1, :-1]                 # drop seed
        cortical = m[N_SUBCORTICAL:, N_SUBCORTICAL:]
        rows.append([os.path.basename(f)[: -len(suffix)], *weighted_degree(cortical)])
    n = len(rows[0]) - 1
    out = pd.DataFrame(rows, columns=["sub"] + [f"ROI_{i + 1}" for i in range(n)])
    os.makedirs(a.out_dir, exist_ok=True)
    out.to_csv(os.path.join(a.out_dir, f"weighted_degree_{a.parcellation}.csv"), index=False)
    print(f"Weighted degree for {len(out)} participants x {n} cortical parcels")


if __name__ == "__main__":
    main()
