"""
NeuroCombat harmonisation.
"""

import argparse
import os

import numpy as np
import pandas as pd
from neuroCombat import neuroCombat
from sklearn.preprocessing import LabelEncoder

CATEGORICAL = ["sex", "APOE_status"]
CONTINUOUS = ["age", "YrsEd", "centiloid", "DSST", "LogicMem", "MMSE", "FCSR"]
BATCH = "DeviceSerialNumber"
REQUIRED = ["age", "sex", "YrsEd", "APOE_status", "centiloid", "DSST", "LogicMem", "MMSE", "FCSR",
            BATCH]
MIN_BATCH_SIZE = 15


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--covariates", required=True)
    ap.add_argument("--subjects", required=True)
    ap.add_argument("--matrix-dir", required=True)
    ap.add_argument("--parcellation", required=True)
    ap.add_argument("--seed", required=True)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    subjects = [s.strip() for s in open(a.subjects) if s.strip()]
    cov = pd.read_csv(a.covariates)
    cov = cov[cov["sub"].isin(subjects)]
    cov = cov[cov[REQUIRED].notna().all(axis=1)].copy()
    matrix_file = lambda sub: os.path.join(a.matrix_dir, f"{sub}_{a.parcellation}_{a.seed}_matrix.csv")
    cov = cov[[os.path.exists(matrix_file(s)) for s in cov["sub"]]]
    # Scanner batches with >= 15 participants
    batch_sizes = cov[BATCH].value_counts()
    cov = cov[cov[BATCH].isin(batch_sizes[batch_sizes >= MIN_BATCH_SIZE].index)]

    # Lower triangle of each matrix
    vectors, n = [], None
    for sub in cov["sub"]:
        m = pd.read_csv(matrix_file(sub), header=None).values
        n = n or m.shape[0]
        if m.shape != (n, n):
            raise ValueError(f"Inconsistent matrix size for {sub}")
        vectors.append(m[np.tril_indices(n, k=-1)])
    data = np.array(vectors).T                                     # edges x participants

    covars = cov[REQUIRED].copy()
    for col in CATEGORICAL + [BATCH]:
        covars[col] = LabelEncoder().fit_transform(covars[col])

    harmonised = neuroCombat(dat=data, covars=covars, batch_col=BATCH,
                             categorical_cols=CATEGORICAL, continuous_cols=CONTINUOUS)["data"]

    for i, sub in enumerate(cov["sub"]):
        m = np.zeros((n, n))
        m[np.tril_indices(n, k=-1)] = harmonised[:, i]
        np.savetxt(os.path.join(a.out_dir, f"{sub}_{a.parcellation}_{a.seed}_matrix.csv"), m + m.T, delimiter=",")
    print(f"Harmonised {len(cov)} participants across {cov[BATCH].nunique()} scanner batches "
          f"({n} x {n} matrices)")


if __name__ == "__main__":
    main()
