"""
Combine per-hemisphere parcel thickness tables and harmonise 

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
REQUIRED = ["age", "sex", "YrsEd", "APOE_status", "centiloid", "DSST", "LogicMem", "MMSE", "FCSR", BATCH]
MIN_BATCH_SIZE = 15
STATS_COLUMNS = ["StructName", "NumVert", "SurfArea", "GrayVol", "ThickAvg",
                 "ThickStd", "MeanCurv", "GausCurv", "FoldInd", "CurvInd"]


def read_thickness(path):
    """ThickAvg column of an mris_anatomical_stats table (lines starting with # skipped)."""
    if not os.path.exists(path):
        return None
    return pd.read_csv(path, comment="#", sep=r"\s+", names=STATS_COLUMNS)["ThickAvg"].values


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--covariates", required=True)
    ap.add_argument("--subjects", required=True)
    ap.add_argument("--tables-dir", required=True)
    ap.add_argument("--parcellation", type=int, default=200)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    subjects = [s.strip() for s in open(a.subjects) if s.strip()]
    cov = pd.read_csv(a.covariates)
    cov = cov[cov["sub"].isin(subjects)]
    cov = cov[cov[REQUIRED].notna().all(axis=1)].copy()

    # Left then right hemisphere thickness of each participant with both tables
    per_hemi = a.parcellation // 2
    table_file = lambda sub, hemi: os.path.join(
        a.tables_dir, f"{sub}_{hemi}.Schaefer2018_{a.parcellation}Parcels_7Networks.table")
    rows, keep = [], []
    for sub in cov["sub"]:
        hemis = [read_thickness(table_file(sub, h)) for h in ("lh", "rh")]
        if all(h is not None and len(h) == per_hemi for h in hemis):
            rows.append(np.concatenate(hemis))
            keep.append(sub)
    cov = cov[cov["sub"].isin(keep)]
    # Scanner batches with >= 15 participants & complete covariates
    batch_sizes = cov[BATCH].value_counts()
    in_batch = cov[BATCH].isin(batch_sizes[batch_sizes >= MIN_BATCH_SIZE].index).values
    cov, keep = cov[in_batch], cov["sub"][in_batch].tolist()
    data = np.array(rows)[in_batch].T                              # parcels x participants

    covars = cov[REQUIRED].copy()
    for col in CATEGORICAL + [BATCH]:
        covars[col] = LabelEncoder().fit_transform(covars[col])

    harmonised = neuroCombat(dat=data, covars=covars, batch_col=BATCH,
                             categorical_cols=CATEGORICAL, continuous_cols=CONTINUOUS)["data"]

    out = pd.DataFrame(harmonised.T, columns=[f"ROI_{i + 1}" for i in range(a.parcellation)])
    out.insert(0, "sub", keep)
    out.to_csv(os.path.join(a.out_dir, f"cortical_thickness_{a.parcellation}.csv"), index=False)
    print(f"Harmonised cortical thickness: {len(out)} participants x {a.parcellation} parcels "
          f"across {cov[BATCH].nunique()} scanner batches")


if __name__ == "__main__":
    main()
