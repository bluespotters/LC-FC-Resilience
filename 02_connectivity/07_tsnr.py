"""
tSNR checks
"""

import argparse
import os

import numpy as np
import pandas as pd
from nilearn.image import load_img
from nilearn.maskers import NiftiMasker

FMRI_FILE = os.path.join("func", "denoised_func_data_nonaggr_wm_csf_tfilt.nii.gz")


def voxelwise_tsnr(fmri_img, mask_img):
    ts = NiftiMasker(mask_img=mask_img, standardize=False).fit_transform(fmri_img)   # time x voxels
    return np.nanmean(ts.mean(axis=0) / ts.std(axis=0))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--subjects-dir", required=True)
    ap.add_argument("--subjects", required=True)
    ap.add_argument("--mask", action="append", required=True, help="NAME=path; repeat for each mask")
    ap.add_argument("--out-csv", required=True)
    a = ap.parse_args()

    masks = {name: load_img(path) for name, path in (m.split("=", 1) for m in a.mask)}
    rows = []
    for sub in (s.strip() for s in open(a.subjects) if s.strip()):
        fmri = load_img(os.path.join(a.subjects_dir, sub, FMRI_FILE))
        rows.append({"sub": sub, **{name: voxelwise_tsnr(fmri, mask) for name, mask in masks.items()}})
    pd.DataFrame(rows).sort_values("sub").to_csv(a.out_csv, index=False)
    print(f"tSNR for {len(rows)} participants written to {a.out_csv}")


if __name__ == "__main__":
    main()
