"""
Parcel and seed time series -> connectivity matrices.
"""

import argparse
import os

import numpy as np
from nilearn.image import load_img
from nilearn.maskers import NiftiLabelsMasker, NiftiMasker

FMRI_FILE = os.path.join("func", "denoised_func_data_nonaggr_wm_csf_tfilt.nii.gz")


def seed_timeseries(fmri_img, seed_mask):
    """Mean signal across the seed mask voxels (n_timepoints x 1)."""
    voxels = NiftiMasker(mask_img=seed_mask, standardize=False).fit_transform(fmri_img)
    return voxels.mean(axis=1, keepdims=True)


def connectivity_matrix(timeseries):
    """Fisher z-transformed Pearson correlation matrix with a zero diagonal."""
    r = np.corrcoef(timeseries, rowvar=False)
    np.fill_diagonal(r, 0)
    return np.arctanh(r)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--subjects-dir", required=True)
    ap.add_argument("--subjects", required=True, help="text file with one participant ID per line")
    ap.add_argument("--atlas", required=True, help="Schaefer-Tian labels")
    ap.add_argument("--parcellation", required=True, help="e.g. 200")
    ap.add_argument("--seed-mask", required=True)
    ap.add_argument("--seed", required=True, help="seed in file names: LC, shifted_LC or VTA")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)
    subjects = [s.strip() for s in open(a.subjects) if s.strip()]
    parcel_masker = NiftiLabelsMasker(labels_img=load_img(a.atlas), standardize=False)
    seed_mask = load_img(a.seed_mask)

    for i, sub in enumerate(subjects, 1):
        fmri_path = os.path.join(a.subjects_dir, sub, FMRI_FILE)
        fmri = load_img(fmri_path)
        timeseries = np.column_stack([parcel_masker.fit_transform(fmri), seed_timeseries(fmri, seed_mask)])
        stem = os.path.join(a.out_dir, f"{sub}_{a.parcellation}_{a.seed}")
        np.savetxt(f"{stem}_timeseries.csv", timeseries, delimiter=",")
        np.savetxt(f"{stem}_matrix.csv", connectivity_matrix(timeseries), delimiter=",")
        print(f"[{i}/{len(subjects)}] {sub}: {timeseries.shape[1]} x {timeseries.shape[1]} matrix")


if __name__ == "__main__":
    main()
