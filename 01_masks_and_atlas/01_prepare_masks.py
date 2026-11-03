"""
Prepare seed masks.
"""

import argparse
import os

import nibabel as nib
import numpy as np
from scipy.ndimage import binary_dilation


def dilate_lc(keren_path, out_path):
    img = nib.load(keren_path)
    dilated = binary_dilation(img.get_fdata() > 0, structure=np.ones((2, 2, 3), dtype=bool))
    nib.save(nib.Nifti1Image(dilated.astype(np.uint8), img.affine, img.header), out_path)
    print(f"LC: {int((img.get_fdata() > 0).sum())} -> {int(dilated.sum())} voxels after dilation: {out_path}")


def bilateral_vta(left_path, right_path, out_path):
    left, right = nib.load(left_path), nib.load(right_path)
    if not np.allclose(left.affine, right.affine):
        raise ValueError("Left and right VTA masks are not on the same grid")
    mask = (left.get_fdata() > 0) | (right.get_fdata() > 0)
    nib.save(nib.Nifti1Image(mask.astype(np.uint8), left.affine, left.header), out_path)
    print(f"VTA: {int(mask.sum())} voxels (bilateral): {out_path}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--keren", required=True, help="Keren 2SD LC template")
    ap.add_argument("--vta-left", help="Brainstem Navigator VTA-PBP, L")
    ap.add_argument("--vta-right", help="Brainstem Navigator VTA-PBP, R")
    ap.add_argument("--out-dir", default="masks/native")
    a = ap.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)
    dilate_lc(a.keren, os.path.join(a.out_dir, "LC_dilated_0.5mm.nii.gz"))
    if a.vta_left and a.vta_right:
        bilateral_vta(a.vta_left, a.vta_right, os.path.join(a.out_dir, "VTA_bilateral_native.nii.gz"))
