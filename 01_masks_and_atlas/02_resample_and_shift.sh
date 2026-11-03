#!/bin/bash
# ===============================================
# 2 mm seed masks used for time-series extraction
# ===============================================

set -euo pipefail

native_dir=${1:-masks/native}
out_dir=${2:-masks/2mm}
ventricle_mask=${3:-masks/brain_excluding_fourth_ventricle_2mm.nii.gz}   
mkdir -p "${out_dir}"

# ---- LC ---------------------------------------------------------------
ResampleImage 3 "${native_dir}/LC_dilated_0.5mm.nii.gz" "${out_dir}/LC_dilated_2mm_unmasked.nii.gz" 2x2x2 0 1
fslmaths "${out_dir}/LC_dilated_2mm_unmasked.nii.gz" -mas "${ventricle_mask}" "${out_dir}/LC_2mm.nii.gz"
echo "LC: $(fslstats "${out_dir}/LC_dilated_2mm_unmasked.nii.gz" -V | awk '{print $1}') voxels ->" \
     "$(fslstats "${out_dir}/LC_2mm.nii.gz" -V | awk '{print $1}') after removing fourth-ventricle voxels"

# ---- Shifted LC ---------------------------------------------------------
shift_matrix=$(mktemp)
printf "1 0 0 0\n0 1 0 12\n0 0 1 -9\n0 0 0 1\n" > "${shift_matrix}"   # +12 mm y (anterior), -9 mm z (caudal)
flirt -in "${out_dir}/LC_2mm.nii.gz" -ref "${out_dir}/LC_2mm.nii.gz" -applyxfm -init "${shift_matrix}" \
      -interp nearestneighbour -out "${out_dir}/shifted_LC_2mm.nii.gz"
rm -f "${shift_matrix}"
echo "Shifted LC: $(fslstats "${out_dir}/shifted_LC_2mm.nii.gz" -V | awk '{print $1}') voxels"

# ---- VTA ------------------------------------------------------------------
if [ -f "${native_dir}/VTA_bilateral_native.nii.gz" ]; then
  ResampleImage 3 "${native_dir}/VTA_bilateral_native.nii.gz" "${out_dir}/VTA_2mm.nii.gz" 2x2x2 0 1
  echo "VTA: $(fslstats "${out_dir}/VTA_2mm.nii.gz" -V | awk '{print $1}') voxels"
fi
