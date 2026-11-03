#!/bin/bash
# ==========================================
# Cortical thickness per Schaefer parcel
# ==========================================

set -euo pipefail

subjects_dir=$1; subject_list=$2; gcs_dir=$3; out_dir=$4; parcellation=${5:-200}
mkdir -p "${out_dir}"

while read -r sub; do
  [ -z "${sub}" ] && continue
  export SUBJECTS_DIR="${subjects_dir}/${sub}"
  for hemi in lh rh; do
    annot="${SUBJECTS_DIR}/freesurfer/label/${hemi}.Schaefer2018_${parcellation}Parcels_7Networks_order.annot"
    mris_ca_label -l "${SUBJECTS_DIR}/freesurfer/label/${hemi}.cortex.label" freesurfer "${hemi}" \
      "${SUBJECTS_DIR}/freesurfer/surf/${hemi}.sphere.reg" \
      "${gcs_dir}/${hemi}.Schaefer2018_${parcellation}Parcels_7Networks.gcs" "${annot}"
    mris_anatomical_stats -f "${out_dir}/${sub}_${hemi}.Schaefer2018_${parcellation}Parcels_7Networks.table" \
      -b -a "${annot}" freesurfer "${hemi}"
  done
  echo "${sub} done"
done < "${subject_list}"
