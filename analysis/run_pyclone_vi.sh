#!/usr/bin/env bash
set -euo pipefail

# Batch-run PyClone-VI on the merged matched-pair/trio inputs used for the
# Figure 4 clone-dominance score plot.
#
# Inputs:
# - one CSV per matched case with columns:
#     mutation_id,sample_id,ref_counts,alt_counts,major_cn,minor_cn,normal_cn
#
# Outputs:
# - one `.h5` PyClone-VI fit object per input CSV
# - one `.results.tsv` table per input CSV
#
# This wrapper follows the official PyClone-VI CLI pattern:
#   pyclone-vi fit ...
#   pyclone-vi write-results-file ...
#
# Typical usage:
#   conda activate pyclone-vi
#   ./analysis/run_pyclone_vi.sh

INPUT_DIR="${INPUT_DIR:-/mnt/myvolume/panel_seq/new_bed_analysis/clone_analysis/pyclone/merged_trios}"
OUTPUT_DIR="${OUTPUT_DIR:-/mnt/myvolume/panel_seq/new_bed_analysis/clone_analysis/pyclone/pyclonevi_runs}"

NUM_CLUSTERS="${NUM_CLUSTERS:-40}"
DENSITY="${DENSITY:-beta-binomial}"
NUM_GRID_POINTS="${NUM_GRID_POINTS:-100}"
NUM_RESTARTS="${NUM_RESTARTS:-10}"

command -v pyclone-vi >/dev/null 2>&1 || {
  echo "ERROR: pyclone-vi not found in PATH. Activate the PyClone-VI environment first." >&2
  exit 1
}

mkdir -p "${OUTPUT_DIR}"
shopt -s nullglob
inputs=("${INPUT_DIR}"/*.csv)

if [[ ${#inputs[@]} -eq 0 ]]; then
  echo "No merged trio CSVs found in ${INPUT_DIR}." >&2
  exit 0
fi

echo "Input directory:  ${INPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Clusters:         ${NUM_CLUSTERS}"
echo "Density:          ${DENSITY}"
echo "Grid points:      ${NUM_GRID_POINTS}"
echo "Restarts:         ${NUM_RESTARTS}"

for in_csv in "${inputs[@]}"; do
  stem="$(basename "${in_csv}" .csv)"
  out_h5="${OUTPUT_DIR}/${stem}.h5"
  out_tsv="${OUTPUT_DIR}/${stem}.results.tsv"

  if [[ -s "${out_tsv}" ]]; then
    echo "Skipping ${stem}: ${out_tsv} already exists."
    continue
  fi

  echo "Running PyClone-VI for ${stem}"
  pyclone-vi fit \
    -i "${in_csv}" \
    -o "${out_h5}" \
    -c "${NUM_CLUSTERS}" \
    -d "${DENSITY}" \
    -g "${NUM_GRID_POINTS}" \
    -r "${NUM_RESTARTS}"

  pyclone-vi write-results-file \
    -i "${out_h5}" \
    -o "${out_tsv}"
done

echo "All PyClone-VI runs completed."
