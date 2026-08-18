#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/myvolume/panel_seq/v4_v5_clincnv/new_clincnv_runs/combined_runs/bed_split_s100_l3_f1"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
V4_SCREEN="clincnv_ssSC_v4LI_s100_l3_f1"
V5_SCREEN="clincnv_ssSC_v5_s100_l3_f1"

python "${SCRIPT_DIR}/build_pair_subsets.py"

if [[ ! -e "${ROOT}/coverage_subsets" ]]; then
  python "${SCRIPT_DIR}/build_subset_coverages.py"
fi

screen -dmS "${V4_SCREEN}" bash -lc "export CLINCNV_THREADS=4 && bash '${SCRIPT_DIR}/run_clincnv_ssSC_v4LI_s100_l3_f1_baf_bedcoverage_wes.sh'"
screen -dmS "${V5_SCREEN}" bash -lc "export CLINCNV_THREADS=4 && bash '${SCRIPT_DIR}/run_clincnv_ssSC_v5_s100_l3_f1_baf_bedcoverage_wes.sh'"

echo "Launched screens:"
echo "  ${V4_SCREEN}"
echo "  ${V5_SCREEN}"
