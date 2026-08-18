#!/usr/bin/env bash
set -euo pipefail

BASE="/mnt/myvolume/panel_seq/new_bed_analysis"
RUN_BASE="/mnt/myvolume/panel_seq/v4_v5_clincnv/new_clincnv_runs/combined_runs/bed_split_s100_l3_f1"
RUN_ROOT="${RUN_BASE}/full_wes_ssSC_v5_s100_l3_f1_baf_bedcoverage_wes"
SCRIPT="/mnt/myvolume/panel_seq/reset_analysis_sample/test_2/ClinCNV/clinCNV.R"
THREADS="${CLINCNV_THREADS:-4}"
PAIR_FILE="${RUN_BASE}/pairs_ssSC_v5.csv"
SUBSET_BASE="${RUN_BASE}/coverage_subsets/ssSC_v5"

mkdir -p "${RUN_ROOT}"

CMD=(
  conda run -n mamba-env
  Rscript "${SCRIPT}"
  --hg38
  --normal "${SUBSET_BASE}/normal.ontarget.wes_mapq0_bedcoverage.cov"
  --tumor "${SUBSET_BASE}/tumor.ontarget.wes_mapq0_bedcoverage.cov"
  --normalOfftarget "${SUBSET_BASE}/normal.offtarget.wes_100kb_mapq10_bedcoverage.cov"
  --tumorOfftarget "${SUBSET_BASE}/tumor.offtarget.wes_100kb_mapq10_bedcoverage.cov"
  --bed "${BASE}/ssSC_v5.gc.genes.bed"
  --bedOfftarget "${BASE}/clincnv_offtarget_wes_100kb_filtered.bed"
  --pair "${PAIR_FILE}"
  --out "${RUN_ROOT}"
  --numberOfThreads "${THREADS}"
  --bafFolder "${BASE}/baf_from_pair_vcfs_clean"
  --colNum 4
  --scoreS 100
  --reanalyseCohort
  --lengthS 3
  --filterStep 1
  --clonePenalty 100
  --purityStep 2.5
  --clonalityForChecking 0.8
)

printf '%q ' "${CMD[@]}"
printf '\n'
"${CMD[@]}" | tee "${RUN_ROOT}/clinCNV.log"
