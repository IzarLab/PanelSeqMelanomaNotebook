#!/usr/bin/env bash
set -euo pipefail

BASE="/mnt/myvolume/panel_seq/new_bed_analysis"
RUN_ROOT="${BASE}/new_clincnv_runs/full_wes_s200_l9_f2_no_baf_bedcoverage_wes"
SCRIPT="/mnt/myvolume/panel_seq/reset_analysis_sample/test_2/ClinCNV/clinCNV.R"
THREADS="${CLINCNV_THREADS:-8}"

mkdir -p "${RUN_ROOT}"

CMD=(
  conda run -n mamba-env
  Rscript "${SCRIPT}"
  --hg38
  --normal "${BASE}/normal.ontarget.wes_mapq0_bedcoverage.cov"
  --tumor "${BASE}/tumor.ontarget.wes_mapq0_bedcoverage.cov"
  --normalOfftarget "${BASE}/normal.offtarget.wes_100kb_mapq10_bedcoverage.cov"
  --tumorOfftarget "${BASE}/tumor.offtarget.wes_100kb_mapq10_bedcoverage.cov"
  --bed "${BASE}/ssSC_v5.gc.genes.bed"
  --bedOfftarget "${BASE}/clincnv_offtarget_wes_100kb_filtered.bed"
  --pair "${BASE}/pairs_df_filtered.csv"
  --out "${RUN_ROOT}"
  --numberOfThreads "${THREADS}"
  --colNum 4
  --scoreS 200
  --reanalyseCohort
  --lengthS 9
  --filterStep 2
)

printf '%q ' "${CMD[@]}"
printf '\n'
"${CMD[@]}" | tee "${RUN_ROOT}/clinCNV.log"
