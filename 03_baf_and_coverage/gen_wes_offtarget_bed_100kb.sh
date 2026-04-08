#!/usr/bin/env bash
set -euo pipefail

BASE="${BASE:-/mnt/myvolume/panel_seq/new_bed_analysis}"
ENV_NAME="${ENV_NAME:-mamba-env}"
ENV_BIN="${ENV_BIN:-/home/ubuntu/anaconda3/envs/${ENV_NAME}/bin}"

TARGET_BED="${TARGET_BED:-${BASE}/ssSC_v5.gc.genes.bed}"
CURRENT_OFFTARGET_BED="${CURRENT_OFFTARGET_BED:-${BASE}/clincnv_offtarget_filtered.bed}"
REFERENCE="${REFERENCE:-/home/ubuntu/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna}"
BIN_SIZE="${BIN_SIZE:-100000}"
MIN_BIN_SIZE="${MIN_BIN_SIZE:-50000}"
OUT_BED="${OUT_BED:-${BASE}/clincnv_offtarget_wes_100kb_filtered.bed}"

BEDTOOLS="${BEDTOOLS:-${ENV_BIN}/bedtools}"
SAMTOOLS="${SAMTOOLS:-${ENV_BIN}/samtools}"

if [[ ! -x "${BEDTOOLS}" ]]; then
  echo "ERROR: bedtools not found at ${BEDTOOLS}" >&2
  exit 1
fi

if [[ ! -x "${SAMTOOLS}" ]]; then
  echo "ERROR: samtools not found at ${SAMTOOLS}" >&2
  exit 1
fi

if [[ ! -f "${TARGET_BED}" ]]; then
  echo "ERROR: target BED not found: ${TARGET_BED}" >&2
  exit 1
fi

if [[ ! -f "${CURRENT_OFFTARGET_BED}" ]]; then
  echo "ERROR: current filtered off-target BED not found: ${CURRENT_OFFTARGET_BED}" >&2
  exit 1
fi

if [[ ! -f "${REFERENCE}" ]]; then
  echo "ERROR: reference FASTA not found: ${REFERENCE}" >&2
  exit 1
fi

TMPDIR_WORK="$(mktemp -d -p "${TMPDIR:-/tmp}" wes_offtarget_100kb.XXXXXX)"
trap 'rm -rf "${TMPDIR_WORK}"' EXIT

if [[ ! -f "${REFERENCE}.fai" ]]; then
  "${SAMTOOLS}" faidx "${REFERENCE}"
fi

python - "${CURRENT_OFFTARGET_BED}" "${TMPDIR_WORK}/allowed_contigs.txt" <<'PY'
from pathlib import Path
import sys

in_bed = Path(sys.argv[1])
out_txt = Path(sys.argv[2])
seen = set()
contigs = []
with in_bed.open() as handle:
    for line in handle:
        if not line.strip():
            continue
        chrom = line.split('\t', 1)[0]
        if chrom not in seen:
            seen.add(chrom)
            contigs.append(chrom)

with out_txt.open("w") as handle:
    for chrom in contigs:
        handle.write(f"{chrom}\n")
PY

cut -f1,2 "${REFERENCE}.fai" \
  | awk 'BEGIN{OFS="\t"} NR==FNR { keep[$1]=1; next } keep[$1] { print $1, 0, $2 }' "${TMPDIR_WORK}/allowed_contigs.txt" - \
  > "${TMPDIR_WORK}/genome.bed"

cut -f1-3 "${TARGET_BED}" \
  | sort -k1,1 -k2,2n \
  | "${BEDTOOLS}" merge \
  > "${TMPDIR_WORK}/targets.merged.bed"

"${BEDTOOLS}" subtract \
  -a "${TMPDIR_WORK}/genome.bed" \
  -b "${TMPDIR_WORK}/targets.merged.bed" \
  > "${TMPDIR_WORK}/offtarget.raw.bed"

"${BEDTOOLS}" makewindows \
  -b "${TMPDIR_WORK}/offtarget.raw.bed" \
  -w "${BIN_SIZE}" \
  > "${TMPDIR_WORK}/offtarget.windows.bed"

awk -v min_size="${MIN_BIN_SIZE}" 'BEGIN{OFS="\t"} ($3-$2) >= min_size { print $1, $2, $3 }' \
  "${TMPDIR_WORK}/offtarget.windows.bed" \
  > "${TMPDIR_WORK}/offtarget.windows.filtered.bed"

"${BEDTOOLS}" nuc \
  -fi "${REFERENCE}" \
  -bed "${TMPDIR_WORK}/offtarget.windows.filtered.bed" \
  > "${TMPDIR_WORK}/offtarget.nuc.tsv"

awk 'BEGIN{OFS="\t"}
  NR==1 {
    gc_col = 0
    for (i=1; i<=NF; ++i) {
      field = tolower($i)
      if (field=="gc" || field=="pct_gc" || field ~ /(^|_)pct_gc$/) {
        gc_col = i
      }
    }
    if (gc_col==0) {
      print "ERROR: unable to find GC column in bedtools nuc output" > "/dev/stderr"
      exit 1
    }
    next
  }
  { print $1, $2, $3, $gc_col }' \
  "${TMPDIR_WORK}/offtarget.nuc.tsv" \
  > "${OUT_BED}"

echo "Wrote ${OUT_BED}"
echo "Intervals: $(wc -l < "${OUT_BED}")"
