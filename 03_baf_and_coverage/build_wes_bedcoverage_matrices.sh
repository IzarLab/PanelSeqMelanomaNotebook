#!/usr/bin/env bash
set -euo pipefail

BASE="${BASE:-/mnt/myvolume/panel_seq/new_bed_analysis}"
ENV_NAME="${ENV_NAME:-mamba-env}"
ENV_BIN="${ENV_BIN:-/home/ubuntu/anaconda3/envs/${ENV_NAME}/bin}"
S3_PREFIX="${S3_PREFIX:-s3://rocken-matched-melanoma-panel-seq/analysis/bam_generation_v3}"
TMP_ROOT="${TMP_ROOT:-/home/ubuntu/temp}"
CACHE_ROOT="${CACHE_ROOT:-${BASE}/bedcoverage_wes_cache}"
LOG_ROOT="${LOG_ROOT:-${BASE}/bedcoverage_wes_logs}"
THREADS="${THREADS:-8}"

TARGET_BED="${TARGET_BED:-${BASE}/ssSC_v5.gc.genes.bed}"
OFFTARGET_BED="${OFFTARGET_BED:-${BASE}/clincnv_offtarget_wes_100kb_filtered.bed}"
TARGET_BED_SORTED="${TARGET_BED_SORTED:-${BASE}/ssSC_v5.gc.genes.bedcoverage.sorted.bed}"
OFFTARGET_BED_SORTED="${OFFTARGET_BED_SORTED:-${BASE}/clincnv_offtarget_wes_100kb_filtered.sorted.bed}"

TUMOR_HEADER_SOURCE="${TUMOR_HEADER_SOURCE:-${BASE}/tumor.ontarget.cov}"
NORMAL_HEADER_SOURCE="${NORMAL_HEADER_SOURCE:-${BASE}/normal.ontarget.cov}"

TUMOR_ONTARGET_OUT="${TUMOR_ONTARGET_OUT:-${BASE}/tumor.ontarget.wes_mapq0_bedcoverage.cov}"
NORMAL_ONTARGET_OUT="${NORMAL_ONTARGET_OUT:-${BASE}/normal.ontarget.wes_mapq0_bedcoverage.cov}"
TUMOR_OFFTARGET_OUT="${TUMOR_OFFTARGET_OUT:-${BASE}/tumor.offtarget.wes_100kb_mapq10_bedcoverage.cov}"
NORMAL_OFFTARGET_OUT="${NORMAL_OFFTARGET_OUT:-${BASE}/normal.offtarget.wes_100kb_mapq10_bedcoverage.cov}"

BEDCOVERAGE="${BEDCOVERAGE:-${ENV_BIN}/BedCoverage}"
AWS="${AWS:-/usr/bin/aws}"
OFFTARGET_SCRIPT="${OFFTARGET_SCRIPT:-${BASE}/gen_wes_offtarget_bed_100kb.sh}"
REFERENCE="${REFERENCE:-/home/ubuntu/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna}"

mkdir -p "${TMP_ROOT}" "${CACHE_ROOT}/ontarget" "${CACHE_ROOT}/offtarget" "${LOG_ROOT}"

if [[ ! -x "${BEDCOVERAGE}" ]]; then
  echo "ERROR: BedCoverage not found at ${BEDCOVERAGE}" >&2
  exit 1
fi

if [[ ! -x "${AWS}" ]]; then
  echo "ERROR: aws CLI not found at ${AWS}" >&2
  exit 1
fi

if [[ ! -f "${TARGET_BED}" ]]; then
  echo "ERROR: target BED not found: ${TARGET_BED}" >&2
  exit 1
fi

if [[ ! -f "${OFFTARGET_BED}" ]]; then
  "${OFFTARGET_SCRIPT}"
fi

prepare_sorted_bed() {
  local in_bed="$1"
  local out_bed="$2"

  python - "${REFERENCE}.fai" "${in_bed}" "${out_bed}" <<'PY'
from pathlib import Path
import sys

fai_path = Path(sys.argv[1])
in_bed = Path(sys.argv[2])
out_bed = Path(sys.argv[3])

chrom_order = {}
with fai_path.open() as handle:
    for idx, line in enumerate(handle):
        chrom = line.split('\t', 1)[0]
        chrom_order[chrom] = idx

rows = []
with in_bed.open() as handle:
    for line in handle:
        if not line.strip():
            continue
        parts = line.rstrip('\n').split('\t')
        chrom = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        rows.append((chrom_order.get(chrom, 10**9), chrom, start, end, parts))

rows.sort(key=lambda item: (item[0], item[2], item[3], item[1]))

with out_bed.open("w") as handle:
    for _, _, _, _, parts in rows:
        handle.write('\t'.join(parts) + '\n')
PY
}

prepare_sorted_bed "${TARGET_BED}" "${TARGET_BED_SORTED}"
prepare_sorted_bed "${OFFTARGET_BED}" "${OFFTARGET_BED_SORTED}"

RUN_TS="$(date -u +%Y%m%dT%H%M%SZ)"
SUMMARY_LOG="${LOG_ROOT}/wes_bedcoverage_${RUN_TS}.log"
MISSING_LOG="${LOG_ROOT}/wes_bedcoverage_missing_samples.tsv"
FAILED_LOG="${LOG_ROOT}/wes_bedcoverage_failed_samples.tsv"

touch "${MISSING_LOG}" "${FAILED_LOG}"

ORDER_DIR="${CACHE_ROOT}/sample_orders"
mkdir -p "${ORDER_DIR}"
TUMOR_ORDER_FILE="${ORDER_DIR}/tumor_samples.txt"
NORMAL_ORDER_FILE="${ORDER_DIR}/normal_samples.txt"
UNION_ORDER_FILE="${ORDER_DIR}/all_samples.txt"

python - "${TUMOR_HEADER_SOURCE}" "${NORMAL_HEADER_SOURCE}" "${TUMOR_ORDER_FILE}" "${NORMAL_ORDER_FILE}" "${UNION_ORDER_FILE}" <<'PY'
from pathlib import Path
import sys

tumor_src = Path(sys.argv[1])
normal_src = Path(sys.argv[2])
tumor_out = Path(sys.argv[3])
normal_out = Path(sys.argv[4])
union_out = Path(sys.argv[5])

def header_samples(path: Path):
    with path.open() as handle:
        return handle.readline().rstrip("\n").split("\t")[3:]

tumors = header_samples(tumor_src)
normals = header_samples(normal_src)

seen = set()
union = []
for sample in tumors + normals:
    if sample not in seen:
        seen.add(sample)
        union.append(sample)

tumor_out.write_text("".join(f"{sample}\n" for sample in tumors))
normal_out.write_text("".join(f"{sample}\n" for sample in normals))
union_out.write_text("".join(f"{sample}\n" for sample in union))
PY

echo "Logging to ${SUMMARY_LOG}" | tee -a "${SUMMARY_LOG}"
echo "Using cache root ${CACHE_ROOT}" | tee -a "${SUMMARY_LOG}"

download_and_cover() {
  local sample="$1"
  local sample_dir="${S3_PREFIX}/${sample}/"
  local bam_s3=""
  local bai_s3=""
  local bam_file=""
  local bai_file=""
  local bam_local="${TMP_ROOT}/${sample}_dedup.bam"
  local bai_local="${TMP_ROOT}/${sample}_dedup.bam.bai"
  local on_tmp="${TMP_ROOT}/${sample}.ontarget.bed"
  local off_tmp="${TMP_ROOT}/${sample}.offtarget.bed"
  local on_cache="${CACHE_ROOT}/ontarget/${sample}.bed"
  local off_cache="${CACHE_ROOT}/offtarget/${sample}.bed"

  if [[ -s "${on_cache}" && -s "${off_cache}" ]]; then
    echo "[skip] ${sample}: cached outputs already exist" | tee -a "${SUMMARY_LOG}"
    return 0
  fi

  bam_file="$("${AWS}" s3 ls "${sample_dir}" | awk '/_dedup\.bam$/ {print $3 "\t" $4}' | sort -nr | head -1 | cut -f2)"
  bai_file="$("${AWS}" s3 ls "${sample_dir}" | awk '/_dedup\.bai$/ {print $4}' | head -1 || true)"

  if [[ -z "${bam_file}" ]]; then
    echo -e "${sample}\tmissing_bam\t${sample_dir}" | tee -a "${MISSING_LOG}" "${SUMMARY_LOG}"
    return 1
  fi

  bam_s3="${sample_dir}${bam_file}"
  echo "[copy] ${sample}: ${bam_s3}" | tee -a "${SUMMARY_LOG}"
  "${AWS}" s3 cp --no-progress "${bam_s3}" "${bam_local}"

  if [[ -n "${bai_file}" ]]; then
    bai_s3="${sample_dir}${bai_file}"
    "${AWS}" s3 cp --no-progress "${bai_s3}" "${bai_local}"
  else
    echo "[warn] ${sample}: no *_dedup.bai found, continuing with BAM only" | tee -a "${SUMMARY_LOG}"
  fi

  if [[ ! -s "${on_cache}" ]]; then
    echo "[cover] ${sample}: on-target MAPQ 0" | tee -a "${SUMMARY_LOG}"
    if ! "${BEDCOVERAGE}" \
      -bam "${bam_local}" \
      -in "${TARGET_BED_SORTED}" \
      -clear \
      -min_mapq 0 \
      -decimals 6 \
      -threads "${THREADS}" \
      -out "${on_tmp}" >> "${SUMMARY_LOG}" 2>&1; then
      echo "[error] ${sample}: on-target BedCoverage failed" | tee -a "${SUMMARY_LOG}"
      return 1
    fi
    if [[ ! -s "${on_tmp}" ]]; then
      echo "[error] ${sample}: expected on-target output missing: ${on_tmp}" | tee -a "${SUMMARY_LOG}"
      return 1
    fi
    mv "${on_tmp}" "${on_cache}"
  else
    echo "[skip] ${sample}: on-target cache already exists" | tee -a "${SUMMARY_LOG}"
  fi

  if [[ ! -s "${off_cache}" ]]; then
    echo "[cover] ${sample}: off-target MAPQ 10" | tee -a "${SUMMARY_LOG}"
    if ! "${BEDCOVERAGE}" \
      -bam "${bam_local}" \
      -in "${OFFTARGET_BED_SORTED}" \
      -clear \
      -min_mapq 10 \
      -decimals 6 \
      -threads "${THREADS}" \
      -out "${off_tmp}" >> "${SUMMARY_LOG}" 2>&1; then
      echo "[error] ${sample}: off-target BedCoverage failed" | tee -a "${SUMMARY_LOG}"
      return 1
    fi
    if [[ ! -s "${off_tmp}" ]]; then
      echo "[error] ${sample}: expected off-target output missing: ${off_tmp}" | tee -a "${SUMMARY_LOG}"
      return 1
    fi
    mv "${off_tmp}" "${off_cache}"
  else
    echo "[skip] ${sample}: off-target cache already exists" | tee -a "${SUMMARY_LOG}"
  fi

  rm -f "${bam_local}" "${bai_local}"

  echo "[done] ${sample}" | tee -a "${SUMMARY_LOG}"
  return 0
}

cleanup_sample_files() {
  local sample="$1"
  rm -f \
    "${TMP_ROOT}/${sample}_dedup.bam" \
    "${TMP_ROOT}/${sample}_dedup.bam.bai" \
    "${TMP_ROOT}/${sample}.ontarget.bed" \
    "${TMP_ROOT}/${sample}.offtarget.bed"
}

while IFS= read -r sample; do
  [[ -z "${sample}" ]] && continue
  cleanup_sample_files "${sample}"
  if ! download_and_cover "${sample}"; then
    echo -e "${sample}\tprocessing_failed" >> "${FAILED_LOG}"
    cleanup_sample_files "${sample}"
    exit 1
  fi
  cleanup_sample_files "${sample}"
done < "${UNION_ORDER_FILE}"

merge_matrix() {
  local bed_path="$1"
  local order_path="$2"
  local cache_dir="$3"
  local out_path="$4"

  python - "${bed_path}" "${order_path}" "${cache_dir}" "${out_path}" <<'PY'
from pathlib import Path
import sys

bed_path = Path(sys.argv[1])
order_path = Path(sys.argv[2])
cache_dir = Path(sys.argv[3])
out_path = Path(sys.argv[4])

with order_path.open() as handle:
    samples = [line.strip() for line in handle if line.strip()]

bed_rows = []
with bed_path.open() as handle:
    for line in handle:
        if not line.strip():
            continue
        if line.startswith("#"):
            continue
        bed_rows.append(tuple(line.rstrip("\n").split("\t")[:3]))

coverage_columns = []
for sample in samples:
    sample_path = cache_dir / f"{sample}.bed"
    if not sample_path.exists():
        raise SystemExit(f"Missing cached coverage for {sample}: {sample_path}")

    rows = []
    with sample_path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                raise SystemExit(f"Coverage file has fewer than 4 columns: {sample_path}")
            rows.append(parts)

    coords = [tuple(parts[:3]) for parts in rows]
    if coords != bed_rows:
        raise SystemExit(f"Coordinate mismatch for {sample}: {sample_path}")

    coverage_columns.append([parts[3] for parts in rows])

with out_path.open("w") as out_handle:
    out_handle.write("\t".join(["chr", "start", "end", *samples]) + "\n")
    for idx, coords in enumerate(bed_rows):
        values = [column[idx] for column in coverage_columns]
        out_handle.write("\t".join([*coords, *values]) + "\n")
PY
}

merge_matrix "${TARGET_BED_SORTED}" "${TUMOR_ORDER_FILE}" "${CACHE_ROOT}/ontarget" "${TUMOR_ONTARGET_OUT}"
merge_matrix "${TARGET_BED_SORTED}" "${NORMAL_ORDER_FILE}" "${CACHE_ROOT}/ontarget" "${NORMAL_ONTARGET_OUT}"
merge_matrix "${OFFTARGET_BED_SORTED}" "${TUMOR_ORDER_FILE}" "${CACHE_ROOT}/offtarget" "${TUMOR_OFFTARGET_OUT}"
merge_matrix "${OFFTARGET_BED_SORTED}" "${NORMAL_ORDER_FILE}" "${CACHE_ROOT}/offtarget" "${NORMAL_OFFTARGET_OUT}"

echo "Finished generating WES BedCoverage matrices:" | tee -a "${SUMMARY_LOG}"
echo "  ${TUMOR_ONTARGET_OUT}" | tee -a "${SUMMARY_LOG}"
echo "  ${NORMAL_ONTARGET_OUT}" | tee -a "${SUMMARY_LOG}"
echo "  ${TUMOR_OFFTARGET_OUT}" | tee -a "${SUMMARY_LOG}"
echo "  ${NORMAL_OFFTARGET_OUT}" | tee -a "${SUMMARY_LOG}"
