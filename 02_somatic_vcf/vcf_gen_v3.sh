#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# User config (edit as needed)
# -----------------------------
PAIRS_FILE="pairs.txt"

# References
REFERENCE_FASTA="/home/ubuntu/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
REFERENCE_FAI="/home/ubuntu/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai"
REFERENCE_DICT="/home/ubuntu/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.dict"

# Targets BED (panel)
TARGETS_BED="/home/ubuntu/panel_bed/Temp_1_Covered.bed"

# Panel of normals
PON_VCF="~/gen_bam/pon.vcf.gz"

# Intervals (uncomment ONE)
#INTERVALS_ARGS="-L ${TARGETS_BED}"                           # no padding
INTERVALS_ARGS="-L ${TARGETS_BED} --interval-padding 100"      # typical padding for panels
#INTERVALS_ARGS=""                                             # no interval restriction

# S3 locations
S3_BAM_ROOT="s3://rocken-matched-melanoma-panel-seq/analysis/bam_generation_v3"
S3_VCF_ROOT="s3://rocken-matched-melanoma-panel-seq/analysis/vcf_gen_v3"

# Local working dir (scratch)
WORKDIR_BASE="${PWD}/mutect2_work"

# GATK binary (ensure in PATH, or set full path)
GATK_BIN="gatk"

# Threads for Mutect2 (optional)
THREADS=8

# -----------------------------
# Pre-flight checks
# -----------------------------
command -v aws >/dev/null 2>&1 || { echo "ERROR: aws CLI not found."; exit 1; }
command -v "${GATK_BIN}" >/dev/null 2>&1 || { echo "ERROR: gatk not found."; exit 1; }

for f in "${REFERENCE_FASTA}" "${REFERENCE_FAI}" "${REFERENCE_DICT}"; do
  [[ -s "${f}" ]] || { echo "ERROR: Missing reference file ${f}"; exit 1; }
done

if [[ -n "${TARGETS_BED}" && "${INTERVALS_ARGS:-}" != "" ]]; then
  [[ -s "${TARGETS_BED}" ]] || { echo "ERROR: Missing TARGETS_BED at ${TARGETS_BED}"; exit 1; }
fi

[[ -s "${PAIRS_FILE}" ]] || { echo "ERROR: pairs.txt not found at ${PAIRS_FILE}"; exit 1; }

# Expand ~ for PON (shell-safe)
PON_VCF_EXPANDED=$(eval echo "${PON_VCF}")
[[ -s "${PON_VCF_EXPANDED}" ]] || { echo "ERROR: PON VCF not found at ${PON_VCF_EXPANDED}"; exit 1; }
[[ -s "${PON_VCF_EXPANDED}.tbi" ]] || { echo "ERROR: PON index not found at ${PON_VCF_EXPANDED}.tbi"; exit 1; }

mkdir -p "${WORKDIR_BASE}"

# -----------------------------
# Main loop
# -----------------------------
while IFS=, read -r TUMOR NORMAL; do
  # Skip empty/comment lines
  [[ -z "${TUMOR}" || "${TUMOR}" =~ ^# ]] && continue
  [[ -z "${NORMAL}" ]] && { echo "WARN: No normal for ${TUMOR}, skipping."; continue; }

  echo "=== Processing tumor=${TUMOR} normal=${NORMAL} ==="

  PAIR_DIR="${WORKDIR_BASE}/${TUMOR}__${NORMAL}"
  mkdir -p "${PAIR_DIR}"
  pushd "${PAIR_DIR}" >/dev/null

  # S3 paths to BAMs
  TUMOR_S3_BAM="${S3_BAM_ROOT}/${TUMOR}/dedup.bam"
  TUMOR_S3_BAI="${S3_BAM_ROOT}/${TUMOR}/dedup.bai"
  NORMAL_S3_BAM="${S3_BAM_ROOT}/${NORMAL}/dedup.bam"
  NORMAL_S3_BAI="${S3_BAM_ROOT}/${NORMAL}/dedup.bai"

  # Local filenames
  TUMOR_BAM="tumor.dedup.bam"
  TUMOR_BAI="tumor.dedup.bai"
  NORMAL_BAM="normal.dedup.bam"
  NORMAL_BAI="normal.dedup.bai"

  echo "[S3->local] Downloading BAMs…"
  aws s3 cp "${TUMOR_S3_BAM}"   "${TUMOR_BAM}"
  aws s3 cp "${TUMOR_S3_BAI}"   "${TUMOR_BAI}"
  aws s3 cp "${NORMAL_S3_BAM}"  "${NORMAL_BAM}"
  aws s3 cp "${NORMAL_S3_BAI}"  "${NORMAL_BAI}"

  # Output basenames
  OUT_PREFIX="${TUMOR}__${NORMAL}"
  OUT_VCF="unfiltered.${OUT_PREFIX}.vcf.gz"
  OUT_F1R2="f1r2.${OUT_PREFIX}.tar.gz"
  OUT_LOG="mutect2.${OUT_PREFIX}.log.txt"

  echo "[Mutect2] Running..."
  # shellcheck disable=SC2086
  "${GATK_BIN}" Mutect2 \
    -R "${REFERENCE_FASTA}" \
    -I "${TUMOR_BAM}"  -tumor "${TUMOR}" \
    -I "${NORMAL_BAM}" -normal "${NORMAL}" \
    --panel-of-normals "${PON_VCF_EXPANDED}" \
    --af-of-alleles-not-in-resource 2.5e-6 \
    --f1r2-tar-gz "${OUT_F1R2}" \
    ${INTERVALS_ARGS:-} \
    --native-pair-hmm-threads "${THREADS}" \
    -O "${OUT_VCF}" \
    2>&1 | tee "${OUT_LOG}"

  echo "[Index check] Ensuring VCF index..."
  if [[ ! -s "${OUT_VCF}.tbi" ]]; then
    "${GATK_BIN}" IndexFeatureFile -I "${OUT_VCF}"
  fi

  # Upload results to S3
  S3_OUT_DIR="${S3_VCF_ROOT}/${OUT_PREFIX}"
  echo "[local->S3] Uploading outputs to ${S3_OUT_DIR}"
  aws s3 cp "${OUT_VCF}"        "${S3_OUT_DIR}/"
  aws s3 cp "${OUT_VCF}.tbi"    "${S3_OUT_DIR}/"
  aws s3 cp "${OUT_F1R2}"       "${S3_OUT_DIR}/"
  aws s3 cp "${OUT_LOG}"        "${S3_OUT_DIR}/"
  aws s3 cp "${OUT_VCF%.vcf}.stats" "${S3_OUT_DIR}/"

  
  aws s3 cp . "${S3_OUT_DIR}/" --recursive \
  --exclude "*" \
  --include "$(basename "${OUT_VCF}")" \
  --include "$(basename "${OUT_VCF}").tbi" \
  --include "$(basename "${OUT_VCF}").stats" \
  --include "$(basename "${OUT_F1R2}")" \
  --include "$(basename "${OUT_LOG}")"


  # Clean up local BAMs
  echo "[cleanup] Removing local BAMs..."
  rm -f "${TUMOR_BAM}" "${TUMOR_BAI}" "${NORMAL_BAM}" "${NORMAL_BAI}"

  popd >/dev/null
  echo "=== Done ${TUMOR}__${NORMAL} ==="
done < "${PAIRS_FILE}"

echo "All pairs processed."
