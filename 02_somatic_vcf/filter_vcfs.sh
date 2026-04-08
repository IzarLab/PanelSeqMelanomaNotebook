#!/usr/bin/env bash
# existing
set -euo pipefail
# add this:
PS4='+ $(date "+%Y-%m-%d %H:%M:%S") \011${BASH_SOURCE##*/}:\011${LINENO}:\011'
set -x

# ====== USER CONFIG (matches what you shared) ======
GATK_BIN="gatk"                               # or absolute path to GATK
THREADS="${THREADS:-8}"

# Reference
REFERENCE_FASTA="/home/ubuntu/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
REFERENCE_FAI="/home/ubuntu/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai"
REFERENCE_DICT="/home/ubuntu/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.dict"

# Panel targets BED
TARGETS_BED="/home/ubuntu/panel_bed/Temp_1_Covered.bed"

# gnomAD common SNPs (on-target) + index
GNOMAD_VCF="/mnt/myvolume/panel_seq/reset_analysis_sample/cnv_kit_pair/gnomad_common_on_target/gnomad.v3.1.1.common.snps.on_target.merged.vcf.gz"
GNOMAD_VCF_TBI="${GNOMAD_VCF}.tbi"

# S3 locations
S3_VCF_PREFIX="s3://rocken-matched-melanoma-panel-seq/analysis/vcf_gen_v3"
S3_BAM_PREFIX="s3://rocken-matched-melanoma-panel-seq/analysis/bam_generation_v3"

# Local working directories
WORKDIR="${PWD}/gatk_post_mutect2_work"
LOCAL_VCF_BASE="${WORKDIR}/vcf_pairs"
LOCAL_BAM_BASE="${WORKDIR}/bams"
LOGDIR="${WORKDIR}/logs"

mkdir -p "$LOCAL_VCF_BASE" "$LOCAL_BAM_BASE" "$LOGDIR"

# ====== sanity checks ======
for f in "$REFERENCE_FASTA" "$REFERENCE_FAI" "$REFERENCE_DICT" "$TARGETS_BED" "$GNOMAD_VCF" "$GNOMAD_VCF_TBI"; do
  [[ -s "$f" ]] || { echo "Missing required file: $f" >&2; exit 1; }
done

command -v aws >/dev/null || { echo "aws CLI not found in PATH" >&2; exit 1; }
command -v "${GATK_BIN}" >/dev/null || { echo "GATK not found in PATH" >&2; exit 1; }

# ====== helper: fetch a single sample's BAM + BAI from its S3 folder ======
fetch_bam_for_sample () {
  local sample="$1"
  local outdir="${LOCAL_BAM_BASE}/${sample}"
  mkdir -p "$outdir"

  # Pull only *_dedup.bam/.bai (as you described)
  aws s3 cp "${S3_BAM_PREFIX}/${sample}/" "$outdir/" \
    --recursive --exclude "*" --include "*_dedup.bam" --include "*_dedup.bai"

  # Resolve paths
  local bam
  bam="$(ls "$outdir"/*_dedup.bam 2>/dev/null | head -n1 || true)"
  [[ -n "${bam}" && -s "${bam}" ]] || { echo "Could not find dedup BAM for ${sample}" >&2; return 1; }

  # Ensure index
  if [[ ! -s "${bam}.bai" ]]; then
    # try to locate a matching .bai with different basename
    local bai
    bai="$(ls "$outdir"/*.bai 2>/dev/null | head -n1 || true)"
    [[ -s "${bai}" ]] || { echo "Could not find BAI for ${bam}" >&2; return 1; }
  fi

  echo "${bam}"
}

# ====== loop over each {tumor}__{normal} folder in the VCF prefix ======
# We assume you stored Mutect2 outputs as: S3_VCF_PREFIX/{tumor}__{normal}/...
# and that the F1R2 tarball + unfiltered VCF live in that folder.
pairs=$(aws s3 ls "${S3_VCF_PREFIX}/" | awk '{print $2}' | sed 's:/$::' | grep '__' || true)

if [[ -z "${pairs}" ]]; then
  echo "No pair folders found under ${S3_VCF_PREFIX}/" >&2
  exit 1
fi

echo "Found pairs:"
echo "${pairs}"

for pair in $(echo "${pairs}" | tr ' ' '\n' | tac); do
  echo "====== Processing pair: ${pair} ======"
  TUMOR="${pair%%__*}"
  NORMAL="${pair##*__}"

  pair_local="${LOCAL_VCF_BASE}/${pair}"
  mkdir -p "${pair_local}"
  aws s3 sync "${S3_VCF_PREFIX}/${pair}/" "${pair_local}/"

  # Find inputs
  UNFILTERED_VCF="$(ls "${pair_local}"/*unfiltered*.vcf.gz 2>/dev/null | head -n1 || true)"
  [[ -n "${UNFILTERED_VCF}" ]] || UNFILTERED_VCF="$(ls "${pair_local}"/*.vcf.gz 2>/dev/null | head -n1 || true)"
  F1R2_TAR="$(ls "${pair_local}"/*f1r2*.tar.gz "${pair_local}"/*F1R2*.tar.gz 2>/dev/null | head -n1 || true)"

  if [[ -z "${UNFILTERED_VCF}" || ! -s "${UNFILTERED_VCF}" ]]; then
    echo "Missing unfiltered VCF for ${pair}. Skipping." >&2
    continue
  fi
  if [[ -z "${F1R2_TAR}" || ! -s "${F1R2_TAR}" ]]; then
    echo "Missing F1R2 tarball for ${pair}. Skipping." >&2
    continue
  fi

  # ---------- Fetch BAMs (quietly) ----------
  TUMOR_DIR="${LOCAL_BAM_BASE}/${TUMOR}"
  NORMAL_DIR="${LOCAL_BAM_BASE}/${NORMAL}"
  mkdir -p "${TUMOR_DIR}" "${NORMAL_DIR}"

  aws s3 cp "${S3_BAM_PREFIX}/${TUMOR}/" "${TUMOR_DIR}/" \
    --recursive --exclude "*" --include "*_dedup.bam" --include "*_dedup.bai" \
    --only-show-errors >/dev/null || true

  aws s3 cp "${S3_BAM_PREFIX}/${NORMAL}/" "${NORMAL_DIR}/" \
    --recursive --exclude "*" --include "*_dedup.bam" --include "*_dedup.bai" \
    --only-show-errors >/dev/null || true

  # Resolve BAM paths
  TUMOR_BAM="$(ls -1 "${TUMOR_DIR}"/*_dedup.bam 2>/dev/null | head -n1 || true)"
  NORMAL_BAM="$(ls -1 "${NORMAL_DIR}"/*_dedup.bam 2>/dev/null | head -n1 || true)"

  if [[ -z "${TUMOR_BAM}" || ! -s "${TUMOR_BAM}" ]]; then
    echo "Could not find tumor BAM for ${TUMOR} in ${TUMOR_DIR}" >&2; continue
  fi
  if [[ -z "${NORMAL_BAM}" || ! -s "${NORMAL_BAM}" ]]; then
    echo "Could not find normal BAM for ${NORMAL} in ${NORMAL_DIR}" >&2; continue
  fi

  # ---------- Ensure/normalize BAI names ----------
  ensure_bai() {
    local bam="$1"
    local dir; dir="$(dirname "$bam")"
    local bai=""

    if   [[ -s "${bam}.bai" ]]; then
      bai="${bam}.bai"                              # sample.bam.bai
    elif [[ -s "${bam%.bam}.bai" ]]; then
      bai="${bam%.bam}.bai"                         # sample.bai (your S3 layout)
    else
      bai="$(ls -1 "${dir}"/*.bai 2>/dev/null | head -n1 || true)"  # any .bai in dir
    fi

    if [[ -z "${bai}" || ! -s "${bai}" ]]; then
      echo "Missing BAI for ${bam} (looked for ${bam}.bai and ${bam%.bam}.bai)" >&2
      return 1
    fi

    # Normalize so ${bam}.bai always exists (best compatibility with HTSJDK)
    if [[ ! -s "${bam}.bai" ]]; then
      cp -f "${bai}" "${bam}.bai"
    fi
  }

  ensure_bai "${TUMOR_BAM}"  || { echo "Skipping ${pair}"; continue; }
  ensure_bai "${NORMAL_BAM}" || { echo "Skipping ${pair}"; continue; }

  OUTDIR="${pair_local}/post_mutect2"
  mkdir -p "${OUTDIR}"

  LOG="${LOGDIR}/${pair}.post_mutect2.log"
  {
    date
    echo "[${pair}] LearnReadOrientationModel"
    ${GATK_BIN} LearnReadOrientationModel \
      -I "${F1R2_TAR}" \
      -O "${OUTDIR}/tumor.artifact-priors.tar.gz"

    echo "[${pair}] GetPileupSummaries (tumor)"
    echo "+ ${GATK_BIN} GetPileupSummaries -I '${TUMOR_BAM}' -V '${GNOMAD_VCF}' -L '${TARGETS_BED}' -O '${OUTDIR}/tumor.pileups.table'"
    ${GATK_BIN} GetPileupSummaries \
      -I "${TUMOR_BAM}" \
      -V "${GNOMAD_VCF}" \
      -L "${TARGETS_BED}" \
      -O "${OUTDIR}/tumor.pileups.table"

    echo "[${pair}] GetPileupSummaries (normal)"
    echo "+ ${GATK_BIN} GetPileupSummaries -I '${NORMAL_BAM}' -V '${GNOMAD_VCF}' -L '${TARGETS_BED}' -O '${OUTDIR}/normal.pileups.table'"
    ${GATK_BIN} GetPileupSummaries \
      -I "${NORMAL_BAM}" \
      -V "${GNOMAD_VCF}" \
      -L "${TARGETS_BED}" \
      -O "${OUTDIR}/normal.pileups.table"


    echo "[${pair}] CalculateContamination"
    "${GATK_BIN}" CalculateContamination \
      -I "${OUTDIR}/tumor.pileups.table" \
      -matched "${OUTDIR}/normal.pileups.table" \
      -O "${OUTDIR}/tumor.contamination.table" \
      -segments "${OUTDIR}/tumor.contamination.segments"

    echo "[${pair}] FilterMutectCalls"
    "${GATK_BIN}" FilterMutectCalls \
      -V "${UNFILTERED_VCF}" \
      -R "${REFERENCE_FASTA}" \
      --contamination-table "${OUTDIR}/tumor.contamination.table" \
      --tumor-segmentation "${OUTDIR}/tumor.contamination.segments" \
      --ob-priors "${OUTDIR}/tumor.artifact-priors.tar.gz" \
      -O "${OUTDIR}/${pair}.filtered.vcf.gz"

    # Ensure index
    if [[ ! -s "${OUTDIR}/${pair}.filtered.vcf.gz.tbi" ]]; then
      echo "[${pair}] Index filtered VCF (tabix)"
      tabix -p vcf "${OUTDIR}/${pair}.filtered.vcf.gz"
    fi
  } | tee "${LOG}"

  echo "[${pair}] Sync results back to S3"
  aws s3 sync "${OUTDIR}/" "${S3_VCF_PREFIX}/${pair}/post_mutect2/"
  date
  echo "[${pair}] DONE"

  # ====== CLEANUP ======
  echo "[${pair}] Cleaning up local files..."
  rm -rf "${pair_local}"
  rm -rf "${LOCAL_BAM_BASE:?}/${TUMOR}"
  rm -rf "${LOCAL_BAM_BASE:?}/${NORMAL}"
done

echo "All pairs processed."
