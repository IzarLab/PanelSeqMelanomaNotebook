#!/usr/bin/env bash
set -euo pipefail

# Annotate filtered panel VCFs with OncoKB without imposing FILTER=PASS.
#
# This consumes the `filtered.*.annotated.vcf[.gz]` files emitted after
# `filter_vcfs.sh`.
#
# The key choice is to keep non-PASS rows that survive the post-Mutect2 VCF
# while applying a modest somatic VAF threshold before OncoKB annotation.
#
# Required environment:
#   export ONCOKB_TOKEN=...
#
# Typical usage:
#   export ONCOKB_TOKEN=...
#   ./run_oncokb_no_pass_filter.sh
#
# Overridable config:
#   VCF_DIR                  directory containing filtered.annotated VCFs
#   OUT_DIR                  directory for per-sample `.oncokb.tsv` outputs
#   ONCOKB_ANNOTATOR_DIR     directory containing `MafAnnotator.py`
#   REF_FASTA                GRCh38 FASTA used for the panel workflow
#   MIN_SOMATIC_VAF          minimum VAF to keep a variant before OncoKB
#   MIN_ALT_READS            optional minimum alternate-read count
#   SLEEP_SECONDS_ON_FAIL    pause after a failed OncoKB API call
#   DEBUG                    set to 1 for shell tracing

VCF_DIR="${VCF_DIR:-/mnt/myvolume/panel_seq/new_bed_analysis/vcfs}"
OUT_DIR="${OUT_DIR:-./oncokb_results_min_no_pass_vaf005}"
ONCOKB_ANNOTATOR_DIR="${ONCOKB_ANNOTATOR_DIR:-/mnt/myvolume/panel_seq/new_bed_analysis/oncokb/oncokb-annotator}"
REF_FASTA="${REF_FASTA:-$HOME/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna}"

MIN_SOMATIC_VAF="${MIN_SOMATIC_VAF:-0.05}"
MIN_ALT_READS="${MIN_ALT_READS:-0}"
SLEEP_SECONDS_ON_FAIL="${SLEEP_SECONDS_ON_FAIL:-8}"
DEBUG="${DEBUG:-0}"

if [[ "${DEBUG}" == "1" ]]; then
  set -x
fi

if [[ -z "${ONCOKB_TOKEN:-}" ]]; then
  echo "ERROR: export ONCOKB_TOKEN before running this script." >&2
  exit 1
fi

command -v bcftools >/dev/null 2>&1 || { echo "ERROR: bcftools not found." >&2; exit 1; }
command -v python >/dev/null 2>&1 || { echo "ERROR: python not found." >&2; exit 1; }
[[ -f "${REF_FASTA}" ]] || { echo "ERROR: REF_FASTA not found at ${REF_FASTA}" >&2; exit 1; }
[[ -f "${REF_FASTA}.fai" ]] || { echo "ERROR: Missing FASTA index ${REF_FASTA}.fai" >&2; exit 1; }
[[ -f "${ONCOKB_ANNOTATOR_DIR}/MafAnnotator.py" ]] || {
  echo "ERROR: MafAnnotator.py not found in ${ONCOKB_ANNOTATOR_DIR}" >&2
  exit 1
}

mkdir -p "${OUT_DIR}"

shopt -s nullglob
mapfile -t VCF_LIST < <(
  find "${VCF_DIR}" -type f \( -name "filtered.*.annotated.vcf" -o -name "filtered.*.annotated.vcf.gz" \) | sort
)

if [[ ${#VCF_LIST[@]} -eq 0 ]]; then
  echo "No filtered annotated VCFs found under ${VCF_DIR}." >&2
  exit 0
fi

echo "Input VCF directory: ${VCF_DIR}"
echo "Output directory:    ${OUT_DIR}"
echo "Min somatic VAF:     ${MIN_SOMATIC_VAF}"
echo "Min alt reads:       ${MIN_ALT_READS}"
echo "FILTER handling:     keep non-PASS rows from the filtered VCF"

for vcf in "${VCF_LIST[@]}"; do
  sample="$(basename "${vcf}")"
  sample="${sample%.annotated.vcf}"
  sample="${sample%.annotated.vcf.gz}"

  maf="${OUT_DIR}/${sample}.minimal.maf"
  output="${OUT_DIR}/${sample}.oncokb.tsv"
  logf="${OUT_DIR}/${sample}.log"
  manifest="${OUT_DIR}/${sample}.variant_manifest.tsv"

  if [[ -s "${output}" ]]; then
    echo "Skipping ${sample}: ${output} already exists."
    continue
  fi

  tumor_from_name="${sample#filtered.}"
  tumor_from_name="${tumor_from_name%%__*}"
  tumor_sample="$(bcftools query -l "${vcf}" | awk -v t="${tumor_from_name}" '$0==t{print; exit}')"
  if [[ -z "${tumor_sample}" ]]; then
    tumor_sample="$(bcftools query -l "${vcf}" | head -n1 || true)"
  fi
  if [[ -z "${tumor_sample}" ]]; then
    echo "Skipping ${sample}: no genotype sample names found." | tee -a "${logf}"
    continue
  fi

  echo "Processing ${sample} using tumor sample ${tumor_sample}"

  total_count="$(bcftools view "${vcf}" | bcftools view -H | wc -l | tr -d ' ')"
  pass_count="$(bcftools view -f PASS "${vcf}" | bcftools view -H | wc -l | tr -d ' ')"
  nonpass_count=$((total_count - pass_count))
  {
    echo "Total variants in filtered VCF (all FILTER states): ${total_count}"
    echo "PASS variants in filtered VCF: ${pass_count}"
    echo "Non-PASS variants retained for this rerun: ${nonpass_count}"
  } > "${logf}"

  hdr="$(bcftools view -h "${vcf}")"
  has_fmt_ad=0
  has_fmt_af=0
  has_info_af=0
  has_info_dp4=0
  grep -q '##FORMAT=<ID=AD' <<<"${hdr}" && has_fmt_ad=1
  grep -q '##FORMAT=<ID=AF' <<<"${hdr}" && has_fmt_af=1
  grep -q '##INFO=<ID=AF' <<<"${hdr}" && has_info_af=1
  grep -q '##INFO=<ID=DP4' <<<"${hdr}" && has_info_dp4=1

  query_fmt='%CHROM\t%POS\t%REF\t%ALT\t%FILTER'
  [[ ${has_fmt_ad} -eq 1 ]] && query_fmt+=$'\t[%AD]'
  [[ ${has_fmt_af} -eq 1 ]] && query_fmt+=$'\t[%AF]'
  [[ ${has_info_af} -eq 1 ]] && query_fmt+=$'\t%INFO/AF'
  [[ ${has_info_dp4} -eq 1 ]] && query_fmt+=$'\t%INFO/DP4'
  query_fmt+=$'\n'

  printf "Chromosome\tStart_Position\tEnd_Position\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tTumor_Sample_Barcode\n" > "${maf}"
  printf "Chromosome\tStart_Position\tEnd_Position\tReference_Allele\tTumor_Seq_Allele2\tTumor_Sample_Barcode\tFILTER\tVAF\n" > "${manifest}"

  if ! bcftools norm -m -both "${vcf}" 2>>"${logf}" \
      | bcftools query -s "${tumor_sample}" -f "${query_fmt}" - 2>>"${logf}" \
      | awk -v OFS='\t' \
             -v sample_name="${sample}" \
             -v min_vaf="${MIN_SOMATIC_VAF}" \
             -v min_alt_reads="${MIN_ALT_READS}" \
             -v maf_path="${maf}" \
             -v manifest_path="${manifest}" \
             -v has_ad="${has_fmt_ad}" \
             -v has_faf="${has_fmt_af}" \
             -v has_iaf="${has_info_af}" \
             -v has_dp4="${has_info_dp4}" '
        function tofloat(x){ if (x=="" || x==".") return -1; return x + 0; }
        function alt_reads_from_ad(ad,   n,a){
          if (!has_ad || ad=="." || ad=="") return -1;
          n = split(ad, a, ",");
          if (n < 2) return -1;
          return a[2] + 0;
        }
        function vaf_from_ad(ad,   n,a,refc,altc){
          if (!has_ad || ad=="." || ad=="") return -1;
          n = split(ad, a, ",");
          if (n < 2) return -1;
          refc = a[1] + 0;
          altc = a[2] + 0;
          if (refc + altc == 0) return -1;
          return altc / (refc + altc);
        }
        function vaf_from_dp4(dp4,   n,a,alt,tot){
          if (!has_dp4 || dp4=="." || dp4=="") return -1;
          n = split(dp4, a, ",");
          if (n < 4) return -1;
          alt = (a[3] + 0) + (a[4] + 0);
          tot = alt + (a[1] + 0) + (a[2] + 0);
          if (tot == 0) return -1;
          return alt / tot;
        }
        {
          chrom = $1;
          pos = $2;
          ref = $3;
          alt = $4;
          filt = $5;
          idx = 6;

          ad = (has_ad ? $(idx) : ".");
          if (has_ad) idx++;
          faf = (has_faf ? $(idx) : ".");
          if (has_faf) idx++;
          iaf = (has_iaf ? $(idx) : ".");
          if (has_iaf) idx++;
          dp4 = (has_dp4 ? $(idx) : ".");

          alt_reads = alt_reads_from_ad(ad);
          vaf = vaf_from_ad(ad);
          if (vaf < 0) vaf = tofloat(faf);
          if (vaf < 0) vaf = tofloat(iaf);
          if (vaf < 0) vaf = vaf_from_dp4(dp4);

          if (vaf < min_vaf) next;
          if (alt_reads >= 0 && alt_reads < min_alt_reads) next;

          start = pos;
          end = pos + length(ref) - 1;
          print chrom, start, end, ref, ref, alt, sample_name >> maf_path;
          print chrom, start, end, ref, alt, sample_name, filt, vaf >> manifest_path;
          kept++;
          if (filt != "PASS") kept_nonpass++;
        }
        END {
          printf("Kept %d variants after VAF >= %.4f and alt reads >= %d\n", kept + 0, min_vaf, min_alt_reads) > "/dev/stderr";
          printf("Of kept variants, %d were non-PASS in the filtered VCF\n", kept_nonpass + 0) > "/dev/stderr";
        }
      ' >> "${logf}" 2>&1
  then
    echo "Failed while building MAF for ${sample}. See ${logf}" >&2
    rm -f "${maf}" "${manifest}"
    continue
  fi

  kept_after="$(($(wc -l < "${maf}") - 1))"
  if [[ "${kept_after}" -le 0 ]]; then
    echo "Skipping ${sample}: no variants survived the VAF/read filters." | tee -a "${logf}"
    rm -f "${maf}" "${manifest}"
    continue
  fi

  if ! python "${ONCOKB_ANNOTATOR_DIR}/MafAnnotator.py" \
        -i "${maf}" \
        -o "${output}" \
        -b "${ONCOKB_TOKEN}" \
        -r GRCh38 \
        -q Genomic_Change \
        -t CM \
        -d >> "${logf}" 2>&1; then
    echo "OncoKB annotation failed for ${sample}. Sleeping ${SLEEP_SECONDS_ON_FAIL}s and continuing." >&2
    rm -f "${output}"
    sleep "${SLEEP_SECONDS_ON_FAIL}"
    continue
  fi

  echo "Finished ${sample} -> ${output}"
done

echo "All samples processed."
