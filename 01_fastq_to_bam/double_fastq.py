#!/usr/bin/env python3
import os, re, sys, shutil, subprocess, datetime
from pathlib import Path
from typing import List, Tuple, Optional
import boto3
import botocore
import pandas as pd
import yaml

# =========================
# ======= SETTINGS ========
# =========================

CSV_PATH = "All2FastQSamples.csv"     # must contain columns: sample, s3_dir, files
SNAKEFILE = "snakemake_wf.py"

WORK_ROOT = "./"
LOGS_ROOT = "./run_logs"              # <— central place to keep logs
S3_OUTPUT_PREFIX = "s3://rocken-matched-melanoma-panel-seq/analysis/bam_generation_v3"
SNAKEMAKE_THREADS = 16

# Keep local logs even if jobs fail; never upload logs to S3.
DELETE_WORKDIR_AFTER = True           # if True, delete per-sample workdir AFTER copying logs out

# Use the explicit two FASTQ URIs from the CSV
FASTQ_LIST_COLUMN = "files"

# Required CSV columns
REQUIRED_COLS = {"sample", "s3_dir", "files"}

# Optional: limit to certain samples (None = all)
ONLY_SAMPLES = None

# Reference & targets
REFERENCE_FASTA = "/home/ubuntu/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
REFERENCE_FAI   = "/home/ubuntu/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai"
REFERENCE_DICT  = "/home/ubuntu/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.dict"

TARGETS_BED = "/home/ubuntu/panel_bed/Temp_1_Covered.bed"
TARGETS_PAD = 100

TOOLS = {
    "gatk": "gatk",
    "bwa": "bwa",
    "samtools": "samtools",
    "picard": "picard",
}

# =========================
# ========= CODE ==========
# =========================

def log(msg: str):
    print(f"[run] {msg}", flush=True)

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def s3_split(uri: str) -> Tuple[str, str]:
    if not uri.startswith("s3://"):
        raise ValueError(f"Not an s3 URI: {uri}")
    rest = uri[5:]
    bucket, *key = rest.split("/", 1)
    return bucket, (key[0] if key else "")

def is_fastq_key(key: str) -> bool:
    k = key.lower()
    return k.endswith(".fastq.gz") or k.endswith(".fq.gz")

def classify_r1_r2(name: str) -> Optional[str]:
    n = name.lower()
    if re.search(r'(^|[_\.\-])r1([_\.\-]|$)', n): return "R1"
    if re.search(r'(^|[_\.\-])r2([_\.\-]|$)', n): return "R2"
    if re.search(r'(^|[_\.\-])1\.f(ast)?q\.gz$', n): return "R1"
    if re.search(r'(^|[_\.\-])2\.f(ast)?q\.gz$', n): return "R2"
    return None

def parse_fastq_list_column(row, _s3_bucket_unused: Optional[str]) -> Optional[Tuple[str, List[str], List[str]]]:
    if FASTQ_LIST_COLUMN is None or FASTQ_LIST_COLUMN not in row or pd.isna(row[FASTQ_LIST_COLUMN]):
        return None
    entries = [e.strip() for e in str(row[FASTQ_LIST_COLUMN]).split(";") if e.strip()]
    if not entries: return None
    r1, r2 = [], []
    bucket = None
    for e in entries:
        if not e.startswith("s3://"):
            raise ValueError(f"Expected absolute s3 URI in '{FASTQ_LIST_COLUMN}', found: {e}")
        b, k = s3_split(e)
        bucket = bucket or b
        if b != bucket:
            raise ValueError(f"Mixed buckets in {FASTQ_LIST_COLUMN}: {bucket} vs {b}")
        if not is_fastq_key(k):
            continue
        mate = classify_r1_r2(os.path.basename(k))
        if mate == "R1": r1.append(k)
        elif mate == "R2": r2.append(k)
    r1.sort(); r2.sort()
    if bucket is None: return None
    return bucket, r1, r2

# ---------- S3 helpers ----------

def s3_object_exists(s3, bucket: str, key: str) -> bool:
    try:
        s3.head_object(Bucket=bucket, Key=key)
        return True
    except botocore.exceptions.ClientError as e:
        code = e.response.get("ResponseMetadata", {}).get("HTTPStatusCode")
        if code in (403, 404):
            return False
        raise

def find_key_by_basename(s3, bucket: str, basename: str, hints: List[str]) -> Optional[str]:
    candidate_prefixes = []
    for h in hints:
        if not h: continue
        candidate_prefixes.extend([
            f"{h}/",
            f"fastq/{h}/",
            f"{h}_01/",
            f"{h}_01_S",
        ])
    candidate_prefixes.append("")  # last resort: bucket root

    paginator = s3.get_paginator("list_objects_v2")
    for pref in candidate_prefixes:
        try:
            for page in paginator.paginate(Bucket=bucket, Prefix=pref):
                for obj in page.get("Contents", []):
                    key = obj["Key"]
                    if key.endswith("/"):  # folder marker
                        continue
                    if os.path.basename(key) == basename:
                        return key
        except botocore.exceptions.ClientError:
            continue
    return None

def resolve_or_search_key(s3, bucket: str, key: str, fo_code: Optional[str]) -> Optional[str]:
    if s3_object_exists(s3, bucket, key):
        return key
    base = os.path.basename(key)
    hints = [fo_code or ""]
    alt = find_key_by_basename(s3, bucket, base, hints)
    return alt

# ---------- IO helpers ----------

def download_to(s3, bucket: str, key: str, dest: Path):
    ensure_dir(dest.parent)
    if dest.exists() and dest.stat().st_size > 0:
        log(f"exists: {dest.name}, skip download")
        return
    log(f"Download: s3://{bucket}/{key} -> {dest}")
    s3.download_file(bucket, key, str(dest))

def write_config_yaml(cfg_path: Path, sample_root: Path, sample: str):
    cfg = {
        "input_dir": str(sample_root),
        "output_dir": str(sample_root),
        "samples": [sample],
        "reference": {
            "fasta": REFERENCE_FASTA,
            "fai":   REFERENCE_FAI,
            "dict":  REFERENCE_DICT,
        },
        "targets": {
            "bed": TARGETS_BED,
            "pad": int(TARGETS_PAD),
        },
        "tools": TOOLS,
    }
    with open(cfg_path, "w") as f:
        yaml.safe_dump(cfg, f, sort_keys=False)

def run_snakemake(snakefile: Path, config_yaml: Path, threads: int, sample_log_dir: Path) -> int:
    # capture full snakemake stdout/stderr per sample
    ensure_dir(sample_log_dir)
    ts = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    sm_log = sample_log_dir / f"snakemake_{ts}.log"
    cmd_txt = sample_log_dir / f"snakemake_cmd_{ts}.txt"

    cmd = [
        "snakemake",
        "-s", str(snakefile),
        "--configfile", str(config_yaml),
        "-j", str(threads),
        "--rerun-incomplete",
    ]
    cmd_str = " ".join(cmd)
    with open(cmd_txt, "w") as f:
        f.write(cmd_str + "\n")

    log("Running: " + cmd_str + f"  (logging to {sm_log})")
    with open(sm_log, "wb") as f:
        # stream both stdout and stderr into the same file
        proc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    return proc.returncode

def upload_file(s3, local: Path, s3_uri_prefix: str, sample: str):
    bucket, base_key = s3_split(s3_uri_prefix.rstrip("/") + f"/{sample}")
    key = base_key.rstrip("/") + "/" + local.name
    log(f"Upload: {local} -> s3://{bucket}/{key}")
    s3.upload_file(str(local), bucket, key)

def ensure_exactly_one_pair(r1_keys: List[str], r2_keys: List[str]) -> bool:
    return len(r1_keys) == 1 and len(r2_keys) == 1

def copy_logs_out(sample_root: Path, sample: str):
    """
    Collect logs to LOGS_ROOT/<sample>/ safely:
      - Copy per-rule logs: <sample_root>/<sample>/logs/*
      - Copy metrics folder (often tiny): <sample_root>/<sample>/metrics/*
      - Copy config.yaml (useful context)
    """
    src_logs    = sample_root / sample / "logs"
    src_metrics = sample_root / sample / "metrics"
    dst_root    = Path(LOGS_ROOT) / sample

    ensure_dir(dst_root)
    # copy logs
    if src_logs.exists():
        dst_logs = dst_root / "rule_logs"
        if dst_logs.exists():
            shutil.rmtree(dst_logs)
        shutil.copytree(src_logs, dst_logs)
    # copy metrics (optional, often small)
    if src_metrics.exists():
        dst_metrics = dst_root / "metrics"
        if dst_metrics.exists():
            shutil.rmtree(dst_metrics)
        shutil.copytree(src_metrics, dst_metrics)
    # copy config
    cfg = sample_root / "config.yaml"
    if cfg.exists():
        shutil.copy2(cfg, dst_root / "config.yaml")

def main():
    s3 = boto3.client("s3")
    df = pd.read_csv(CSV_PATH)

    if not REQUIRED_COLS.issubset(df.columns):
        raise SystemExit(f"CSV must contain columns: {', '.join(sorted(REQUIRED_COLS))}")

    if ONLY_SAMPLES:
        df = df[df["sample"].astype(str).isin(ONLY_SAMPLES)]
        if df.empty:
            raise SystemExit("No samples match ONLY_SAMPLES.")

    ensure_dir(Path(WORK_ROOT))
    ensure_dir(Path(LOGS_ROOT))

    for _, row in df.iterrows():
        sample = str(row["sample"]).strip()
        s3_dir = str(row["s3_dir"]).strip()
        fo_code = str(row.get("fo_code", "")).strip()
        if not sample or not s3_dir.startswith("s3://"):
            continue

        log(f"=== Processing {sample} ===")
        sample_root = Path(WORK_ROOT) / sample
        raw_dir     = sample_root / "raw"
        sample_in   = sample_root / sample
        sample_log_dir = Path(LOGS_ROOT) / sample  # where snakemake stdout/stderr goes

        try:
            ensure_dir(raw_dir)
            ensure_dir(sample_in / "logs")
            ensure_dir(sample_in / "metrics")
            ensure_dir(sample_log_dir)

            # Read the two FASTQs from 'files'
            s3_bucket_from_prefix, _ = s3_split(s3_dir)
            chosen = parse_fastq_list_column(row, s3_bucket_from_prefix)
            if chosen is None:
                log(f"{sample}: no FASTQs found in '{FASTQ_LIST_COLUMN}' column; skipping.")
                continue

            bucket, r1_keys, r2_keys = chosen
            log(f"{sample}: CSV-provided FASTQs -> R1:{len(r1_keys)} R2:{len(r2_keys)} (bucket={bucket})")

            if not ensure_exactly_one_pair(r1_keys, r2_keys):
                log(f"{sample}: expected exactly one R1 and one R2; got R1:{len(r1_keys)} R2:{len(r2_keys)}; skipping.")
                continue

            # Resolve keys (verify existence; search by basename if missing)
            r1_key = resolve_or_search_key(s3, bucket, r1_keys[0], fo_code)
            r2_key = resolve_or_search_key(s3, bucket, r2_keys[0], fo_code)
            if r1_key is None or r2_key is None:
                log(f"{sample}: could not locate FASTQs in bucket (R1 found: {r1_key is not None}, R2 found: {r2_key is not None}); skipping.")
                continue

            r1_name = os.path.basename(r1_key)
            r2_name = os.path.basename(r2_key)

            dst_r1 = sample_in / f"{sample}_R1.fastq.gz"
            dst_r2 = sample_in / f"{sample}_R2.fastq.gz"

            # Download to raw/
            download_to(s3, bucket, r1_key, raw_dir / r1_name)
            download_to(s3, bucket, r2_key, raw_dir / r2_name)

            # Copy into canonical filenames expected by pipeline
            if not dst_r1.exists():
                shutil.copy2(raw_dir / r1_name, dst_r1)
            if not dst_r2.exists():
                shutil.copy2(raw_dir / r2_name, dst_r2)

            # Config & run
            cfg_path = sample_root / "config.yaml"
            write_config_yaml(cfg_path, sample_root, sample)

            rc = run_snakemake(Path(SNAKEFILE), cfg_path, SNAKEMAKE_THREADS, sample_log_dir)

            # Regardless of success, collect per-rule logs + config
            copy_logs_out(sample_root, sample)

            if rc != 0:
                log(f"{sample}: Snakemake failed with exit code {rc}; logs retained under {sample_log_dir}")
                # no uploads on failure
            else:
                # Upload only BAM/BAI (not logs)
                out_bam = sample_root / f"{sample}_dedup.bam"
                out_bai = sample_root / f"{sample}_dedup.bai"
                if out_bam.exists():
                    upload_file(s3, out_bam, S3_OUTPUT_PREFIX, sample)
                else:
                    log(f"{sample}: expected BAM not found at {out_bam}")
                if out_bai.exists():
                    upload_file(s3, out_bai, S3_OUTPUT_PREFIX, sample)
                else:
                    log(f"{sample}: expected BAI not found at {out_bai}")

        finally:
            # Cleanup: optionally remove bulky workdir, but logs have been copied to LOGS_ROOT/<sample> already.
            try:
                if DELETE_WORKDIR_AFTER and sample_root.exists():
                    shutil.rmtree(sample_root)
                    log(f"{sample}: cleaned up {sample_root} (logs preserved in {LOGS_ROOT}/{sample})")
                else:
                    log(f"{sample}: kept workdir at {sample_root} (logs also in {LOGS_ROOT}/{sample})")
            except Exception as e:
                log(f"{sample}: cleanup failed: {e}")

        log(f"=== Done {sample} ===")

    log("All samples processed.")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(130)
