#!/usr/bin/env python3
import os, re, sys, shutil, subprocess
from pathlib import Path
from typing import List, Tuple, Optional
import boto3
import pandas as pd
import yaml

# =========================
# ======= SETTINGS ========
# =========================

CSV_PATH = "AllDoubleFastQSamples.csv"     # must contain columns: sample, s3_dir
SNAKEFILE = "snakemake_wf.py"

WORK_ROOT = "./"
S3_OUTPUT_PREFIX = "s3://rocken-matched-melanoma-panel-seq/analysis/bam_generation_v3"
SNAKEMAKE_THREADS = 16

# If your CSV includes a column with a semicolon-separated list of FASTQ keys/URIs,
# set its name here. Examples:
#   - absolute URIs:  s3://bucket/…/R1.fastq.gz;s3://bucket/…/R2.fastq.gz
#   - relative keys:  fastq/…/R1.fastq.gz;fastq/…/R2.fastq.gz  (bucket inferred from s3_dir)
FASTQ_LIST_COLUMN = None  # e.g. "fastq_keys" ; leave as None to ignore

# Required CSV columns (always needed)
REQUIRED_COLS = {"sample", "s3_dir"}

# Optional: limit to certain samples (None = all)
ONLY_SAMPLES = None  # e.g. ["QHPRG001B4", "QHPRG002BC"]

# Reference & targets (align with your Snakefile)
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
    """
    Return 'R1' or 'R2' if recognized; else None.
    Handles *_R1.fastq.gz, *_L001_R2.fastq.gz, *_1.fastq.gz, *_2.fastq.gz, BamToFastq_R1_001.fastq.gz, etc.
    """
    n = name.lower()
    if re.search(r'(^|[_\.\-])r1([_\.\-]|$)', n):
        return "R1"
    if re.search(r'(^|[_\.\-])r2([_\.\-]|$)', n):
        return "R2"
    if re.search(r'(^|[_\.\-])1\.f(ast)?q\.gz$', n):
        return "R1"
    if re.search(r'(^|[_\.\-])2\.f(ast)?q\.gz$', n):
        return "R2"
    # common BamToFastq pattern handled by R1/R2 above
    return None

def list_fastqs_under_prefix(s3, s3_prefix: str) -> Tuple[str, List[str], List[str]]:
    """Return (bucket, r1_keys, r2_keys) listed under a given s3 prefix."""
    bucket, prefix = s3_split(s3_prefix.rstrip("/") + "/")
    paginator = s3.get_paginator("list_objects_v2")
    r1, r2 = [], []
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            key = obj["Key"]
            if not is_fastq_key(key):
                continue
            mate = classify_r1_r2(os.path.basename(key))
            if mate == "R1":
                r1.append(key)
            elif mate == "R2":
                r2.append(key)
    r1.sort()
    r2.sort()
    return bucket, r1, r2

def parse_fastq_list_column(row, s3_dir_bucket: Optional[str]) -> Optional[Tuple[str, List[str], List[str]]]:
    """
    If FASTQ_LIST_COLUMN is set and present, parse it.
    Supports absolute s3:// URIs or bare keys (bucket inferred from s3_dir).
    Returns (bucket, r1_keys, r2_keys) or None if unavailable.
    """
    if FASTQ_LIST_COLUMN is None:
        return None
    if FASTQ_LIST_COLUMN not in row or pd.isna(row[FASTQ_LIST_COLUMN]):
        return None

    entries = [e.strip() for e in str(row[FASTQ_LIST_COLUMN]).split(";") if e.strip()]
    if not entries:
        return None

    r1, r2 = [], []
    bucket_for_bare = s3_dir_bucket
    bucket = None

    for e in entries:
        if e.startswith("s3://"):
            b, k = s3_split(e)
            bucket = bucket or b
            if b != bucket:
                raise ValueError(f"Mixed buckets in {FASTQ_LIST_COLUMN}: {bucket} vs {b}")
            key = k
        else:
            # bare key -> needs a bucket from s3_dir
            if not bucket_for_bare:
                raise ValueError(f"Bare key provided but cannot infer bucket: {e}")
            bucket = bucket or bucket_for_bare
            key = e

        if not is_fastq_key(key):
            continue
        mate = classify_r1_r2(os.path.basename(key))
        if mate == "R1":
            r1.append(key)
        elif mate == "R2":
            r2.append(key)

    r1.sort(); r2.sort()
    if bucket is None:
        return None
    return bucket, r1, r2

def download_to(s3, bucket: str, key: str, dest: Path):
    ensure_dir(dest.parent)
    if dest.exists() and dest.stat().st_size > 0:
        log(f"exists: {dest.name}, skip download")
        return
    log(f"Download: s3://{bucket}/{key} -> {dest}")
    s3.download_file(bucket, key, str(dest))

def concat_gz_files(srcs: List[Path], dest: Path):
    """Concatenate gz members as raw bytes (valid multi-member gzip)."""
    ensure_dir(dest.parent)
    if len(srcs) == 1:
        # Fast path: single file -> copy
        shutil.copy2(srcs[0], dest)
        return
    log(f"Concatenating {len(srcs)} files -> {dest}")
    with open(dest, "wb") as w:
        for s in srcs:
            with open(s, "rb") as r:
                shutil.copyfileobj(r, w)

def write_config_yaml(cfg_path: Path, sample_root: Path, sample: str):
    """Write per-sample config.yaml with input_dir=output_dir=<sample_root> and samples:[sample]."""
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

def run_snakemake(snakefile: Path, config_yaml: Path, threads: int) -> int:
    cmd = [
        "snakemake",
        "-s", str(snakefile),
        "--configfile", str(config_yaml),
        "-j", str(threads),
        "--rerun-incomplete",
    ]
    log("Running: " + " ".join(cmd))
    return subprocess.run(cmd).returncode

def upload_file(s3, local: Path, s3_uri_prefix: str, sample: str):
    bucket, base_key = s3_split(s3_uri_prefix.rstrip("/") + f"/{sample}")
    key = base_key.rstrip("/") + "/" + local.name
    log(f"Upload: {local} -> s3://{bucket}/{key}")
    s3.upload_file(str(local), bucket, key)

def ensure_two_mates(r1_keys: List[str], r2_keys: List[str]) -> bool:
    """Return True if we have at least 1 R1 and 1 R2; otherwise False."""
    return len(r1_keys) >= 1 and len(r2_keys) >= 1

def main():
    s3 = boto3.client("s3")

    df = pd.read_csv(CSV_PATH)
    # If your CSV has intro lines to skip, set your own .iloc slicing here.
    # df = df.iloc[2:].reset_index(drop=True)

    if not REQUIRED_COLS.issubset(df.columns):
        raise SystemExit(f"CSV must contain columns: {', '.join(sorted(REQUIRED_COLS))}")

    if ONLY_SAMPLES:
        df = df[df["sample"].astype(str).isin(ONLY_SAMPLES)]
        if df.empty:
            raise SystemExit("No samples match ONLY_SAMPLES.")

    ensure_dir(Path(WORK_ROOT))

    for _, row in df.iterrows():
        sample = str(row["sample"]).strip()
        s3_dir = str(row["s3_dir"]).strip()
        if not sample or not s3_dir.startswith("s3://"):
            continue

        log(f"=== Processing {sample} ===")
        sample_root = Path(WORK_ROOT) / sample
        raw_dir     = sample_root / "raw"
        merged_dir  = sample_root / "merged"
        sample_in   = sample_root / sample  # IN_DIR/<sample> (Snakefile expects fastqs here)
        try:
            ensure_dir(raw_dir)
            ensure_dir(merged_dir)
            ensure_dir(sample_in / "logs")
            ensure_dir(sample_in / "metrics")

            # --- 1) Determine FASTQ keys (CSV list preferred, else list under prefix) ---
            s3_bucket_from_prefix, _ = s3_split(s3_dir)
            chosen = parse_fastq_list_column(row, s3_bucket_from_prefix)
            if chosen is not None:
                bucket, r1_keys, r2_keys = chosen
                log(f"{sample}: using FASTQ list from '{FASTQ_LIST_COLUMN}' (bucket={bucket})")
            else:
                bucket, r1_keys, r2_keys = list_fastqs_under_prefix(s3, s3_dir)
                log(f"{sample}: discovered under prefix -> R1:{len(r1_keys)} R2:{len(r2_keys)}")

            if not ensure_two_mates(r1_keys, r2_keys):
                log(f"{sample}: missing mates (R1:{len(r1_keys)}, R2:{len(r2_keys)}); skipping.")
                continue

            # --- 2) Download all found FASTQs ---
            all_keys = r1_keys + r2_keys
            log(f"{sample}: downloading {len(all_keys)} FASTQs")
            for k in all_keys:
                fname = os.path.basename(k)
                download_to(s3, bucket, k, raw_dir / fname)

            # --- 3) Prepare per-mate inputs (skip merge if exactly one each) ---
            r1_parts = [raw_dir / os.path.basename(k) for k in r1_keys]
            r2_parts = [raw_dir / os.path.basename(k) for k in r2_keys]

            dst_r1 = sample_in / f"{sample}_R1.fastq.gz"
            dst_r2 = sample_in / f"{sample}_R2.fastq.gz"

            if len(r1_parts) == 1 and len(r2_parts) == 1:
                # Fast path: exactly one R1 and one R2 -> just copy to destination
                log(f"{sample}: single-pair FASTQs detected; skipping merge.")
                if not dst_r1.exists():
                    shutil.copy2(r1_parts[0], dst_r1)
                if not dst_r2.exists():
                    shutil.copy2(r2_parts[0], dst_r2)
            else:
                # Merge multiple parts per mate into canonical files
                merged_r1 = merged_dir / f"{sample}_R1.fastq.gz"
                merged_r2 = merged_dir / f"{sample}_R2.fastq.gz"
                concat_gz_files(r1_parts, merged_r1)
                concat_gz_files(r2_parts, merged_r2)
                if not dst_r1.exists():
                    shutil.copy2(merged_r1, dst_r1)
                if not dst_r2.exists():
                    shutil.copy2(merged_r2, dst_r2)

            # --- 4) per-sample Snakemake config ---
            cfg_path = sample_root / "config.yaml"
            write_config_yaml(cfg_path, sample_root, sample)

            # --- 5) run snakemake ---
            rc = run_snakemake(Path(SNAKEFILE), cfg_path, SNAKEMAKE_THREADS)
            if rc != 0:
                log(f"{sample}: Snakemake failed with exit code {rc}; skipping upload.")
                continue

            # --- 6) upload final outputs ---
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
            # --- 7) cleanup ---
            try:
                if sample_root.exists():
                    shutil.rmtree(sample_root)
                log(f"{sample}: cleaned up {sample_root}")
            except Exception as e:
                log(f"{sample}: cleanup failed: {e}")

        log(f"=== Done {sample} ===")

    log("All samples processed.")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(130)
