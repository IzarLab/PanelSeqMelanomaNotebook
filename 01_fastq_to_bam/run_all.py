#!/usr/bin/env python3
import os
import re
import sys
import shutil
import subprocess
from pathlib import Path
from typing import List, Tuple, Optional

import boto3
import pandas as pd
import yaml

# =========================
# ====== USER CONFIG ======
# =========================
CSV_PATH = "../bam_generation_v3/All4FastQSamples_fix.csv"     # must contain columns: sample, s3_dir
SNAKEFILE = "snakemake_wf.py"

WORK_ROOT = "./"
S3_OUTPUT_PREFIX = "s3://rocken-matched-melanoma-panel-seq/analysis/bam_generation_v3"
SNAKEMAKE_THREADS = 16

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

# Optional: limit to certain samples (None = all)
ONLY_SAMPLES = None  # e.g. ["QHPRG001B4", "QHPRG002BC"]

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

def is_fastq_path(p: Path) -> bool:
    n = p.name.lower()
    return n.endswith(".fastq") or n.endswith(".fastq.gz") or n.endswith(".fq.gz") or n.endswith(".fq")

def is_fastq_key(key: str) -> bool:
    k = key.lower()
    return k.endswith(".fastq.gz") or k.endswith(".fq.gz") or k.endswith(".fastq") or k.endswith(".fq")

def classify_r1_r2(name: str) -> Optional[str]:
    """
    Return 'R1' or 'R2' if recognized; else None.
    Handles patterns like *_R1.fastq.gz, *_L001_R2.fastq.gz, *_1.fastq.gz, *_2.fastq.gz.
    """
    n = name.lower()
    if re.search(r'[_\.]r1([_\.\-]|\.f(ast)?q(\.gz)?$)', n):
        return "R1"
    if re.search(r'[_\.]r2([_\.\-]|\.f(ast)?q(\.gz)?$)', n):
        return "R2"
    if re.search(r'[_\-\.]1\.f(ast)?q(\.gz)?$', n):
        return "R1"
    if re.search(r'[_\-\.]2\.f(ast)?q(\.gz)?$', n):
        return "R2"
    return None

def list_fastqs(s3, s3_prefix: str):
    """Return (bucket, r1_keys, r2_keys) under a given prefix."""
    bucket, prefix = s3_split(s3_prefix.rstrip("/") + "/")
    paginator = s3.get_paginator("list_objects_v2")
    r1, r2 = [], []
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            key = obj["Key"]
            if not is_fastq_key(key):
                continue
            mate = classify_r1_r2(key)
            if mate == "R1":
                r1.append(key)
            elif mate == "R2":
                r2.append(key)
    r1.sort()
    r2.sort()
    return bucket, r1, r2

def download_to(s3, bucket: str, key: str, dest: Path):
    ensure_dir(dest.parent)
    if dest.exists() and dest.stat().st_size > 0:
        log(f"exists: {dest.name}, skip download")
        return
    s3.download_file(bucket, key, str(dest))

def concat_gz_files(srcs: List[Path], dest: Path):
    """
    Concatenate gz members as raw bytes (valid multi-member gzip).
    """
    ensure_dir(dest.parent)
    log(f"Concatenating {len(srcs)} files -> {dest}")
    with open(dest, "wb") as w:
        for s in srcs:
            with open(s, "rb") as r:
                shutil.copyfileobj(r, w)

def write_config_yaml(cfg_path: Path, sample_root: Path, sample: str):
    """
    Write per-sample config.yaml with input_dir=output_dir=<sample_root> and samples:[sample].
    Note: Snakefile you shared does not use targets.bed, but we record it anyway.
    """
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

def run_snakemake(snakefile: Path, config_yaml: Path, threads: int, log_path: Path) -> int:
    cmd = [
        "snakemake",
        "-s", str(snakefile),
        "--configfile", str(config_yaml),
        "-j", str(threads),
        "--rerun-incomplete",
        "--printshellcmds",
    ]
    log("Running: " + " ".join(cmd))
    ensure_dir(log_path.parent)
    with open(log_path, "w") as lf:
        proc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
    return proc.returncode

def upload_file(s3, local: Path, s3_uri_prefix: str, sample: str, dest_subdir: Optional[str] = None):
    """
    Upload local file to s3://<prefix>/<sample>[/<dest_subdir>]/<filename>
    """
    bucket, base_key = s3_split(s3_uri_prefix.rstrip("/") + f"/{sample}")
    if dest_subdir:
        key = base_key.rstrip("/") + f"/{dest_subdir.strip('/')}/{local.name}"
    else:
        key = base_key.rstrip("/") + "/" + local.name
    log(f"Upload: {local} -> s3://{bucket}/{key}")
    s3.upload_file(str(local), bucket, key)

def upload_tree_nonfastq(s3, root: Path, s3_uri_prefix: str, sample: str, dest_subdir: str = "run_artifacts"):
    """
    Recursively upload all files under 'root' EXCEPT FASTQs to
    s3://<prefix>/<sample>/<dest_subdir>/<relative_path>.
    """
    bucket, base_key = s3_split(s3_uri_prefix.rstrip("/") + f"/{sample}")
    for p in sorted(root.rglob("*")):
        if p.is_dir():
            continue
        if is_fastq_path(p):
            continue
        rel = p.relative_to(root)
        key = f"{base_key.rstrip('/')}/{dest_subdir.strip('/')}/{str(rel)}"
        log(f"Upload: {p} -> s3://{bucket}/{key}")
        ensure_dir(p.parent)
        s3.upload_file(str(p), bucket, key)

def safe_rmdir(p: Path):
    if p.exists():
        shutil.rmtree(p)

def main():
    s3 = boto3.client("s3")

    df = pd.read_csv(CSV_PATH)

    need = {"sample", "s3_dir"}
    if not need.issubset(df.columns):
        raise SystemExit(f"CSV must contain columns: {', '.join(sorted(need))}")

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

        snakemake_log = sample_root / "run_snakemake.log"   # capture snakemake stdout+stderr
        driver_log    = sample_root / "driver_run.log"       # we’ll also write a tiny manifest here
        cfg_path      = sample_root / "config.yaml"

        # Make dirs
        ensure_dir(raw_dir)
        ensure_dir(merged_dir)
        ensure_dir(sample_in / "logs")
        ensure_dir(sample_in / "metrics")

        bucket = None
        r1_keys: List[str] = []
        r2_keys: List[str] = []

        try:
            # 1) S3 listing & downloads
            bucket, r1_keys, r2_keys = list_fastqs(s3, s3_dir)
            if not r1_keys or not r2_keys:
                log(f"{sample}: missing mates (R1:{len(r1_keys)}, R2:{len(r2_keys)}); skipping snakemake.")

            else:
                log(f"{sample}: downloading {len(r1_keys)+len(r2_keys)} FASTQs")
                for k in r1_keys + r2_keys:
                    fname = os.path.basename(k)
                    download_to(s3, bucket, k, raw_dir / fname)

                # 2) merge per mate
                r1_parts = [raw_dir / Path(k).name for k in r1_keys]
                r2_parts = [raw_dir / Path(k).name for k in r2_keys]
                merged_r1 = merged_dir / f"{sample}_R1.fastq.gz"
                merged_r2 = merged_dir / f"{sample}_R2.fastq.gz"
                concat_gz_files(r1_parts, merged_r1)
                concat_gz_files(r2_parts, merged_r2)

                # >>> Delete RAW right after creating merged (per request)
                try:
                    log(f"{sample}: deleting raw dir {raw_dir}")
                    safe_rmdir(raw_dir)
                except Exception as e:
                    log(f"{sample}: warning: failed to delete raw dir: {e}")

                # 3) place merged where Snakefile expects them
                dst_r1 = sample_in / f"{sample}_R1.fastq.gz"
                dst_r2 = sample_in / f"{sample}_R2.fastq.gz"
                if not dst_r1.exists():
                    shutil.copy2(merged_r1, dst_r1)
                if not dst_r2.exists():
                    shutil.copy2(merged_r2, dst_r2)

            # 4) per-sample config.yaml (always write)
            write_config_yaml(cfg_path, sample_root, sample)

            # 5) run snakemake (even if FASTQs missing, it will fail early and we capture logs)
            rc = run_snakemake(Path(SNAKEFILE), cfg_path, SNAKEMAKE_THREADS, snakemake_log)

            # 6) upload final outputs (match rule all)
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

            # Always write a tiny manifest/driver log
            with open(driver_log, "w") as dfh:
                dfh.write(f"sample={sample}\n")
                dfh.write(f"s3_dir={s3_dir}\n")
                dfh.write(f"snakemake_exit_code={rc}\n")

            # 7) BEFORE ANY DELETION: upload all NON-FASTQ artifacts to S3 (logs, metrics, yaml, bam intermediates, etc.)
            try:
                upload_tree_nonfastq(s3, sample_root, S3_OUTPUT_PREFIX, sample, dest_subdir="run_artifacts")
            except Exception as e:
                log(f"{sample}: WARNING: upload of non-FASTQ artifacts failed: {e}")

            # 8) Clean up FASTQs written locally (both merged/ and <sample>/*.fastq.gz)
            try:
                # Remove merged fastqs
                if merged_dir.exists():
                    for p in merged_dir.glob("*.fastq*"):
                        try:
                            p.unlink()
                        except Exception as e:
                            log(f"{sample}: warning: could not remove {p}: {e}")

                # Remove sample_in fastqs
                if sample_in.exists():
                    for p in sample_in.glob("*.fastq*"):
                        try:
                            p.unlink()
                        except Exception as e:
                            log(f"{sample}: warning: could not remove {p}: {e}")
            except Exception as e:
                log(f"{sample}: warning during FASTQ cleanup: {e}")

            # 9) Decide whether to keep or delete the whole sample directory
            #    Your request says: "before deleting folder, save logs ... upload all files besides fastqs ... then delete"
            #    We've uploaded the non-FASTQ artifacts already, so now we can safely delete the entire sample_root.
        finally:
            # Final cleanup of the working folder
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
