from __future__ import annotations

import re
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd


# Build the Figure 1 aneuploidy score inputs from the combined s100_l3_f1
# ClinCNV run.
PANEL_SEQ_ROOT = Path("/mnt/myvolume/panel_seq")
ROOT = PANEL_SEQ_ROOT / "new_bed_analysis/hrd_input"
COMBINED_RUN_DIR = (
    PANEL_SEQ_ROOT
    / "v4_v5_clincnv/new_clincnv_runs/combined_runs/combined_bed_split_s100_l3_f1"
)
COMBINED_MANIFEST = COMBINED_RUN_DIR / "combined_file_manifest.tsv"
OUTPUT_DIR = ROOT / "hrd_inputs"
SEGMENT_DIR = OUTPUT_DIR / "aneuploidy_segments"
RESULTS_DIR = OUTPUT_DIR / "aneuploidy_results"
CYTOBAND_FILE = PANEL_SEQ_ROOT / "gistic/hg38_cytoBandIdeo.txt"
AUTOSOMES = {f"chr{i}" for i in range(1, 23)}
THRESHOLD = 0.5


def parse_header_value(lines, key):
    pattern = re.compile(rf"^##{re.escape(key)}:\s*(.+?)\s*$")
    for line in lines:
        match = pattern.match(line)
        if match:
            return match.group(1)
    return None


def clean_chr(value):
    value = str(value).strip()
    return value if value.startswith("chr") else f"chr{value}"


def read_clincnv_table(path):
    lines = path.read_text().splitlines()
    meta = {
        "ploidy": parse_header_value(lines, "ploidy"),
        "sex": parse_header_value(lines, "gender of sample"),
    }

    header_idx = None
    for idx, line in enumerate(lines):
        if line.startswith("#chr\t") or line.startswith("#chrom\t"):
            header_idx = idx
            break

    if header_idx is None:
        raise ValueError(f"No ClinCNV header row found in {path}")

    header = lines[header_idx].lstrip("#")
    data_lines = [header] + lines[header_idx + 1 :]
    data_lines = [line for line in data_lines if line.strip()]
    if len(data_lines) <= 1:
        return meta, pd.DataFrame()

    df = pd.read_csv(StringIO("\n".join(data_lines)), sep="\t", dtype=str)
    return meta, df


def build_arm_table(cytoband_file=CYTOBAND_FILE):
    cytoband = pd.read_csv(
        cytoband_file,
        sep="\t",
        header=None,
        names=["chrom", "chromStart", "chromEnd", "band", "stain"],
        dtype={"chrom": str, "band": str},
    )
    cytoband = cytoband[cytoband["chrom"].isin(AUTOSOMES)].copy()
    cytoband["arm_simple"] = cytoband["band"].str[0]
    cytoband = cytoband[cytoband["arm_simple"].isin(["p", "q"])].copy()

    arm_df = (
        cytoband.groupby(["chrom", "arm_simple"], as_index=False)
        .agg(start=("chromStart", "min"), end=("chromEnd", "max"))
        .rename(columns={"chrom": "chr"})
    )
    arm_df["start"] = arm_df["start"].astype(int) + 1
    arm_df["end"] = arm_df["end"].astype(int)
    arm_df["arm"] = arm_df["chr"] + arm_df["arm_simple"]
    arm_df["chr_order"] = arm_df["chr"].str.replace("chr", "", regex=False).astype(int)
    arm_df["arm_order"] = arm_df["arm_simple"].map({"p": 0, "q": 1})
    arm_df = arm_df.sort_values(["chr_order", "arm_order"]).reset_index(drop=True)
    return arm_df[["chr", "start", "end", "arm"]]


def weighted_median(values, weights):
    mask = (~np.isnan(values)) & (~np.isnan(weights)) & (weights > 0)
    if not np.any(mask):
        return np.nan
    values = values[mask]
    weights = weights[mask]
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cum_weights = np.cumsum(weights) / weights.sum()
    return float(values[np.searchsorted(cum_weights, 0.5, side="left")])


def classify_arm(arm_tcn, ploidy, threshold=THRESHOLD):
    if pd.isna(arm_tcn) or pd.isna(ploidy):
        return None
    if arm_tcn > ploidy + threshold:
        return 1
    if arm_tcn < ploidy - threshold:
        return -1
    return 0


def prepare_segments(
    combined_manifest=COMBINED_MANIFEST,
    output_dir=OUTPUT_DIR,
    segment_dir=SEGMENT_DIR,
    results_dir=RESULTS_DIR,
):
    combined_manifest_df = pd.read_csv(combined_manifest, sep="\t")
    segment_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    skipped = []
    for record in combined_manifest_df.itertuples(index=False):
        sample = str(record.pair_folder)
        source_file = Path(record.output_cna_file)
        meta, df = read_clincnv_table(source_file)
        if df.empty:
            skipped.append((sample, "empty_cna_table"))
            continue

        required = {"chr", "start", "end", "tumor_CN_change"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"{source_file} missing required columns: {sorted(missing)}")

        seg = df.loc[:, ["chr", "start", "end", "tumor_CN_change"]].copy()
        seg["chr"] = seg["chr"].map(clean_chr)
        seg["start"] = pd.to_numeric(seg["start"], errors="coerce")
        seg["end"] = pd.to_numeric(seg["end"], errors="coerce")
        seg["TCN"] = pd.to_numeric(seg["tumor_CN_change"], errors="coerce")
        seg = seg.dropna(subset=["chr", "start", "end", "TCN"]).copy()
        seg = seg[seg["end"] > seg["start"]].copy()
        seg = seg[seg["chr"].isin(AUTOSOMES)].copy()
        if seg.empty:
            skipped.append((sample, "no_valid_autosomal_segments"))
            continue

        seg["start"] = seg["start"].astype(int)
        seg["end"] = seg["end"].astype(int)
        seg["chr_order"] = seg["chr"].str.replace("chr", "", regex=False).astype(int)
        seg = seg.sort_values(["chr_order", "start", "end"]).drop(columns="chr_order")

        out_file = segment_dir / f"{sample}.aneuploidy_input.tsv"
        seg[["chr", "start", "end", "TCN"]].to_csv(out_file, sep="\t", index=False)

        ploidy = meta["ploidy"]
        rows.append(
            {
                "sample": sample,
                "source_file": str(source_file),
                "prepared_file": str(out_file.relative_to(output_dir.parent)),
                "ploidy_header": float(ploidy) if ploidy not in [None, ""] else np.nan,
                "sex": meta["sex"],
                "n_segments": int(len(seg)),
            }
        )

    manifest_df = pd.DataFrame(rows).sort_values("sample").reset_index(drop=True)
    manifest_df.to_csv(output_dir / "aneuploidy_manifest.tsv", sep="\t", index=False)
    if skipped:
        skipped_df = pd.DataFrame(skipped, columns=["sample", "reason"]).sort_values("sample")
        print("Skipped combined-run samples:")
        print(skipped_df.to_string(index=False))
    return manifest_df


def compute_scores(manifest_df, arm_df, root=ROOT):
    summary_rows = []
    arm_rows = []

    for record in manifest_df.itertuples(index=False):
        sample = str(record.sample)
        ploidy = float(record.ploidy_header) if pd.notna(record.ploidy_header) else np.nan
        seg_file = root / str(record.prepared_file)
        seg = pd.read_csv(seg_file, sep="\t")

        if seg.empty:
            summary_rows.append(
                {
                    "sample": sample,
                    "ploidy": ploidy,
                    "aneuploidy_score": np.nan,
                    "n_arms_called": np.nan,
                    "status": "no_valid_segments",
                }
            )
            continue

        if pd.isna(ploidy):
            summary_rows.append(
                {
                    "sample": sample,
                    "ploidy": np.nan,
                    "aneuploidy_score": np.nan,
                    "n_arms_called": np.nan,
                    "status": "missing_ploidy",
                }
            )
            continue

        sample_arms = []
        for arm in arm_df.itertuples(index=False):
            chr_seg = seg[seg["chr"] == arm.chr]
            if chr_seg.empty:
                continue

            overlap_start = np.maximum(chr_seg["start"].to_numpy(), arm.start)
            overlap_end = np.minimum(chr_seg["end"].to_numpy(), arm.end)
            overlap_width = overlap_end - overlap_start + 1
            overlap_mask = overlap_width > 0
            if not np.any(overlap_mask):
                continue

            arm_wmedian = weighted_median(
                chr_seg.loc[overlap_mask, "TCN"].to_numpy(dtype=float),
                overlap_width[overlap_mask].astype(float),
            )
            arm_call = classify_arm(arm_wmedian, ploidy)
            sample_arms.append(
                {
                    "sample": sample,
                    "ploidy": ploidy,
                    "arm": arm.arm,
                    "arm_wMedian": arm_wmedian,
                    "arm_call": arm_call,
                }
            )

        if not sample_arms:
            summary_rows.append(
                {
                    "sample": sample,
                    "ploidy": ploidy,
                    "aneuploidy_score": np.nan,
                    "n_arms_called": np.nan,
                    "status": "no_arm_overlaps",
                }
            )
            continue

        arm_calls_df = pd.DataFrame(sample_arms).sort_values("arm").reset_index(drop=True)
        arm_rows.extend(arm_calls_df.to_dict(orient="records"))
        summary_rows.append(
            {
                "sample": sample,
                "ploidy": ploidy,
                "aneuploidy_score": int(np.abs(arm_calls_df["arm_call"]).sum()),
                "n_arms_called": int(len(arm_calls_df)),
                "status": "ok",
            }
        )

    summary_df = pd.DataFrame(summary_rows).sort_values("sample").reset_index(drop=True)
    arm_calls_df = pd.DataFrame(arm_rows).sort_values(["sample", "arm"]).reset_index(drop=True)
    return summary_df, arm_calls_df


def write_outputs(summary_df, arm_calls_df, results_dir=RESULTS_DIR):
    summary_tsv = results_dir / "aneuploidy_score_summary.tsv"
    summary_csv = results_dir / "aneuploidy_score_summary.csv"
    arm_tsv = results_dir / "aneuploidy_arm_calls.tsv"
    arm_csv = results_dir / "aneuploidy_arm_calls.csv"

    summary_df.to_csv(summary_tsv, sep="\t", index=False, na_rep="NA")
    summary_df.to_csv(summary_csv, index=False, na_rep="NA")
    arm_calls_df.to_csv(arm_tsv, sep="\t", index=False, na_rep="NA")
    arm_calls_df.to_csv(arm_csv, index=False, na_rep="NA")


def main():
    combined_manifest_rows = len(pd.read_csv(COMBINED_MANIFEST, sep="\t"))
    arm_df = build_arm_table()
    manifest_df = prepare_segments()
    summary_df, arm_calls_df = compute_scores(manifest_df, arm_df)
    write_outputs(summary_df, arm_calls_df)

    print(f"Combined manifest rows: {combined_manifest_rows}")
    print(f"Prepared segment files: {len(manifest_df)}")
    print(f"Summary rows: {len(summary_df)}")
    print(f"Arm-call rows: {len(arm_calls_df)}")
    print("Status counts:")
    print(summary_df["status"].value_counts(dropna=False).sort_index().to_string())
    print("Outputs written:")
    print(OUTPUT_DIR / "aneuploidy_manifest.tsv")
    print(RESULTS_DIR / "aneuploidy_score_summary.tsv")
    print(RESULTS_DIR / "aneuploidy_arm_calls.tsv")


if __name__ == "__main__":
    main()
