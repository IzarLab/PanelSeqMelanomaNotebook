#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd


ROOT = Path("/mnt/myvolume/panel_seq/v4_v5_clincnv/new_clincnv_runs/combined_runs/bed_split_s100_l3_f1")

RUN_CONFIGS = {
    "ssSC_v4LI": {
        "pair_file": ROOT / "pairs_ssSC_v4LI.csv",
        "input_base": Path("/mnt/myvolume/panel_seq/v4_v5_clincnv"),
        "output_base": ROOT / "coverage_subsets" / "ssSC_v4LI",
    },
    "ssSC_v5": {
        "pair_file": ROOT / "pairs_ssSC_v5.csv",
        "input_base": Path("/mnt/myvolume/panel_seq/new_bed_analysis"),
        "output_base": ROOT / "coverage_subsets" / "ssSC_v5",
    },
}

FILES_TO_SUBSET = {
    "tumor": [
        "tumor.ontarget.wes_mapq0_bedcoverage.cov",
        "tumor.offtarget.wes_100kb_mapq10_bedcoverage.cov",
    ],
    "normal": [
        "normal.ontarget.wes_mapq0_bedcoverage.cov",
        "normal.offtarget.wes_100kb_mapq10_bedcoverage.cov",
    ],
}


def subset_coverage_file(input_path: Path, output_path: Path, keep_samples: list[str]) -> int:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with input_path.open("r", newline="") as infile:
        reader = csv.reader(infile, delimiter="\t")
        header = next(reader)

        sample_to_idx = {sample: idx for idx, sample in enumerate(header[3:], start=3)}
        missing = [sample for sample in keep_samples if sample not in sample_to_idx]
        if missing:
            raise ValueError(f"{input_path} is missing requested sample columns: {missing[:10]}")

        keep_indices = [0, 1, 2] + [sample_to_idx[sample] for sample in keep_samples]

        with output_path.open("w", newline="") as outfile:
            writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
            writer.writerow(header[:3] + keep_samples)
            row_count = 0
            for row in reader:
                writer.writerow([row[idx] for idx in keep_indices])
                row_count += 1

    return row_count


def main() -> int:
    for panel_name, config in RUN_CONFIGS.items():
        pair_df = pd.read_csv(config["pair_file"], header=None, names=["tumor", "normal"]).astype(str)
        pair_df["tumor"] = pair_df["tumor"].str.strip()
        pair_df["normal"] = pair_df["normal"].str.strip()

        tumor_samples = sorted(pair_df["tumor"].unique().tolist())
        normal_samples = sorted(pair_df["normal"].unique().tolist())

        print(
            f"{panel_name}: {len(pair_df)} pairs, "
            f"{len(tumor_samples)} tumor samples, {len(normal_samples)} normal samples"
        )

        for cohort_key, filenames in FILES_TO_SUBSET.items():
            keep_samples = tumor_samples if cohort_key == "tumor" else normal_samples
            for filename in filenames:
                input_path = config["input_base"] / filename
                output_path = config["output_base"] / filename
                row_count = subset_coverage_file(input_path, output_path, keep_samples)
                print(
                    f"  wrote {output_path} with {len(keep_samples)} samples and {row_count} regions"
                )

    print("Coverage subset generation complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
