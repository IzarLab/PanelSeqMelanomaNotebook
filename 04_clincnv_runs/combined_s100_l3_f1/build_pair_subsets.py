#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path("/mnt/myvolume/panel_seq/v4_v5_clincnv/new_clincnv_runs/combined_runs/bed_split_s100_l3_f1")
ANNOTATION_XLSX = Path("/mnt/myvolume/panel_seq/v4_v5_clincnv/MEL_Annotation.xlsx")
ANNOTATION_SHEET = "MFT_SINNBERG_QHPRG_MeasTabelle9"
PAIR_CSV = Path("/mnt/myvolume/panel_seq/new_bed_analysis/pairs_df_filtered.csv")

EXPECTED_COUNTS = {
    "ssSC_v4LI": 232,
    "ssSC_v5": 108,
}

OUTPUT_FILES = {
    "ssSC_v4LI": ROOT / "pairs_ssSC_v4LI.csv",
    "ssSC_v5": ROOT / "pairs_ssSC_v5.csv",
}


def main() -> int:
    annotation_df = pd.read_excel(ANNOTATION_XLSX, sheet_name=ANNOTATION_SHEET)
    sample_to_panel = (
        annotation_df[["code", "system_name_short"]]
        .dropna(subset=["code", "system_name_short"])
        .drop_duplicates(subset=["code"], keep="last")
        .assign(
            code=lambda df: df["code"].astype(str).str.strip(),
            system_name_short=lambda df: df["system_name_short"].astype(str).str.strip(),
        )
        .set_index("code")["system_name_short"]
        .to_dict()
    )

    pair_df = pd.read_csv(PAIR_CSV, header=None, names=["tumor", "normal"]).astype(str)
    pair_df["tumor"] = pair_df["tumor"].str.strip()
    pair_df["normal"] = pair_df["normal"].str.strip()
    pair_df["system_name_short"] = pair_df["tumor"].map(sample_to_panel)

    if pair_df["system_name_short"].isna().any():
        missing = pair_df.loc[pair_df["system_name_short"].isna(), "tumor"].tolist()
        raise ValueError(f"Missing panel assignments for tumor samples: {missing}")

    ROOT.mkdir(parents=True, exist_ok=True)

    for panel_key, out_path in OUTPUT_FILES.items():
        subset_df = pair_df.loc[pair_df["system_name_short"].eq(panel_key), ["tumor", "normal"]].copy()
        expected = EXPECTED_COUNTS[panel_key]
        if len(subset_df) != expected:
            raise ValueError(f"{panel_key} subset has {len(subset_df)} pairs, expected {expected}")
        subset_df.to_csv(out_path, header=False, index=False)
        print(f"Wrote {len(subset_df)} pairs to {out_path}")

    total_pairs = sum(EXPECTED_COUNTS.values())
    if len(pair_df) != total_pairs:
        raise ValueError(f"Expected {total_pairs} total pairs, found {len(pair_df)}")

    print("Pair subset generation complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
