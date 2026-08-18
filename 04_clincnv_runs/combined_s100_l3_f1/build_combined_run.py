#!/usr/bin/env python3
from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd


SPLIT_ROOT = Path("/mnt/myvolume/panel_seq/v4_v5_clincnv/new_clincnv_runs/combined_runs/bed_split_s100_l3_f1")
OUTPUT_ROOT = Path("/mnt/myvolume/panel_seq/v4_v5_clincnv/new_clincnv_runs/combined_runs/combined_bed_split_s100_l3_f1")
ANNOTATION_XLSX = Path("/mnt/myvolume/panel_seq/v4_v5_clincnv/MEL_Annotation.xlsx")
ANNOTATION_SHEET = "MFT_SINNBERG_QHPRG_MeasTabelle9"
PAIR_CSV = Path("/mnt/myvolume/panel_seq/new_bed_analysis/pairs_df_filtered.csv")

PANEL_CONFIG = {
    "ssSC_v4LI": {
        "system_name": "Sure Select XT LI Somatic Cancer Panel v4",
        "bed_assignment": "/mnt/myvolume/panel_seq/v4_v5_clincnv/ssSC_v4.gc.genes.bedcoverage.sorted.bed",
        "source_run": SPLIT_ROOT / "full_wes_ssSC_v4LI_s100_l3_f1_baf_bedcoverage_wes",
    },
    "ssSC_v5": {
        "system_name": "Sure Select Somatic Cancer Panel v5",
        "bed_assignment": "/mnt/myvolume/panel_seq/new_bed_analysis/ssSC_v5.gc.genes.bed",
        "source_run": SPLIT_ROOT / "full_wes_ssSC_v5_s100_l3_f1_baf_bedcoverage_wes",
    },
}


def read_lines(path: Path) -> list[str]:
    return path.read_text(encoding="utf-8", errors="ignore").splitlines(keepends=True)


def write_lines(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines), encoding="utf-8")


def add_provenance(lines: list[str], provenance_lines: list[str]) -> list[str]:
    insert_at = 0
    while insert_at < len(lines) and lines[insert_at].startswith("##"):
        insert_at += 1
    return lines[:insert_at] + provenance_lines + lines[insert_at:]


def copy_if_present(source: Path, destination: Path) -> bool:
    if not source.exists():
        return False
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    return True


def load_sample_panel_table() -> dict[str, dict[str, str]]:
    annotation_df = pd.read_excel(ANNOTATION_XLSX, sheet_name=ANNOTATION_SHEET)
    annotation_df = (
        annotation_df[["code", "system_name_short", "system_name"]]
        .dropna(subset=["code", "system_name_short"])
        .drop_duplicates(subset=["code"], keep="last")
        .assign(
            code=lambda df: df["code"].astype(str).str.strip(),
            system_name_short=lambda df: df["system_name_short"].astype(str).str.strip(),
            system_name=lambda df: df["system_name"].fillna("").astype(str).str.strip(),
        )
    )

    sample_to_panel: dict[str, dict[str, str]] = {}
    for row in annotation_df.itertuples(index=False):
        sample_to_panel[row.code] = {
            "system_name_short": row.system_name_short,
            "system_name": row.system_name or PANEL_CONFIG[row.system_name_short]["system_name"],
        }
    return sample_to_panel


def load_pairs() -> list[tuple[str, str]]:
    pair_df = pd.read_csv(PAIR_CSV, header=None, names=["tumor", "normal"]).astype(str)
    pair_df["tumor"] = pair_df["tumor"].str.strip()
    pair_df["normal"] = pair_df["normal"].str.strip()
    return list(pair_df.itertuples(index=False, name=None))


def main() -> int:
    somatic_root = OUTPUT_ROOT / "somatic"
    manifest_path = OUTPUT_ROOT / "combined_file_manifest.tsv"
    readme_path = OUTPUT_ROOT / "README.txt"

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    somatic_root.mkdir(parents=True, exist_ok=True)

    sample_to_panel = load_sample_panel_table()
    pairs = load_pairs()
    manifest_rows: list[dict[str, str]] = []
    selection_counts = {key: 0 for key in PANEL_CONFIG}
    missing_pairs = 0

    for tumor_sample, normal_sample in pairs:
        pair_folder = f"{tumor_sample}-{normal_sample}"
        panel_info = sample_to_panel.get(tumor_sample)
        if panel_info is None:
            raise ValueError(f"Missing panel annotation for tumor sample {tumor_sample}")

        panel_key = panel_info["system_name_short"]
        if panel_key not in PANEL_CONFIG:
            raise ValueError(f"Unsupported panel annotation {panel_key!r} for {tumor_sample}")

        config = PANEL_CONFIG[panel_key]
        source_run = config["source_run"]
        source_pair_dir = source_run / "somatic" / pair_folder
        source_cna = source_pair_dir / f"CNAs_{pair_folder}.txt"
        annotated_candidates = [
            source_pair_dir / f"Annotated_ClinCNV_{pair_folder}.txt",
            source_pair_dir / f"Annotated_CNA_CNAs_{pair_folder}.txt",
        ]
        source_annotated = next((path for path in annotated_candidates if path.exists()), None)

        output_pair_dir = somatic_root / pair_folder
        output_cna = output_pair_dir / f"CNAs_{pair_folder}.txt"
        output_annotated = output_pair_dir / f"Annotated_ClinCNV_{pair_folder}.txt"
        cov_seg_source = source_pair_dir / f"{pair_folder}_cov.seg"
        cnvs_seg_source = source_pair_dir / f"{pair_folder}_cnvs.seg"
        cov_seg_output = output_pair_dir / f"{pair_folder}_cov.seg"
        cnvs_seg_output = output_pair_dir / f"{pair_folder}_cnvs.seg"

        row = {
            "pair_folder": pair_folder,
            "tumor_sample": tumor_sample,
            "normal_sample": normal_sample,
            "system_name_short": panel_key,
            "system_name": panel_info["system_name"],
            "bed_assignment": config["bed_assignment"],
            "source_split_root": str(SPLIT_ROOT),
            "source_run": str(source_run),
            "source_pair_dir": str(source_pair_dir),
            "source_cna_file": str(source_cna),
            "source_annotated_file": str(source_annotated) if source_annotated else "",
            "output_pair_dir": str(output_pair_dir),
            "output_cna_file": str(output_cna),
            "output_annotated_file": str(output_annotated),
            "cov_seg_source_file": str(cov_seg_source),
            "cov_seg_written": "False",
            "cnvs_seg_source_file": str(cnvs_seg_source),
            "cnvs_seg_written": "False",
        }

        if not source_cna.exists():
            missing_pairs += 1
            row["status"] = "missing_cna_source"
            manifest_rows.append(row)
            continue

        if source_annotated is None:
            missing_pairs += 1
            row["status"] = "missing_annotated_source"
            manifest_rows.append(row)
            continue

        selection_counts[panel_key] += 1
        provenance = [
            "##combined_run_name=combined_bed_split_s100_l3_f1\n",
            f"##combined_assignment_basis={ANNOTATION_XLSX} sheet={ANNOTATION_SHEET}\n",
            f"##combined_pair_folder={pair_folder}\n",
            f"##combined_tumor_sample={tumor_sample}\n",
            f"##combined_system_name_short={panel_key}\n",
            f"##combined_system_name={panel_info['system_name']}\n",
            f"##combined_bed_assignment={config['bed_assignment']}\n",
            f"##combined_source_run={source_run}\n",
            f"##combined_source_pair_dir={source_pair_dir}\n",
        ]

        cna_provenance = provenance + [f"##combined_source_file={source_cna}\n"]
        write_lines(output_cna, add_provenance(read_lines(source_cna), cna_provenance))

        annotated_provenance = provenance + [
            f"##combined_source_file={source_annotated}\n",
            "##combined_annotation_mode=copied_annotated_source\n",
        ]
        write_lines(output_annotated, add_provenance(read_lines(source_annotated), annotated_provenance))

        row["cov_seg_written"] = str(copy_if_present(cov_seg_source, cov_seg_output))
        row["cnvs_seg_written"] = str(copy_if_present(cnvs_seg_source, cnvs_seg_output))
        row["status"] = "written"
        manifest_rows.append(row)

    manifest_df = pd.DataFrame(manifest_rows).sort_values(["system_name_short", "pair_folder"]).reset_index(drop=True)
    manifest_df.to_csv(manifest_path, sep="\t", index=False)

    written_pairs = int((manifest_df["status"] == "written").sum()) if not manifest_df.empty else 0
    readme = "\n".join(
        [
            "combined_bed_split_s100_l3_f1",
            "",
            "This directory contains the combined ClinCNV output assembled from the panel-specific runs.",
            "Each written pair folder under somatic/ contains:",
            "- CNAs_<pair>.txt",
            "- Annotated_ClinCNV_<pair>.txt",
            "- <pair>_cov.seg",
            "- <pair>_cnvs.seg",
            "",
            "Source runs:",
            f"- {PANEL_CONFIG['ssSC_v4LI']['source_run']}",
            f"- {PANEL_CONFIG['ssSC_v5']['source_run']}",
            "",
            "Settings:",
            "- scoreS 100",
            "- lengthS 3",
            "- filterStep 1",
            "",
            "Counts:",
            f"- total_pairs_expected={len(pairs)}",
            f"- total_pairs_written={written_pairs}",
            f"- total_pairs_missing={missing_pairs}",
            f"- ssSC_v4LI={selection_counts['ssSC_v4LI']}",
            f"- ssSC_v5={selection_counts['ssSC_v5']}",
            "",
            f"Manifest: {manifest_path}",
            "",
        ]
    )
    readme_path.write_text(readme, encoding="utf-8")

    print(f"Wrote combined run to {OUTPUT_ROOT}")
    print(f"Written pairs: {written_pairs}")
    print(f"Missing pairs: {missing_pairs}")
    print(f"Selection counts: {selection_counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
