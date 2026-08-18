#!/usr/bin/env python3
"""Build the merged OncoKB + ClinCNV summary table used by Figures 2 and 3.

This script builds `all_oncokb_genes_annotated.txt` from two staged inputs:

1. Per-sample OncoKB TSV files produced from filtered annotated somatic VCFs.
2. Per-sample annotated ClinCNV CNA tables.

The output columns are the ones consumed by the figure notebooks, including the
derived `Biallelic_Loss_Genes`, `Functional_Oncogenic_Complete_LOF`, and
`Functional_Oncogenic_GOF` fields.
"""

from __future__ import annotations

import argparse
import glob
import io
import re
from pathlib import Path

import pandas as pd


QHPRG_RE = re.compile(r"(QHPRG[0-9A-Z]+)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--oncokb-dir",
        type=Path,
        default=Path("/mnt/myvolume/panel_seq/new_bed_analysis/oncokb/oncokb_results_min_no_pass_vaf005"),
        help="Directory containing per-sample `.oncokb.tsv` files.",
    )
    parser.add_argument(
        "--clincnv-root",
        type=Path,
        default=Path("/mnt/myvolume/panel_seq/v4_v5_clincnv/new_clincnv_runs/combined_runs/combined_bed_split_s100_l3_f1"),
        help="Root directory containing pair-level annotated ClinCNV CNA tables.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("all_oncokb_genes_annotated.txt"),
        help="Output TSV path.",
    )
    return parser.parse_args()


def extract_qhprg(text: str) -> str:
    match = QHPRG_RE.search(str(text))
    return match.group(1) if match else ""


def split_genes(cell: str) -> set[str]:
    if not isinstance(cell, str) or not cell.strip():
        return set()
    return {token.strip() for token in cell.split(",") if token.strip()}


def only_genes(label_str: str) -> str:
    if not isinstance(label_str, str) or not label_str.strip():
        return ""
    out: list[str] = []
    seen: set[str] = set()
    for item in (x.strip() for x in label_str.split(",") if x.strip()):
        gene = item.split()[0]
        if gene and gene not in seen:
            seen.add(gene)
            out.append(gene)
    return ", ".join(out)


def only_genes_from_labels(cell: str) -> set[str]:
    if not isinstance(cell, str) or not cell.strip():
        return set()
    genes: set[str] = set()
    for item in (x.strip() for x in cell.split(",") if x.strip()):
        parts = item.split()
        if parts:
            genes.add(parts[0])
    return genes


def count_items(value: str) -> int:
    if not isinstance(value, str) or not value.strip():
        return 0
    return sum(1 for token in value.split(",") if token.strip())


def get_col(df: pd.DataFrame, names: list[str]) -> pd.Series:
    for name in names:
        if name in df.columns:
            return df[name]
    lower_map = {col.lower(): col for col in df.columns}
    for name in names:
        if name.lower() in lower_map:
            return df[lower_map[name.lower()]]
    return pd.Series([""] * len(df), index=df.index)


def normalize_cna_type(series: pd.Series) -> pd.Series:
    mapping = {
        "amplification": "amplification",
        "amp": "amplification",
        "gain": "gain",
        "deletion": "deletion",
        "del": "deletion",
        "loss": "loss",
    }
    return series.astype(str).str.strip().str.lower().map(lambda x: mapping.get(x, x))


def is_gof(effect: str) -> bool:
    return "gain-of-function" in str(effect).strip().lower()


def is_lof(effect: str) -> bool:
    return "loss-of-function" in str(effect).strip().lower()


def find_oncokb_files(root: Path) -> list[Path]:
    return sorted(root.glob("*.oncokb.tsv"))


def find_clincnv_files(root: Path) -> list[Path]:
    patterns = [
        "**/Annotated_ClinCNV_*.txt",
        "**/*_Annotated_ClinCNV.txt",
        "**/Annotated_CNA_CNAs_*.txt",
    ]
    hits: set[str] = set()
    for pattern in patterns:
        hits.update(glob.glob(str(root / pattern), recursive=True))
    return [Path(path) for path in sorted(hits)]


def load_clincnv_table(path: Path) -> pd.DataFrame:
    """Read a ClinCNV table while skipping injected provenance headers."""
    with path.open() as handle:
        lines = handle.readlines()

    start_idx = 0
    for idx, line in enumerate(lines):
        if line.startswith("#chr\t"):
            lines[idx] = line[1:]
            start_idx = idx
            break
        if line.startswith("chr\t"):
            start_idx = idx
            break

    text = "".join(lines[start_idx:])
    return pd.read_csv(io.StringIO(text), sep="\t", dtype=str, engine="python").fillna("")


def summarize_oncokb_sample(path: Path) -> dict[str, object]:
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    sample = extract_qhprg(path.name)
    if not sample and "Tumor_Sample_Barcode" in df.columns and not df.empty:
        sample = extract_qhprg(df["Tumor_Sample_Barcode"].iloc[0])
    if not sample:
        raise ValueError(f"Could not extract QHPRG code from {path}")

    genes = sorted(
        gene
        for gene in df.get("ONCOKB_HUGO_SYMBOL", pd.Series([], dtype=str)).astype(str).unique().tolist()
        if gene.strip()
    )

    gof_mask = df.get("MUTATION_EFFECT", pd.Series([], dtype=str)).str.contains("Gain-of-function", case=False, na=False)
    lof_mask = df.get("MUTATION_EFFECT", pd.Series([], dtype=str)).str.contains("Loss-of-function", case=False, na=False)
    oncogenic_mask = df.get("ONCOGENIC", pd.Series([], dtype=str)).str.contains("Oncogenic", case=False, na=False)

    df["variant_label"] = (
        df.get("ONCOKB_HUGO_SYMBOL", pd.Series([], dtype=str)).astype(str).str.strip()
        + " "
        + df.get("ONCOKB_PROTEIN_CHANGE", pd.Series([], dtype=str)).astype(str).str.strip()
    ).str.strip()
    df["variant_label"] = df["variant_label"].str.replace(r"\s+", " ", regex=True)

    gof_variants = df.loc[gof_mask & oncogenic_mask, "variant_label"].replace("", pd.NA).dropna().unique().tolist()
    lof_one = df.loc[lof_mask & oncogenic_mask, "variant_label"].replace("", pd.NA).dropna().unique().tolist()

    lof_counts = (
        df.loc[lof_mask & oncogenic_mask, "ONCOKB_HUGO_SYMBOL"]
        .astype(str)
        .str.strip()
        .replace("", pd.NA)
        .dropna()
        .value_counts()
    )
    lof_both = lof_counts[lof_counts > 1].index.tolist()

    protein_change = df.get("ONCOKB_PROTEIN_CHANGE", pd.Series([], dtype=str)).astype(str)
    gene_col = df.get("ONCOKB_HUGO_SYMBOL", pd.Series([], dtype=str)).astype(str)

    return {
        "Sample": sample,
        "GOF_Oncogenic_Variants": ", ".join(gof_variants) if gof_variants else "",
        "LOF_Oncogenic_One_Allele": ", ".join(lof_one) if lof_one else "",
        "LOF_Oncogenic_Both_Alleles": ", ".join(lof_both) if lof_both else "",
        "Has_V600_variant": bool(protein_change.str.contains("V600", case=False, na=False).any()),
        "BRAF_GOF": bool(((gene_col == "BRAF") & gof_mask).any()),
        "NRAS_GOF": bool(((gene_col == "NRAS") & gof_mask).any()),
        "NF1_LOF": bool(((gene_col == "NF1") & lof_mask).any()),
        "TERT_GOF": bool(((gene_col == "TERT") & gof_mask).any()),
        "KIT_GOF": bool(((gene_col == "KIT") & gof_mask).any()),
        "All_Variant_Genes": ", ".join(genes),
    }


def build_oncokb_variant_summary(root: Path) -> pd.DataFrame:
    rows = [summarize_oncokb_sample(path) for path in find_oncokb_files(root)]
    if not rows:
        raise FileNotFoundError(f"No `.oncokb.tsv` files found in {root}")

    summary = pd.DataFrame(rows)
    summary["Variant_Richness"] = (
        summary["GOF_Oncogenic_Variants"].map(count_items)
        + summary["LOF_Oncogenic_One_Allele"].map(count_items)
        + summary["LOF_Oncogenic_Both_Alleles"].map(count_items)
    )
    summary = (
        summary.sort_values(["Sample", "Variant_Richness"], ascending=[True, False])
        .drop_duplicates(subset="Sample", keep="first")
        .sort_values("Sample")
        .reset_index(drop=True)
    )
    summary["GOF_Genes"] = summary["GOF_Oncogenic_Variants"].map(only_genes)
    summary["LOF_OneAllele_Genes"] = summary["LOF_Oncogenic_One_Allele"].map(only_genes)
    summary["LOF_Biallelic_Genes"] = summary["LOF_Oncogenic_Both_Alleles"].map(only_genes)
    return summary


def summarize_clincnv_sample(path: Path) -> dict[str, object]:
    df = load_clincnv_table(path)
    sample = extract_qhprg(path.name)
    if not sample:
        raise ValueError(f"Could not extract QHPRG code from {path}")

    hugo = get_col(df, ["Hugo_Symbol", "HUGO_SYMBOL"]).astype(str).str.strip()
    alt = normalize_cna_type(get_col(df, ["Copy_Number_Alteration", "ALTERATION"]))
    onco = get_col(df, ["ONCOGENIC"]).astype(str).str.strip().str.lower()
    effect = get_col(df, ["MUTATION_EFFECT"]).astype(str)

    keep = (hugo != "") & alt.isin(["amplification", "gain", "deletion", "loss"])
    work = pd.DataFrame(
        {
            "HUGO": hugo[keep],
            "ALT": alt[keep],
            "ONC": onco[keep],
            "EFF": effect[keep],
        }
    )
    if work.empty:
        return {
            "Sample": sample,
            "Oncogenic_GOF_AmpGain": "",
            "Oncogenic_GOF_AmpGain_Genes": "",
            "Oncogenic_LOF_Loss": "",
            "Oncogenic_LOF_Loss_Genes": "",
            "Oncogenic_LOF_Deletions": "",
            "Oncogenic_LOF_Deletions_Genes": "",
            "All_GainAmp_Genes": "",
            "All_Loss_Genes": "",
            "All_Deletion_Genes": "",
            "Total_CNA_Genes": 0,
        }

    onc_mask = work["ONC"].isin(["oncogenic", "likely oncogenic"])
    gof_mask = work["EFF"].map(is_gof)
    lof_mask = work["EFF"].map(is_lof)
    amp_mask = work["ALT"].eq("amplification")
    gain_mask = work["ALT"].eq("gain")
    del_mask = work["ALT"].eq("deletion")
    loss_mask = work["ALT"].eq("loss")

    work["label"] = work["HUGO"] + " " + work["ALT"].str.title()

    gof_ampgain = work.loc[onc_mask & gof_mask & (amp_mask | gain_mask), "label"].dropna().unique().tolist()
    lof_loss = work.loc[onc_mask & lof_mask & loss_mask, "label"].dropna().unique().tolist()
    lof_del = work.loc[onc_mask & lof_mask & del_mask, "label"].dropna().unique().tolist()

    all_ampgain = work.loc[amp_mask | gain_mask, "label"].dropna().unique().tolist()
    all_loss = work.loc[loss_mask, "label"].dropna().unique().tolist()
    all_del = work.loc[del_mask, "label"].dropna().unique().tolist()

    return {
        "Sample": sample,
        "Oncogenic_GOF_AmpGain": ", ".join(gof_ampgain) if gof_ampgain else "",
        "Oncogenic_GOF_AmpGain_Genes": only_genes(", ".join(gof_ampgain)) if gof_ampgain else "",
        "Oncogenic_LOF_Loss": ", ".join(lof_loss) if lof_loss else "",
        "Oncogenic_LOF_Loss_Genes": only_genes(", ".join(lof_loss)) if lof_loss else "",
        "Oncogenic_LOF_Deletions": ", ".join(lof_del) if lof_del else "",
        "Oncogenic_LOF_Deletions_Genes": only_genes(", ".join(lof_del)) if lof_del else "",
        "All_GainAmp_Genes": only_genes(", ".join(all_ampgain)) if all_ampgain else "",
        "All_Loss_Genes": only_genes(", ".join(all_loss)) if all_loss else "",
        "All_Deletion_Genes": only_genes(", ".join(all_del)) if all_del else "",
        "Total_CNA_Genes": sum(
            count_items(value)
            for value in [
                only_genes(", ".join(gof_ampgain)) if gof_ampgain else "",
                only_genes(", ".join(lof_loss)) if lof_loss else "",
                only_genes(", ".join(lof_del)) if lof_del else "",
            ]
        ),
    }


def build_clincnv_summary(root: Path) -> pd.DataFrame:
    rows = [summarize_clincnv_sample(path) for path in find_clincnv_files(root)]
    if not rows:
        raise FileNotFoundError(f"No annotated ClinCNV files found under {root}")
    summary = pd.DataFrame(rows)
    summary = (
        summary.sort_values(["Sample", "Total_CNA_Genes"], ascending=[True, False])
        .drop_duplicates(subset="Sample", keep="first")
        .sort_values("Sample")
        .reset_index(drop=True)
    )
    return summary


def compute_biallelic_loss(row: pd.Series) -> str:
    all_loss = split_genes(row.get("All_Loss_Genes", ""))
    all_deletion = split_genes(row.get("All_Deletion_Genes", ""))
    all_variant = split_genes(row.get("LOF_OneAllele_Genes", "")) | split_genes(row.get("LOF_Biallelic_Genes", ""))
    gof_variants = split_genes(row.get("GOF_Oncogenic_Variants", ""))
    cond1 = (all_loss & all_variant) - gof_variants
    cond2 = all_deletion
    final = cond1 | cond2
    return ",".join(sorted(final)) if final else ""


def add_functional_columns(merged: pd.DataFrame) -> pd.DataFrame:
    merged = merged.copy()
    merged["Biallelic_Loss_Genes"] = merged.apply(compute_biallelic_loss, axis=1)

    one_series = merged["LOF_Oncogenic_One_Allele"].fillna("").map(only_genes_from_labels)
    both_series = merged["LOF_Oncogenic_Both_Alleles"].fillna("").map(only_genes_from_labels)

    cohort_lof_genes: set[str] = set()
    if len(one_series):
        cohort_lof_genes = cohort_lof_genes.union(*one_series.tolist())
    if len(both_series):
        cohort_lof_genes = cohort_lof_genes.union(*both_series.tolist())

    def functional_complete_lof_row(row: pd.Series) -> str:
        genes_one = only_genes_from_labels(row.get("LOF_Oncogenic_One_Allele", ""))
        genes_both = only_genes_from_labels(row.get("LOF_Oncogenic_Both_Alleles", ""))
        loss_genes = split_genes(row.get("All_Loss_Genes", ""))
        del_genes = split_genes(row.get("All_Deletion_Genes", ""))
        del_oncogenes = split_genes(row.get("Oncogenic_LOF_Deletions_Genes", ""))
        rule1 = genes_one & loss_genes
        rule2 = genes_both
        rule3 = del_genes & cohort_lof_genes
        return ", ".join(sorted(rule1 | rule2 | rule3 | del_oncogenes))

    if "GOF_Genes" in merged.columns:
        gof_genes_series = merged["GOF_Genes"].fillna("").map(split_genes)
    else:
        gof_genes_series = merged.get("GOF_Oncogenic_Variants", pd.Series([], dtype=str)).fillna("").map(only_genes_from_labels)

    cohort_gof_genes: set[str] = set()
    if len(gof_genes_series):
        cohort_gof_genes = cohort_gof_genes.union(*gof_genes_series.tolist())

    def functional_oncogenic_gof_row(row: pd.Series) -> str:
        row_gof_genes = split_genes(row.get("GOF_Genes", "")) if "GOF_Genes" in merged.columns else only_genes_from_labels(row.get("GOF_Oncogenic_Variants", ""))
        gainamp_genes = split_genes(row.get("All_GainAmp_Genes", ""))
        rule2 = cohort_gof_genes & gainamp_genes
        rule3 = split_genes(row.get("Oncogenic_GOF_AmpGain_Genes", ""))
        final = row_gof_genes | rule2 | rule3
        return ", ".join(sorted(final))

    merged["Functional_Oncogenic_Complete_LOF"] = merged.apply(functional_complete_lof_row, axis=1)
    merged["Functional_Oncogenic_GOF"] = merged.apply(functional_oncogenic_gof_row, axis=1)
    return merged


def main() -> None:
    args = parse_args()
    variant_summary = build_oncokb_variant_summary(args.oncokb_dir)
    cnv_summary = build_clincnv_summary(args.clincnv_root)

    merged = variant_summary.merge(cnv_summary, on="Sample", how="outer")
    for col in merged.columns:
        if merged[col].dtype == object:
            merged[col] = merged[col].fillna("")
    merged = add_functional_columns(merged)
    merged = merged.sort_values("Sample").reset_index(drop=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(merged)} rows to {args.output}")


if __name__ == "__main__":
    main()
