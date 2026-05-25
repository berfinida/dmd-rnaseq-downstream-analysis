#!/usr/bin/env python3
from pathlib import Path

import pandas as pd


def main() -> None:
    project_root = Path(__file__).resolve().parents[1]
    summary_path = project_root / "results" / "matrix" / "dmd_vs_wt_summary.tsv"
    mapping_path = project_root / "annotation" / "transcript_to_gene_symbol.tsv"
    output_path = project_root / "results" / "matrix" / "dmd_vs_wt_summary_annotated.tsv"

    summary = pd.read_csv(summary_path, sep="\t")
    if "Name" not in summary.columns:
        raise ValueError("Missing 'Name' column in dmd_vs_wt_summary.tsv")

    if not mapping_path.exists():
        print(
            "No mapping file found at annotation/transcript_to_gene_symbol.tsv. "
            "Skipping annotation."
        )
        return

    mapping = pd.read_csv(mapping_path, sep="\t")
    required = {"transcript_id", "gene_symbol"}
    missing = required - set(mapping.columns)
    if missing:
        raise ValueError(
            "Mapping file is missing required columns: "
            + ", ".join(sorted(missing))
        )

    mapping = (
        mapping[["transcript_id", "gene_symbol"]]
        .drop_duplicates(subset=["transcript_id"])
        .copy()
    )
    merged = summary.merge(
        mapping,
        how="left",
        left_on="Name",
        right_on="transcript_id",
    )
    merged = merged.drop(columns=["transcript_id"])
    merged.to_csv(output_path, sep="\t", index=False)

    hit_count = int(merged["gene_symbol"].notna().sum())
    print(f"Saved: {output_path}")
    print(f"Annotated rows with gene_symbol: {hit_count}/{len(merged)}")


if __name__ == "__main__":
    main()
