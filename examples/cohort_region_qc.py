#!/usr/bin/env python3
"""
Cohort region QC — Bamboo's killer workflow.

Indexed region → columnar Arrow scan → Polars aggregation.

This is the path Bamboo is built for: skip per-record Python objects, stay in
Rust until you have a dataframe ready for filtering and summary stats.

Usage:
    python examples/cohort_region_qc.py aligned.bam --region chr1:1000000-5000000
    python examples/cohort_region_qc.py sample.cram --region chr1:1-5000000 \\
        --reference ref.fasta

Requires: bamboo, pyarrow. Optional: polars (falls back to printing Arrow table).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import bamboo as bm


def _load_region_table(
    path: Path,
    *,
    region: str,
    min_mapq: int,
    reference: Path | None,
) -> object:
    columns = ["qname", "rname", "pos", "mapq", "flag", "cigar"]
    suffix = path.suffix.lower()

    if suffix == ".cram":
        if reference is None:
            raise SystemExit("CRAM requires --reference FASTA")
        return bm.read_cram_columns(
            str(path),
            columns=columns,
            region=region,
            min_mapq=min_mapq,
            reference_filename=str(reference),
        )

    return bm.read_columns(
        str(path),
        columns=columns,
        region=region,
        min_mapq=min_mapq,
    )


def _summarize(table) -> dict[str, object]:
    try:
        import polars as pl
    except ImportError:
        return {"arrow_rows": table.num_rows, "arrow_columns": table.column_names}

    df = bm.to_polars(table)
    return {
        "reads": df.height,
        "mean_mapq": df["mapq"].mean(),
        "median_mapq": df["mapq"].median(),
        "unique_references": df["rname"].n_unique(),
        "top_positions": (
            df.group_by("pos")
            .len()
            .sort("len", descending=True)
            .head(5)
            .to_dicts()
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Region QC via Bamboo columnar scan (pysam-compatible inputs)",
    )
    parser.add_argument("alignment", type=Path, help="BAM or CRAM path (local or cloud URI)")
    parser.add_argument(
        "--region",
        default="chr1:1000000-5000000",
        help="samtools-style region (default: chr1:1000000-5000000)",
    )
    parser.add_argument("--min-mapq", type=int, default=0, help="Minimum mapping quality")
    parser.add_argument(
        "--reference",
        type=Path,
        default=None,
        help="Reference FASTA for CRAM",
    )
    args = parser.parse_args()

    table = _load_region_table(
        args.alignment,
        region=args.region,
        min_mapq=args.min_mapq,
        reference=args.reference,
    )

    summary = _summarize(table)
    print(f"Region: {args.region}")
    print(f"File:   {args.alignment}")
    for key, value in summary.items():
        print(f"{key}: {value}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())