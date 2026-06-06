#!/usr/bin/env python3
"""Minimal Bamboo BAM reader example."""

from __future__ import annotations

import argparse
import sys

import bamboo as bm


def main() -> int:
    parser = argparse.ArgumentParser(description="Read a BAM with Bamboo")
    parser.add_argument("bam", help="Path to a coordinate-sorted BAM file")
    parser.add_argument(
        "--region",
        help="samtools-style region, e.g. chr1:1000-2000",
    )
    parser.add_argument(
        "--to-arrow",
        action="store_true",
        help="Print a PyArrow table instead of iterating records",
    )
    args = parser.parse_args()

    with bm.AlignmentFile(args.bam) as bam:
        print(f"References: {list(zip(bam.references(), bam.reference_lengths))}")

        if args.to_arrow:
            table = bam.to_arrow(region=args.region)
            print(table)
            return 0

        if args.region:
            reads = bam.fetch(region=args.region)
        else:
            reads = bam

        for read in reads:
            print(
                read.query_name,
                read.reference_name,
                read.reference_start,
                read.mapping_quality,
                read.cigarstring,
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())