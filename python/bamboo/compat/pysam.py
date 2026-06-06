"""
pysam-compatible entry points for Bamboo.

Use this module when migrating pipelines that import pysam for BAM reading::

    # before
    import pysam

    # after
    from bamboo.compat import pysam

The object API mirrors common pysam patterns. For analytics workloads prefer
``bamboo.read_columns()`` or ``AlignmentFile.read_table()`` — they stay in
Rust end-to-end and avoid per-record Python overhead.
"""

from __future__ import annotations

from bamboo import AlignedSegment, AlignmentFile, read_bam_table, read_columns

# Historical pysam alias
Samfile = AlignmentFile
AlignedRead = AlignedSegment


def alignment_count(filename: str) -> int:
    """Count alignments in a BAM (``samtools view -c`` equivalent)."""
    with AlignmentFile(filename) as bam:
        return bam.count()


def read_alignment_columns(
    filename: str,
    *,
    columns: list[str] | None = None,
    region: str | None = None,
    min_mapq: int | None = None,
):
    """Fast columnar read — Bamboo's preferred path over record iteration."""
    return read_columns(
        filename,
        columns=columns,
        region=region,
        min_mapq=min_mapq,
    )


__all__ = [
    "AlignedRead",
    "AlignedSegment",
    "AlignmentFile",
    "Samfile",
    "alignment_count",
    "read_alignment_columns",
    "read_bam_table",
    "read_columns",
]