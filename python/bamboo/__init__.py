"""
Bamboo
======

High-performance, modern Python library for high-throughput sequencing (HTS) data.
"""

from __future__ import annotations

__version__ = "0.1.0"

try:
    from bamboo._bamboo import (  # type: ignore[attr-defined]
        AlignedSegment,
        AlignmentFile,
        AlignmentIterator,
        CramFile,
        VariantFile,
        VariantRecord,
        __version__ as _rust_version,
        read_bam_table,
        read_bcf_table,
        read_columns,
        read_cram_columns,
        read_vcf_table,
        scan_bam_table,
        pileup_available,
        htslib_available,
        primary_backend,
    )
    from bamboo import compat
    from bamboo.arrow import to_pandas, to_polars
except ImportError:
    AlignedSegment = None  # type: ignore[assignment,misc]
    AlignmentFile = None  # type: ignore[assignment,misc]
    AlignmentIterator = None  # type: ignore[assignment,misc]
    CramFile = None  # type: ignore[assignment,misc]
    VariantFile = None  # type: ignore[assignment,misc]
    VariantRecord = None  # type: ignore[assignment,misc]
    read_bam_table = None  # type: ignore[assignment,misc]
    read_bcf_table = None  # type: ignore[assignment,misc]
    read_columns = None  # type: ignore[assignment,misc]
    read_cram_columns = None  # type: ignore[assignment,misc]
    read_vcf_table = None  # type: ignore[assignment,misc]
    scan_bam_table = None  # type: ignore[assignment,misc]
    compat = None  # type: ignore[assignment,misc]
    to_polars = None  # type: ignore[assignment,misc]
    to_pandas = None  # type: ignore[assignment,misc]
    pileup_available = None  # type: ignore[assignment,misc]
    htslib_available = None  # type: ignore[assignment,misc]
    primary_backend = None  # type: ignore[assignment,misc]
    _rust_version = None

if _rust_version is not None:
    __version__ = _rust_version

__all__ = [
    "__version__",
    "AlignedSegment",
    "AlignmentFile",
    "AlignmentIterator",
    "CramFile",
    "VariantFile",
    "VariantRecord",
    "compat",
    "read_bam_table",
    "read_bcf_table",
    "read_columns",
    "read_cram_columns",
    "read_vcf_table",
    "scan_bam_table",
    "to_polars",
    "to_pandas",
    "pileup_available",
    "htslib_available",
    "primary_backend",
]