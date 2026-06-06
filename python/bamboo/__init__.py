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
        __version__ as _rust_version,
        read_bam_table,
        scan_bam_table,
    )
    from bamboo.arrow import to_polars
except ImportError:
    AlignedSegment = None  # type: ignore[assignment,misc]
    AlignmentFile = None  # type: ignore[assignment,misc]
    AlignmentIterator = None  # type: ignore[assignment,misc]
    read_bam_table = None  # type: ignore[assignment,misc]
    scan_bam_table = None  # type: ignore[assignment,misc]
    to_polars = None  # type: ignore[assignment,misc]
    _rust_version = None

if _rust_version is not None:
    __version__ = _rust_version

__all__ = [
    "__version__",
    "AlignedSegment",
    "AlignmentFile",
    "AlignmentIterator",
    "read_bam_table",
    "scan_bam_table",
    "to_polars",
]