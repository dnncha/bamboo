"""
Bamboo
======

High-performance, modern Python library for high-throughput sequencing (HTS) data.

Bamboo is designed as the successor to pysam for the 2020s/2030s Python data science world.

Key features (planned / in progress):
- Blazing fast Rust core
- Native Arrow / Polars / pandas integration
- First-class cloud object store support (S3, GCS, etc.)
- Clean, Pythonic, well-typed API
- Streaming and zero-copy friendly

Example usage (aspirational):

    import bamboo as bm
    import polars as pl

    with bm.AlignmentFile("data.bam") as bam:
        for read in bam.fetch("chr1", 1_000_000, 1_001_000):
            ...

        # Or get a proper DataFrame
        df = bam.to_polars()

Current status: Very early skeleton. The real work is just beginning.
"""

from __future__ import annotations

__version__ = "0.1.0"

# The Rust extension module will be built as bamboo._bamboo
try:
    from bamboo._bamboo import hello, __version__ as _rust_version  # type: ignore
except ImportError:
    def hello() -> str:
        return "Bamboo Rust extension not built yet. Run `maturin develop` or `pip install -e .`"

    _rust_version = None

__all__ = ["__version__", "hello"]
