"""Smoke tests that do not require external tools."""

from __future__ import annotations

import bamboo


def test_public_api_exports() -> None:
    import bamboo as bm

    assert bm.AlignmentFile is not None
    assert bm.CramFile is not None
    assert bm.VariantFile is not None
    assert bm.read_vcf_table is not None
    assert bm.read_bcf_table is not None
    assert bm.read_cram_columns is not None


def test_public_api_exports_bam() -> None:
    assert bamboo.__version__ == "0.1.0"
    assert bamboo.AlignmentFile is not None
    assert bamboo.AlignedSegment is not None
    assert bamboo.read_bam_table is not None
    assert bamboo.scan_bam_table is not None