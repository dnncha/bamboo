"""Smoke tests that do not require external tools."""

from __future__ import annotations

import bamboo


def test_public_api_exports() -> None:
    assert bamboo.__version__ == "0.1.0"
    assert bamboo.AlignmentFile is not None
    assert bamboo.AlignedSegment is not None
    assert bamboo.read_bam_table is not None
    assert bamboo.scan_bam_table is not None