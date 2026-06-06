"""Roadmap phase gates — lightweight checks that each phase's deliverables exist."""

from __future__ import annotations

from pathlib import Path

import pytest

bamboo = pytest.importorskip("bamboo")


def test_phase1_bam_read_write_api() -> None:
    assert hasattr(bamboo, "AlignmentFile")
    assert hasattr(bamboo, "read_columns")
    assert hasattr(bamboo, "read_bam_table")


def test_phase2_indexed_region_columnar(tiny_bam_with_index: Path) -> None:
    pa = pytest.importorskip("pyarrow")

    table = bamboo.read_columns(
        str(tiny_bam_with_index),
        columns=["qname", "pos"],
        region="chr1:100-101",
    )
    assert isinstance(table, pa.Table)
    assert table.num_rows >= 1


def test_phase2_packaging_metadata() -> None:
    pyproject = Path(__file__).resolve().parents[1] / "pyproject.toml"
    text = pyproject.read_text()
    assert 'name = "bamboo"' in text
    assert "pyarrow" in text


def test_phase3_htslib_stub_present() -> None:
    lib = Path(__file__).resolve().parents[1] / "crates" / "bamboo-htslib" / "src" / "lib.rs"
    assert lib.exists()
    assert "phase-2-stub" in lib.read_text()


def test_phase4_polars_adapter(tiny_bam_path: Path) -> None:
    pl = pytest.importorskip("polars")

    table = bamboo.read_bam_table(str(tiny_bam_path), columns=["qname", "mapq"])
    df = bamboo.to_polars(table)
    assert len(df) == 2
    assert isinstance(df, pl.DataFrame)