"""Roadmap phase gates — lightweight checks that each phase's deliverables exist."""

from __future__ import annotations

from pathlib import Path

import pytest

bamboo = pytest.importorskip("bamboo")


def test_phase1_bam_read_write_api() -> None:
    assert hasattr(bamboo, "AlignmentFile")
    assert hasattr(bamboo, "read_columns")
    assert hasattr(bamboo, "read_bam_table")


@pytest.fixture(scope="session")
def tiny_vcf_path() -> Path:
    path = Path(__file__).resolve().parent / "data" / "tiny.vcf"
    if not path.exists():
        pytest.skip("missing tests/data/tiny.vcf")
    return path


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


def test_phase2_indexed_vcf_fetch(tiny_vcf_gz_with_index: Path) -> None:
    with bamboo.VariantFile(str(tiny_vcf_gz_with_index)) as vcf:
        assert vcf.has_index()
        records = vcf.fetch(region="chr1:100-200")
    assert len(records) == 1


@pytest.fixture(scope="session")
def tiny_vcf_gz_with_index() -> Path:
    path = Path(__file__).resolve().parent / "data" / "tiny.vcf.gz"
    if not path.exists() or not Path(f"{path}.tbi").exists():
        pytest.skip("missing tests/data/tiny.vcf.gz or .tbi")
    return path


def test_phase3_vcf_reader(tiny_vcf_path: Path) -> None:
    pa = pytest.importorskip("pyarrow")

    table = bamboo.read_vcf_table(str(tiny_vcf_path), columns=["chrom", "pos", "ref", "alt"])
    assert isinstance(table, pa.Table)
    assert table.num_rows == 2


def test_phase2_cram_reader(tiny_cram_path: Path) -> None:
    with bamboo.CramFile(str(tiny_cram_path)) as cram:
        assert cram.count() == 2
        assert sum(1 for _ in cram) == 2


@pytest.fixture(scope="session")
def tiny_cram_path() -> Path:
    path = Path(__file__).resolve().parent / "data" / "tiny.cram"
    if not path.exists():
        pytest.skip("missing tests/data/tiny.cram")
    return path


def test_phase3_htslib_stub_present() -> None:
    lib = Path(__file__).resolve().parents[1] / "crates" / "bamboo-htslib" / "src" / "lib.rs"
    assert lib.exists()
    text = lib.read_text()
    assert "phase-2" in text
    assert "noodles" in text


def test_phase4_polars_adapter(tiny_bam_path: Path) -> None:
    pl = pytest.importorskip("polars")

    table = bamboo.read_bam_table(str(tiny_bam_path), columns=["qname", "mapq"])
    df = bamboo.to_polars(table)
    assert len(df) == 2
    assert isinstance(df, pl.DataFrame)