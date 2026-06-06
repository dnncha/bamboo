"""Integration tests for the Bamboo VCF reader."""

from __future__ import annotations

from pathlib import Path

import pytest

bamboo = pytest.importorskip("bamboo")


@pytest.fixture(scope="session")
def tiny_vcf_path() -> Path:
    path = Path(__file__).parent / "data" / "tiny.vcf"
    if not path.exists():
        pytest.skip(
            "missing tests/data/tiny.vcf — run "
            "`cargo run -p bamboo-noodles --example generate_fixtures`"
        )
    return path


def test_variant_file_count_and_references(tiny_vcf_path: Path) -> None:
    with bamboo.VariantFile(str(tiny_vcf_path)) as vcf:
        assert vcf.count() == 2
        assert vcf.references() == ["chr1", "chr2"]


def test_variant_file_fetch_region(tiny_vcf_path: Path) -> None:
    with bamboo.VariantFile(str(tiny_vcf_path)) as vcf:
        records = vcf.fetch(region="chr1:100-101")
    assert len(records) == 1
    assert records[0].chrom == "chr1"
    assert records[0].pos == 100
    assert records[0].ref_ == "A"
    assert records[0].alt == "G"


@pytest.fixture(scope="session")
def tiny_vcf_gz_with_index() -> Path:
    path = Path(__file__).parent / "data" / "tiny.vcf.gz"
    index_path = Path(f"{path}.tbi")
    if not path.exists() or not index_path.exists():
        pytest.skip("missing tests/data/tiny.vcf.gz or .tbi")
    return path


def test_variant_file_has_index(tiny_vcf_gz_with_index: Path) -> None:
    with bamboo.VariantFile(str(tiny_vcf_gz_with_index)) as vcf:
        assert vcf.has_index()


def test_indexed_vcf_fetch(tiny_vcf_gz_with_index: Path) -> None:
    with bamboo.VariantFile(str(tiny_vcf_gz_with_index)) as vcf:
        records = vcf.fetch(region="chr1:100-200")
    assert len(records) == 1
    assert records[0].chrom == "chr1"
    assert records[0].pos == 100


def test_read_vcf_table_exports_arrow(tiny_vcf_path: Path) -> None:
    pa = pytest.importorskip("pyarrow")

    table = bamboo.read_vcf_table(
        str(tiny_vcf_path),
        columns=["chrom", "pos", "ref", "alt"],
    )
    assert isinstance(table, pa.Table)
    assert table.num_rows == 2
    assert table.column("chrom")[0].as_py() == "chr1"