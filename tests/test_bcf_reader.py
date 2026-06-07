"""Integration tests for the Bamboo BCF reader."""

from __future__ import annotations

from pathlib import Path

import pytest

bamboo = pytest.importorskip("bamboo")
pysam = pytest.importorskip("pysam")


def _variant_tuple(record: object) -> tuple[str, int, str, str, str | None]:
    return (
        record.chrom,
        record.pos,
        record.ref_,
        record.alt,
        record.qual,
    )


@pytest.fixture(scope="session")
def tiny_bcf_path() -> Path:
    path = Path(__file__).parent / "data" / "tiny.bcf"
    if not path.exists():
        pytest.skip(
            "missing tests/data/tiny.bcf — run "
            "`cargo run -p bamboo-noodles --example generate_fixtures`"
        )
    return path


@pytest.fixture(scope="session")
def tiny_bcf_with_index(tiny_bcf_path: Path) -> Path:
    if not Path(f"{tiny_bcf_path}.csi").exists():
        pytest.skip("missing tests/data/tiny.bcf.csi")
    return tiny_bcf_path


@pytest.fixture(scope="session")
def tiny_vcf_path() -> Path:
    path = Path(__file__).parent / "data" / "tiny.vcf"
    if not path.exists():
        pytest.skip("missing tests/data/tiny.vcf")
    return path


def test_bcf_count_matches_pysam(tiny_bcf_path: Path) -> None:
    with bamboo.VariantFile(str(tiny_bcf_path)) as bcf:
        bamboo_count = bcf.count()
    with pysam.VariantFile(str(tiny_bcf_path)) as bcf:
        pysam_count = sum(1 for _ in bcf)
    assert bamboo_count == pysam_count


def test_bcf_arrow_matches_vcf(tiny_bcf_path: Path, tiny_vcf_path: Path) -> None:
    pa = pytest.importorskip("pyarrow")
    columns = ["chrom", "pos", "ref", "alt", "qual"]

    bcf_table = bamboo.read_bcf_table(str(tiny_bcf_path), columns=columns)
    vcf_table = bamboo.read_vcf_table(str(tiny_vcf_path), columns=columns)
    assert isinstance(bcf_table, pa.Table)
    assert bcf_table.to_pydict() == vcf_table.to_pydict()


def test_bcf_fetch_matches_pysam(tiny_bcf_path: Path) -> None:
    region = "chr1:100-200"
    contig, interval = region.split(":", 1)
    start_s, end_s = interval.split("-", 1)
    start = int(start_s)
    end = int(end_s)

    with bamboo.VariantFile(str(tiny_bcf_path)) as bcf:
        bamboo_rows = [_variant_tuple(record) for record in bcf.fetch(region=region)]

    with pysam.VariantFile(str(tiny_bcf_path)) as bcf:
        pysam_rows = [
            (rec.chrom, rec.pos, rec.ref, rec.alts[0] if rec.alts else "", rec.qual)
            for rec in bcf
            if rec.chrom == contig and start <= rec.pos <= end
        ]

    assert bamboo_rows == pysam_rows


def test_bcf_has_csi_index(tiny_bcf_with_index: Path) -> None:
    with bamboo.VariantFile(str(tiny_bcf_with_index)) as bcf:
        assert bcf.has_index()


def test_bcf_indexed_fetch_matches_pysam(tiny_bcf_with_index: Path) -> None:
    region = "chr1:100-200"
    contig, interval = region.split(":", 1)
    start_s, end_s = interval.split("-", 1)
    start = max(int(start_s) - 1, 0)
    end = int(end_s)

    with bamboo.VariantFile(str(tiny_bcf_with_index)) as bcf:
        assert bcf.has_index()
        bamboo_rows = [_variant_tuple(record) for record in bcf.fetch(region=region)]

    with pysam.VariantFile(str(tiny_bcf_with_index)) as bcf:
        pysam_rows = [
            (rec.chrom, rec.pos, rec.ref, rec.alts[0] if rec.alts else "", rec.qual)
            for rec in bcf.fetch(contig, start, end)
        ]

    assert bamboo_rows == pysam_rows