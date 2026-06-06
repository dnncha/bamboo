"""Cross-validate Bamboo VCF outputs against pysam."""

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
def tiny_vcf_path() -> Path:
    path = Path(__file__).parent / "data" / "tiny.vcf"
    if not path.exists():
        pytest.skip("missing tests/data/tiny.vcf")
    return path


def test_vcf_count_matches_pysam(tiny_vcf_path: Path) -> None:
    with bamboo.VariantFile(str(tiny_vcf_path)) as vcf:
        bamboo_count = vcf.count()
    with pysam.VariantFile(str(tiny_vcf_path)) as vcf:
        pysam_count = sum(1 for _ in vcf)
    assert bamboo_count == pysam_count


@pytest.fixture(scope="session")
def tiny_vcf_gz_with_index() -> Path:
    path = Path(__file__).parent / "data" / "tiny.vcf.gz"
    if not path.exists() or not Path(f"{path}.tbi").exists():
        pytest.skip("missing indexed tiny.vcf.gz fixture")
    return path


def test_vcf_fetch_matches_pysam(tiny_vcf_gz_with_index: Path) -> None:
    region = "chr1:100-200"

    with bamboo.VariantFile(str(tiny_vcf_gz_with_index)) as vcf:
        bamboo_rows = [_variant_tuple(record) for record in vcf.fetch(region=region)]

    with pysam.VariantFile(str(tiny_vcf_gz_with_index)) as vcf:
        pysam_rows = [
            (rec.chrom, rec.pos, rec.ref, rec.alts[0] if rec.alts else "", rec.qual)
            for rec in vcf.fetch(region=region)
        ]

    assert bamboo_rows == pysam_rows


def test_vcf_arrow_matches_pysam(tiny_vcf_path: Path) -> None:
    pa = pytest.importorskip("pyarrow")
    columns = ["chrom", "pos", "ref", "alt", "qual"]

    bamboo_table = bamboo.read_vcf_table(str(tiny_vcf_path), columns=columns)
    assert isinstance(bamboo_table, pa.Table)

    expected = {name: [] for name in columns}
    with pysam.VariantFile(str(tiny_vcf_path)) as vcf:
        for rec in vcf:
            expected["chrom"].append(rec.chrom)
            expected["pos"].append(rec.pos)
            expected["ref"].append(rec.ref)
            expected["alt"].append(rec.alts[0] if rec.alts else "")
            expected["qual"].append(rec.qual)

    for name in columns:
        assert bamboo_table.column(name).to_pylist() == expected[name]