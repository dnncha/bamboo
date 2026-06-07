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
    assert 'name = "bamboo-hts"' in text
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


def test_phase3_bcf_reader(tiny_bcf_path: Path) -> None:
    pa = pytest.importorskip("pyarrow")

    table = bamboo.read_bcf_table(str(tiny_bcf_path), columns=["chrom", "pos", "ref", "alt"])
    assert isinstance(table, pa.Table)
    assert table.num_rows == 2


def test_phase3_bcf_indexed_fetch(tiny_bcf_with_index: Path) -> None:
    with bamboo.VariantFile(str(tiny_bcf_with_index)) as bcf:
        assert bcf.has_index()
        records = bcf.fetch(region="chr1:100-200")
    assert len(records) == 1
    assert records[0].chrom == "chr1"
    assert records[0].pos == 100


@pytest.fixture(scope="session")
def tiny_bcf_path() -> Path:
    path = Path(__file__).resolve().parent / "data" / "tiny.bcf"
    if not path.exists():
        pytest.skip("missing tests/data/tiny.bcf")
    return path


@pytest.fixture(scope="session")
def tiny_bcf_with_index(tiny_bcf_path: Path) -> Path:
    if not Path(f"{tiny_bcf_path}.csi").exists():
        pytest.skip("missing tests/data/tiny.bcf.csi")
    return tiny_bcf_path


def test_phase3_vcf_reader(tiny_vcf_path: Path) -> None:
    pa = pytest.importorskip("pyarrow")

    table = bamboo.read_vcf_table(str(tiny_vcf_path), columns=["chrom", "pos", "ref", "alt"])
    assert isinstance(table, pa.Table)
    assert table.num_rows == 2


def test_phase2_cram_reader(tiny_cram_with_index: Path) -> None:
    with bamboo.CramFile(str(tiny_cram_with_index)) as cram:
        assert cram.count() == 2
        assert cram.has_index()
        assert sum(1 for _ in cram) == 2
        assert len(list(cram.fetch(region="chr1:100-101"))) == 1


def test_phase2_cram_columnar_with_indexed_fasta(
    tiny_cram_with_index: Path,
    tiny_fasta_with_index: Path,
) -> None:
    pa = pytest.importorskip("pyarrow")

    table = bamboo.read_cram_columns(
        str(tiny_cram_with_index),
        columns=["qname", "pos"],
        region="chr1:100-101",
        reference_filename=str(tiny_fasta_with_index),
    )
    assert isinstance(table, pa.Table)
    assert table.num_rows == 1


@pytest.fixture(scope="session")
def tiny_fasta_with_index() -> Path:
    path = Path(__file__).resolve().parent / "data" / "tiny.fasta"
    if not path.exists() or not Path(f"{path}.fai").exists():
        pytest.skip("missing tests/data/tiny.fasta or .fai")
    return path


def test_phase2_cram_columnar(tiny_cram_with_index: Path) -> None:
    pa = pytest.importorskip("pyarrow")

    table = bamboo.read_cram_columns(
        str(tiny_cram_with_index),
        columns=["qname", "pos"],
        region="chr1:100-101",
    )
    assert isinstance(table, pa.Table)
    assert table.num_rows == 1


@pytest.fixture(scope="session")
def tiny_cram_with_index() -> Path:
    path = Path(__file__).resolve().parent / "data" / "tiny.cram"
    if not path.exists() or not Path(f"{path}.crai").exists():
        pytest.skip("missing tests/data/tiny.cram or .crai")
    return path


def test_phase3_htslib_stub_present() -> None:
    lib = Path(__file__).resolve().parents[1] / "crates" / "bamboo-htslib" / "src" / "lib.rs"
    assert lib.exists()
    text = lib.read_text()
    assert "phase-2" in text
    assert "noodles" in text
    assert "pileup_available" in text


def test_phase4_polars_adapter(tiny_bam_path: Path) -> None:
    pl = pytest.importorskip("polars")

    table = bamboo.read_bam_table(str(tiny_bam_path), columns=["qname", "mapq"])
    df = bamboo.to_polars(table)
    assert len(df) == 2
    assert isinstance(df, pl.DataFrame)


def test_phase4_pandas_adapter(tiny_bam_path: Path) -> None:
    pd = pytest.importorskip("pandas")

    table = bamboo.read_bam_table(str(tiny_bam_path), columns=["qname", "mapq"])
    df = bamboo.to_pandas(table)
    assert len(df) == 2
    assert isinstance(df, pd.DataFrame)


def test_phase3_htslib_capability_probes() -> None:
    assert bamboo.primary_backend() == "noodles"
    assert isinstance(bamboo.htslib_available(), bool)
    assert isinstance(bamboo.pileup_available(), bool)


def test_phase3_cram_pileup_with_reference(
    tiny_cram_with_index: Path,
    tiny_fasta_with_index: Path,
) -> None:
    if not bamboo.pileup_available():
        pytest.skip("pileup requires htslib build")

    with bamboo.CramFile(
        str(tiny_cram_with_index),
        reference_filename=str(tiny_fasta_with_index),
    ) as cram:
        columns = list(cram.pileup(region="chr1:100-101"))
    assert len(columns) >= 1
    assert columns[0].n >= 1