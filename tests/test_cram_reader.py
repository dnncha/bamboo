"""Integration tests for the Bamboo CRAM reader."""

from __future__ import annotations

from pathlib import Path

import pytest

bamboo = pytest.importorskip("bamboo")

COLUMNS = ["qname", "rname", "pos", "mapq", "cigar"]


def _read_tuple(read: object) -> tuple[str | None, str | None, int | None, int | None, str]:
    return (
        read.query_name,
        read.reference_name,
        read.reference_start,
        read.mapping_quality,
        read.cigarstring,
    )


@pytest.fixture(scope="session")
def tiny_cram_path() -> Path:
    path = Path(__file__).parent / "data" / "tiny.cram"
    if not path.exists():
        pytest.skip(
            "missing tests/data/tiny.cram — run "
            "`cargo run -p bamboo-noodles --example generate_fixtures`"
        )
    return path


@pytest.fixture(scope="session")
def tiny_cram_with_index(tiny_cram_path: Path) -> Path:
    if not Path(f"{tiny_cram_path}.crai").exists():
        pytest.skip("missing tests/data/tiny.cram.crai")
    return tiny_cram_path


@pytest.fixture(scope="session")
def tiny_fasta_path() -> Path:
    path = Path(__file__).parent / "data" / "tiny.fasta"
    if not path.exists():
        pytest.skip("missing tests/data/tiny.fasta")
    return path


@pytest.fixture(scope="session")
def tiny_bam_path() -> Path:
    path = Path(__file__).parent / "data" / "tiny.bam"
    if not path.exists():
        pytest.skip("missing tests/data/tiny.bam")
    return path


def test_cram_file_count_and_references(tiny_cram_path: Path) -> None:
    with bamboo.CramFile(str(tiny_cram_path)) as cram:
        assert cram.count() == 2
        assert cram.references() == ["chr1", "chr2"]
        assert cram.reference_lengths == [1000, 1000]
        assert cram.header() == {"chr1": 1000, "chr2": 1000}


def test_cram_has_index(tiny_cram_with_index: Path) -> None:
    with bamboo.CramFile(str(tiny_cram_with_index)) as cram:
        assert cram.has_index()


def test_cram_iteration_matches_bam(tiny_cram_path: Path, tiny_bam_path: Path) -> None:
    with bamboo.CramFile(str(tiny_cram_path)) as cram:
        cram_rows = [_read_tuple(read) for read in cram]
    with bamboo.AlignmentFile(str(tiny_bam_path)) as bam:
        bam_rows = [_read_tuple(read) for read in bam]
    assert cram_rows == bam_rows


def test_cram_fetch_region(tiny_cram_path: Path) -> None:
    with bamboo.CramFile(str(tiny_cram_path)) as cram:
        rows = [_read_tuple(read) for read in cram.fetch(region="chr1:100-101")]
    assert len(rows) == 1
    assert rows[0][0] == "read1"


def test_cram_indexed_fetch_region(tiny_cram_with_index: Path) -> None:
    with bamboo.CramFile(str(tiny_cram_with_index)) as cram:
        assert cram.has_index()
        indexed_rows = [_read_tuple(read) for read in cram.fetch(region="chr1:100-101")]
        linear_rows = [
            _read_tuple(read)
            for read in cram
            if read.reference_name == "chr1" and read.reference_start == 99
        ]
    assert indexed_rows == linear_rows
    assert len(indexed_rows) == 1


def test_cram_reads_with_external_reference(
    tiny_cram_path: Path,
    tiny_fasta_path: Path,
) -> None:
    with bamboo.CramFile(str(tiny_cram_path), reference_filename=str(tiny_fasta_path)) as cram:
        rows = [_read_tuple(read) for read in cram]
    assert len(rows) == 2