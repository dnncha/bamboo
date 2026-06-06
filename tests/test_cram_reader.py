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