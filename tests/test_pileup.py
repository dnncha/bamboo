"""Pileup parity against pysam when htslib support is compiled in."""

from __future__ import annotations

from pathlib import Path

import pytest

bamboo = pytest.importorskip("bamboo")
pysam = pytest.importorskip("pysam")

pytestmark = pytest.mark.skipif(
    not bamboo.pileup_available(),
    reason="Bamboo built without htslib pileup support",
)


def _bamboo_pileup(path: Path, contig: str, start: int, end: int) -> list[tuple[int, int, list[str | None]]]:
    rows: list[tuple[int, int, list[str | None]]] = []
    with bamboo.AlignmentFile(str(path)) as bam:
        for column in bam.pileup(contig, start, end):
            names = [read.query_name for read in column.pileups]
            rows.append((column.pos, column.n, names))
    return rows


def _pysam_pileup(path: Path, contig: str, start: int, end: int) -> list[tuple[int, int, list[str | None]]]:
    rows: list[tuple[int, int, list[str | None]]] = []
    with pysam.AlignmentFile(str(path), "rb") as bam:
        for column in bam.pileup(contig, start, end):
            names = [pileupread.alignment.query_name for pileupread in column.pileups]
            rows.append((column.pos, column.n, names))
    return rows


def test_pileup_reads_false_reports_depth_in_loop(tiny_bam_with_index: Path) -> None:
    depths: list[int] = []
    with bamboo.AlignmentFile(str(tiny_bam_with_index)) as bam:
        for column in bam.pileup("chr1", 99, 101, reads=False):
            depths.append(column.n)
    assert depths == [1, 1, 1, 1, 1, 1]


def test_pileup_matches_pysam_over_tiny_region(tiny_bam_with_index: Path) -> None:
    bamboo_rows = _bamboo_pileup(tiny_bam_with_index, "chr1", 99, 101)
    pysam_rows = _pysam_pileup(tiny_bam_with_index, "chr1", 99, 101)
    assert bamboo_rows == pysam_rows


def test_pileup_region_kwarg_matches_pysam(tiny_bam_with_index: Path) -> None:
    with bamboo.AlignmentFile(str(tiny_bam_with_index)) as bam:
        bamboo_rows = [
            (column.pos, column.n, [read.query_name for read in column.pileups])
            for column in bam.pileup(region="chr1:100-101")
        ]

    contig, interval = "chr1:100-101".split(":", 1)
    start_s, end_s = interval.split("-", 1)
    start = max(int(start_s) - 1, 0)
    end = int(end_s)

    with pysam.AlignmentFile(str(tiny_bam_with_index), "rb") as bam:
        pysam_rows = [
            (column.pos, column.n, [pileupread.alignment.query_name for pileupread in column.pileups])
            for column in bam.pileup(contig, start, end)
        ]

    assert bamboo_rows == pysam_rows


@pytest.fixture(scope="session")
def tiny_cram_with_index() -> Path:
    path = Path(__file__).resolve().parent / "data" / "tiny.cram"
    if not path.exists() or not Path(f"{path}.crai").exists():
        pytest.skip("missing tests/data/tiny.cram or .crai")
    return path


@pytest.fixture(scope="session")
def tiny_fasta_path() -> Path:
    path = Path(__file__).resolve().parent / "data" / "tiny.fasta"
    if not path.exists():
        pytest.skip("missing tests/data/tiny.fasta")
    return path


def _bamboo_cram_pileup(
    cram_path: Path,
    fasta_path: Path,
    contig: str,
    start: int,
    end: int,
) -> list[tuple[int, int, list[str | None]]]:
    rows: list[tuple[int, int, list[str | None]]] = []
    with bamboo.CramFile(str(cram_path), reference_filename=str(fasta_path)) as cram:
        for column in cram.pileup(contig, start, end):
            names = [read.query_name for read in column.pileups]
            rows.append((column.pos, column.n, names))
    return rows


def _pysam_cram_pileup(
    cram_path: Path,
    fasta_path: Path,
    contig: str,
    start: int,
    end: int,
) -> list[tuple[int, int, list[str | None]]]:
    rows: list[tuple[int, int, list[str | None]]] = []
    with pysam.AlignmentFile(
        str(cram_path),
        "rc",
        reference_filename=str(fasta_path),
    ) as cram:
        for column in cram.pileup(contig, start, end):
            names = [pileupread.alignment.query_name for pileupread in column.pileups]
            rows.append((column.pos, column.n, names))
    return rows


def test_cram_pileup_matches_pysam(
    tiny_cram_with_index: Path,
    tiny_fasta_path: Path,
) -> None:
    bamboo_rows = _bamboo_cram_pileup(tiny_cram_with_index, tiny_fasta_path, "chr1", 99, 101)
    pysam_rows = _pysam_cram_pileup(tiny_cram_with_index, tiny_fasta_path, "chr1", 99, 101)
    assert bamboo_rows == pysam_rows


def test_cram_pileup_region_kwarg_matches_pysam(
    tiny_cram_with_index: Path,
    tiny_fasta_path: Path,
) -> None:
    with bamboo.CramFile(
        str(tiny_cram_with_index),
        reference_filename=str(tiny_fasta_path),
    ) as cram:
        bamboo_rows = [
            (column.pos, column.n, [read.query_name for read in column.pileups])
            for column in cram.pileup(region="chr1:100-101")
        ]

    contig, interval = "chr1:100-101".split(":", 1)
    start_s, end_s = interval.split("-", 1)
    start = max(int(start_s) - 1, 0)
    end = int(end_s)

    with pysam.AlignmentFile(
        str(tiny_cram_with_index),
        "rc",
        reference_filename=str(tiny_fasta_path),
    ) as cram:
        pysam_rows = [
            (column.pos, column.n, [pileupread.alignment.query_name for pileupread in column.pileups])
            for column in cram.pileup(contig, start, end)
        ]

    assert bamboo_rows == pysam_rows