"""Cross-validate Bamboo outputs against pysam (scientific accuracy gate)."""

from __future__ import annotations

from pathlib import Path

import pytest

bamboo = pytest.importorskip("bamboo")
pysam = pytest.importorskip("pysam")

COLUMNS = ["qname", "rname", "pos", "mapq", "cigar"]


def _read_tuple(read: object) -> tuple[str | None, str | None, int | None, int | None, str]:
    return (
        read.query_name,
        read.reference_name,
        read.reference_start,
        read.mapping_quality,
        read.cigarstring,
    )


def _bamboo_rows(path: Path, *, region: str | None = None, min_mapq: int | None = None) -> list[tuple]:
    rows: list[tuple] = []
    with bamboo.AlignmentFile(str(path)) as bam:
        iterator = (
            bam.fetch(region=region, min_mapq=min_mapq)
            if region is not None or min_mapq is not None
            else bam
        )
        for read in iterator:
            rows.append(_read_tuple(read))
    return rows


def _pysam_rows(
    path: Path,
    *,
    region: str | None = None,
    min_mapq: int | None = None,
) -> list[tuple]:
    rows: list[tuple] = []
    with pysam.AlignmentFile(str(path), "rb") as bam:
        if region is not None:
            contig, interval = region.split(":", 1)
            start_s, end_s = interval.split("-", 1)
            start = max(int(start_s) - 1, 0)
            end = int(end_s)
            iterator = bam.fetch(contig, start, end)
        else:
            iterator = bam

        for read in iterator:
            if min_mapq is not None and (read.mapping_quality or 0) < min_mapq:
                continue
            rows.append(_read_tuple(read))
    return rows


def _arrow_dict(table) -> dict[str, list]:
    return {name: table.column(name).to_pylist() for name in table.column_names}


@pytest.mark.parametrize(
    "region",
    [
        None,
        "chr1:100-101",
        "chr1:50-150",
    ],
)
def test_fetch_rows_match_pysam(tiny_bam_with_index: Path, region: str | None) -> None:
    bamboo_rows = _bamboo_rows(tiny_bam_with_index, region=region)
    pysam_rows = _pysam_rows(tiny_bam_with_index, region=region)
    assert bamboo_rows == pysam_rows


def test_count_matches_pysam(tiny_bam_path: Path) -> None:
    with bamboo.AlignmentFile(str(tiny_bam_path)) as bam:
        bamboo_count = bam.count()
    with pysam.AlignmentFile(str(tiny_bam_path), "rb") as bam:
        pysam_count = sum(1 for _ in bam)
    assert bamboo_count == pysam_count


def test_mapq_filter_matches_pysam(tiny_bam_path: Path) -> None:
    bamboo_rows = _bamboo_rows(tiny_bam_path, min_mapq=30)
    pysam_rows = _pysam_rows(tiny_bam_path, min_mapq=30)
    assert bamboo_rows == pysam_rows


def test_region_fetch_arrow_matches_pysam(tiny_bam_with_index: Path) -> None:
    pa = pytest.importorskip("pyarrow")

    region = "chr1:100-101"
    with bamboo.AlignmentFile(str(tiny_bam_with_index)) as bam:
        bamboo_table = bam.fetch_arrow(columns=COLUMNS, region=region)
    assert isinstance(bamboo_table, pa.Table)

    contig, interval = region.split(":", 1)
    start_s, end_s = interval.split("-", 1)
    start = max(int(start_s) - 1, 0)
    end = int(end_s)

    columns: dict[str, list] = {name: [] for name in COLUMNS}
    with pysam.AlignmentFile(str(tiny_bam_with_index), "rb") as bam:
        for read in bam.fetch(contig, start, end):
            columns["qname"].append(read.query_name)
            columns["rname"].append(read.reference_name)
            columns["pos"].append(read.reference_start)
            columns["mapq"].append(read.mapping_quality)
            columns["cigar"].append(read.cigarstring)
    pysam_table = pa.table(columns)

    assert _arrow_dict(bamboo_table) == _arrow_dict(pysam_table)


def test_arrow_columns_match_pysam(tiny_bam_path: Path) -> None:
    pa = pytest.importorskip("pyarrow")

    bamboo_table = bamboo.read_columns(str(tiny_bam_path), columns=COLUMNS)
    assert isinstance(bamboo_table, pa.Table)

    pysam_table = pa.table(
        {
            "qname": [],
            "rname": [],
            "pos": [],
            "mapq": [],
            "cigar": [],
        }
    )
    columns: dict[str, list] = {name: [] for name in COLUMNS}
    with pysam.AlignmentFile(str(tiny_bam_path), "rb") as bam:
        for read in bam:
            columns["qname"].append(read.query_name)
            columns["rname"].append(read.reference_name)
            columns["pos"].append(read.reference_start)
            columns["mapq"].append(read.mapping_quality)
            columns["cigar"].append(read.cigarstring)
    pysam_table = pa.table(columns)

    assert _arrow_dict(bamboo_table) == _arrow_dict(pysam_table)


def test_write_roundtrip_matches_pysam(tiny_bam_path: Path, tmp_path: Path) -> None:
    bamboo_out = tmp_path / "bamboo_copy.bam"
    pysam_out = tmp_path / "pysam_copy.bam"

    with bamboo.AlignmentFile(str(tiny_bam_path)) as src:
        with bamboo.AlignmentFile(str(bamboo_out), "wb", template=src) as out:
            for read in src:
                out.write(read)

    with pysam.AlignmentFile(str(tiny_bam_path), "rb") as template_src:
        header = template_src.header.to_dict()
    with pysam.AlignmentFile(str(pysam_out), "wb", header=header) as out:
        with pysam.AlignmentFile(str(tiny_bam_path), "rb") as src:
            for read in src:
                out.write(read)

    assert _bamboo_rows(bamboo_out) == _pysam_rows(pysam_out)
    assert _bamboo_rows(bamboo_out) == _bamboo_rows(tiny_bam_path)


def test_bench_file_count_matches_pysam(bench_bam_with_index: Path) -> None:
    with bamboo.AlignmentFile(str(bench_bam_with_index)) as bam:
        bamboo_count = bam.count()
    with pysam.AlignmentFile(str(bench_bam_with_index), "rb") as bam:
        pysam_count = sum(1 for _ in bam)
    assert bamboo_count == pysam_count


def test_bench_region_fetch_count_matches_pysam(bench_bam_with_index: Path) -> None:
    region = "chr1:1000000-5000000"
    bamboo_rows = _bamboo_rows(bench_bam_with_index, region=region)
    pysam_rows = _pysam_rows(bench_bam_with_index, region=region)
    assert len(bamboo_rows) == len(pysam_rows)
    assert bamboo_rows == pysam_rows


def test_tiny_bam_flag_helpers_match_pysam(tiny_bam_path: Path) -> None:
    with bamboo.AlignmentFile(str(tiny_bam_path)) as bam:
        bamboo_flags = [(read.flag, read.is_paired, read.is_unmapped, read.is_reverse) for read in bam]
    with pysam.AlignmentFile(str(tiny_bam_path), "rb") as bam:
        pysam_flags = [(read.flag, read.is_paired, read.is_unmapped, read.is_reverse) for read in bam]
    assert bamboo_flags == pysam_flags