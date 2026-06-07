"""pysam parity on realistic synthetic cohort data (50k reads, indexed region scans)."""

from __future__ import annotations

from pathlib import Path

import pytest

bamboo = pytest.importorskip("bamboo")
pysam = pytest.importorskip("pysam")

REGION = "chr1:1000000-5000000"
FULL_COLUMNS = ["qname", "flag", "rname", "pos", "mapq", "cigar", "seq", "qual"]


def _parse_region(region: str) -> tuple[str, int, int]:
    contig, interval = region.split(":", 1)
    start_s, end_s = interval.split("-", 1)
    return contig, max(int(start_s) - 1, 0), int(end_s)


def _pysam_column_dict(
    path: Path,
    columns: list[str],
    *,
    region: str | None = None,
    min_mapq: int | None = None,
    reference_filename: str | None = None,
) -> dict[str, list]:
    data: dict[str, list] = {name: [] for name in columns}
    mode = "rc" if path.suffix.lower() == ".cram" else "rb"
    open_kwargs: dict[str, str] = {}
    if reference_filename is not None:
        open_kwargs["reference_filename"] = reference_filename
    with pysam.AlignmentFile(str(path), mode, **open_kwargs) as bam:
        if region is None:
            iterator = bam
        else:
            contig, start, end = _parse_region(region)
            iterator = bam.fetch(contig, start, end)

        for read in iterator:
            if min_mapq is not None and (read.mapping_quality or 0) < min_mapq:
                continue
            for name in columns:
                if name == "qname":
                    data[name].append(read.query_name)
                elif name == "flag":
                    data[name].append(read.flag)
                elif name == "rname":
                    data[name].append(read.reference_name)
                elif name == "pos":
                    data[name].append(read.reference_start)
                elif name == "mapq":
                    data[name].append(read.mapping_quality)
                elif name == "cigar":
                    data[name].append(read.cigarstring)
                elif name == "seq":
                    data[name].append(read.query_sequence)
                elif name == "qual":
                    data[name].append(
                        read.query_qualities.tobytes().decode("ascii")
                        if read.query_qualities is not None
                        else None
                    )
                else:
                    raise AssertionError(f"unsupported column {name}")
    return data


def _arrow_dict(table) -> dict[str, list]:
    return {name: table.column(name).to_pylist() for name in table.column_names}


def test_bench_region_columnar_matches_pysam(bench_bam_with_index: Path) -> None:
    bamboo_table = bamboo.read_columns(
        str(bench_bam_with_index),
        columns=FULL_COLUMNS,
        region=REGION,
    )
    pysam_data = _pysam_column_dict(bench_bam_with_index, FULL_COLUMNS, region=REGION)
    assert _arrow_dict(bamboo_table) == pysam_data


def test_bench_region_mapq_filter_matches_pysam(bench_bam_with_index: Path) -> None:
    columns = ["qname", "pos", "mapq"]
    bamboo_table = bamboo.read_columns(
        str(bench_bam_with_index),
        columns=columns,
        region=REGION,
        min_mapq=30,
    )
    pysam_data = _pysam_column_dict(
        bench_bam_with_index,
        columns,
        region=REGION,
        min_mapq=30,
    )
    assert _arrow_dict(bamboo_table) == pysam_data


def test_bench_region_fetch_iterator_matches_pysam(bench_bam_with_index: Path) -> None:
    bamboo_rows: list[tuple] = []
    with bamboo.AlignmentFile(str(bench_bam_with_index)) as bam:
        for read in bam.fetch(region=REGION, min_mapq=20):
            bamboo_rows.append(
                (
                    read.query_name,
                    read.reference_name,
                    read.reference_start,
                    read.mapping_quality,
                    read.flag,
                    read.is_unmapped,
                )
            )

    pysam_rows: list[tuple] = []
    contig, start, end = _parse_region(REGION)
    with pysam.AlignmentFile(str(bench_bam_with_index), "rb") as bam:
        for read in bam.fetch(contig, start, end):
            if (read.mapping_quality or 0) < 20:
                continue
            pysam_rows.append(
                (
                    read.query_name,
                    read.reference_name,
                    read.reference_start,
                    read.mapping_quality,
                    read.flag,
                    read.is_unmapped,
                )
            )

    assert bamboo_rows == pysam_rows


def test_bench_cram_region_columnar_matches_pysam(
    bench_cram_with_index: Path,
    bench_fasta_with_index: Path,
) -> None:
    bamboo_table = bamboo.read_cram_columns(
        str(bench_cram_with_index),
        columns=FULL_COLUMNS,
        region=REGION,
        reference_filename=str(bench_fasta_with_index),
    )
    pysam_data = _pysam_column_dict(
        bench_cram_with_index,
        FULL_COLUMNS,
        region=REGION,
        reference_filename=str(bench_fasta_with_index),
    )
    assert _arrow_dict(bamboo_table) == pysam_data


def test_compat_import_region_count_matches_pysam(bench_bam_with_index: Path) -> None:
    from bamboo.compat import pysam as bamboo_pysam

    contig, start, end = _parse_region(REGION)
    with bamboo_pysam.AlignmentFile(str(bench_bam_with_index)) as bam:
        bamboo_count = sum(1 for _ in bam.fetch(contig, start, end))
    with pysam.AlignmentFile(str(bench_bam_with_index), "rb") as bam:
        pysam_count = sum(1 for _ in bam.fetch(contig, start, end))
    assert bamboo_count == pysam_count