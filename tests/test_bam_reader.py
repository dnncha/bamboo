"""Integration tests for the Bamboo BAM reader."""

from __future__ import annotations

from pathlib import Path

import pytest

bamboo = pytest.importorskip("bamboo")


def test_alignment_file_count_and_fetch(tiny_bam_path: Path) -> None:
    with bamboo.AlignmentFile(str(tiny_bam_path)) as bam:
        assert bam.count() == 2
        assert bam.references() == ["chr1", "chr2"]
        assert bam.reference_lengths == [1000, 1000]
        assert bam.header() == {"chr1": 1000, "chr2": 1000}

        reads = list(bam.fetch(contig="chr1", start=50, stop=150))
        assert len(reads) == 1
        assert reads[0].query_name == "read1"
        assert reads[0].reference_name == "chr1"
        assert reads[0].reference_start == 99
        assert reads[0].cigarstring == "6M"


def test_alignment_file_indexed_fetch(tiny_bam_with_index: Path) -> None:
    with bamboo.AlignmentFile(str(tiny_bam_with_index)) as bam:
        assert bam.has_index()
        reads = list(bam.fetch(region="chr1:100-101"))
        assert len(reads) == 1
        assert reads[0].query_name == "read1"


def test_read_bam_table_exports_arrow(tiny_bam_path: Path) -> None:
    pa = pytest.importorskip("pyarrow")

    table = bamboo.read_bam_table(
        str(tiny_bam_path),
        columns=["qname", "rname", "pos", "mapq"],
        region="chr1:100-101",
    )
    assert isinstance(table, pa.Table)
    assert table.num_rows == 1
    assert table.column("qname")[0].as_py() == "read1"


def test_alignment_file_to_arrow(tiny_bam_path: Path) -> None:
    pa = pytest.importorskip("pyarrow")

    with bamboo.AlignmentFile(str(tiny_bam_path)) as bam:
        table = bam.to_arrow(columns=["qname", "rname", "pos"], min_mapq=30)

    assert table.num_rows == 1
    assert table.column("qname")[0].as_py() == "read1"


def test_to_polars_helper(tiny_bam_path: Path) -> None:
    pl = pytest.importorskip("polars")
    pa = pytest.importorskip("pyarrow")

    table = bamboo.read_bam_table(str(tiny_bam_path), columns=["qname", "mapq"])
    df = bamboo.to_polars(table)
    assert len(df) == 2
    assert df.columns == ["qname", "mapq"]