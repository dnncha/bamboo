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


def test_alignment_file_write_round_trip(tiny_bam_path: Path, tmp_path: Path) -> None:
    out_path = tmp_path / "copy.bam"

    with bamboo.AlignmentFile(str(tiny_bam_path)) as src:
        with bamboo.AlignmentFile(str(out_path), "wb", template=src) as out:
            assert out.mode == "wb"
            for read in src:
                out.write(read)

    with bamboo.AlignmentFile(str(out_path)) as copied:
        assert copied.count() == 2
        assert copied.references() == ["chr1", "chr2"]
        reads = list(copied.fetch(contig="chr1", start=50, stop=150))
        assert len(reads) == 1
        assert reads[0].query_name == "read1"


def test_alignment_file_opens_file_uri(tiny_bam_with_index: Path) -> None:
    uri = f"file://{tiny_bam_with_index}"
    with bamboo.AlignmentFile(uri) as bam:
        assert bam.filename() == uri
        assert bam.count() == 2
        assert bam.has_index()


def test_read_columns_fast_path(tiny_bam_path: Path) -> None:
    pa = pytest.importorskip("pyarrow")

    table = bamboo.read_columns(
        str(tiny_bam_path),
        columns=["qname", "rname", "pos", "mapq", "cigar"],
    )
    assert isinstance(table, pa.Table)
    assert table.num_rows == 2
    assert table.column("qname")[0].as_py() == "read1"


def test_pysam_compat_imports(tiny_bam_path: Path) -> None:
    from bamboo.compat import pysam as pysam_shim

    assert pysam_shim.AlignmentFile is bamboo.AlignmentFile
    assert pysam_shim.alignment_count(str(tiny_bam_path)) == 2


def test_to_polars_helper(tiny_bam_path: Path) -> None:
    pl = pytest.importorskip("polars")
    pa = pytest.importorskip("pyarrow")

    table = bamboo.read_bam_table(str(tiny_bam_path), columns=["qname", "mapq"])
    df = bamboo.to_polars(table)
    assert len(df) == 2
    assert df.columns == ["qname", "mapq"]