"""Benchmark task implementations for Bamboo and competitors."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Callable

import pyarrow as pa

REGION = "chr1:1000000-5000000"
ARROW_COLUMNS = ["qname", "rname", "pos", "mapq", "cigar"]


def available_backends() -> list[str]:
    backends = ["bamboo"]
    try:
        import pysam  # noqa: F401
    except ImportError:
        pass
    else:
        backends.append("pysam")
    if shutil.which("samtools"):
        backends.append("samtools")
    return backends


def _materialize_fields(read: object) -> tuple[str | None, str | None, int | None, int | None, str]:
    if hasattr(read, "query_name"):
        qname = read.query_name
        rname = read.reference_name
        pos = read.reference_start
        mapq = read.mapping_quality
        cigar = read.cigarstring
    else:
        qname = read.query_name if read.query_name else None
        rname = read.reference_name if read.reference_name else None
        pos = read.reference_start
        mapq = read.mapping_quality
        cigar = read.cigarstring
    return qname, rname, pos, mapq, cigar


def task_count_records(path: Path, backend: str) -> Callable[[], int]:
    if backend == "bamboo":
        import bamboo as bm

        def run() -> int:
            with bm.AlignmentFile(str(path)) as bam:
                return bam.count()

        return run

    if backend == "pysam":
        import pysam

        def run() -> int:
            total = 0
            with pysam.AlignmentFile(str(path), "rb") as bam:
                for _ in bam:
                    total += 1
            return total

        return run

    if backend == "samtools":
        def run() -> int:
            completed = subprocess.run(
                ["samtools", "view", "-c", str(path)],
                check=True,
                capture_output=True,
                text=True,
            )
            return int(completed.stdout.strip())

        return run

    raise ValueError(f"unsupported backend '{backend}'")


def task_iterate_materialize(path: Path, backend: str) -> Callable[[], int]:
    if backend == "bamboo":
        import bamboo as bm

        def run() -> int:
            rows: list[tuple[str | None, str | None, int | None, int | None, str]] = []
            with bm.AlignmentFile(str(path)) as bam:
                for read in bam:
                    rows.append(_materialize_fields(read))
            return len(rows)

        return run

    if backend == "pysam":
        import pysam

        def run() -> int:
            rows: list[tuple[str | None, str | None, int | None, int | None, str]] = []
            with pysam.AlignmentFile(str(path), "rb") as bam:
                for read in bam:
                    rows.append(_materialize_fields(read))
            return len(rows)

        return run

    if backend == "samtools":
        def run() -> int:
            completed = subprocess.run(
                ["samtools", "view", str(path)],
                check=True,
                capture_output=True,
                text=True,
            )
            return len(completed.stdout.splitlines())

        return run

    raise ValueError(f"unsupported backend '{backend}'")


def task_region_fetch(path: Path, backend: str, region: str = REGION) -> Callable[[], int]:
    if backend == "bamboo":
        import bamboo as bm

        def run() -> int:
            rows: list[tuple[str | None, str | None, int | None, int | None, str]] = []
            with bm.AlignmentFile(str(path)) as bam:
                for read in bam.fetch(region=region):
                    rows.append(_materialize_fields(read))
            return len(rows)

        return run

    if backend == "pysam":
        import pysam

        contig, interval = region.split(":", 1)
        start_s, end_s = interval.split("-", 1)
        start = max(int(start_s) - 1, 0)
        end = int(end_s)

        def run() -> int:
            rows: list[tuple[str | None, str | None, int | None, int | None, str]] = []
            with pysam.AlignmentFile(str(path), "rb") as bam:
                for read in bam.fetch(contig, start, end):
                    rows.append(_materialize_fields(read))
            return len(rows)

        return run

    if backend == "samtools":
        def run() -> int:
            completed = subprocess.run(
                ["samtools", "view", str(path), region],
                check=True,
                capture_output=True,
                text=True,
            )
            return len(completed.stdout.splitlines())

        return run

    raise ValueError(f"unsupported backend '{backend}'")


def task_arrow_export(path: Path, backend: str) -> Callable[[], int]:
    if backend == "bamboo":
        import bamboo as bm

        def run() -> int:
            with bm.AlignmentFile(str(path)) as bam:
                table = bam.to_arrow(columns=ARROW_COLUMNS)
            if not isinstance(table, pa.Table):
                raise TypeError("expected pyarrow.Table from bamboo.to_arrow")
            return table.num_rows

        return run

    if backend == "pysam":
        import pysam

        def run() -> int:
            qname: list[str | None] = []
            rname: list[str | None] = []
            pos: list[int | None] = []
            mapq: list[int | None] = []
            cigar: list[str] = []
            with pysam.AlignmentFile(str(path), "rb") as bam:
                for read in bam:
                    qname.append(read.query_name)
                    rname.append(read.reference_name)
                    pos.append(read.reference_start)
                    mapq.append(read.mapping_quality)
                    cigar.append(read.cigarstring)
            table = pa.table(
                {
                    "qname": qname,
                    "rname": rname,
                    "pos": pos,
                    "mapq": mapq,
                    "cigar": cigar,
                }
            )
            return table.num_rows

        return run

    if backend == "samtools":
        def run() -> int:
            raise RuntimeError("samtools does not provide a native Arrow export path")

        return run

    raise ValueError(f"unsupported backend '{backend}'")


def task_write_roundtrip(path: Path, backend: str) -> Callable[[], int]:
    if backend == "bamboo":
        import bamboo as bm

        def run() -> int:
            with tempfile.TemporaryDirectory() as tmp:
                out_path = Path(tmp) / "roundtrip.bam"
                with bm.AlignmentFile(str(path)) as src:
                    with bm.AlignmentFile(str(out_path), "wb", template=src) as out:
                        for read in src:
                            out.write(read)
                with bm.AlignmentFile(str(out_path)) as copied:
                    return copied.count()

        return run

    if backend == "pysam":
        import pysam

        def run() -> int:
            with tempfile.TemporaryDirectory() as tmp:
                out_path = Path(tmp) / "roundtrip.bam"
                with pysam.AlignmentFile(str(path), "rb") as template_src:
                    header = template_src.header
                with pysam.AlignmentFile(str(out_path), "wb", header=header) as out:
                    with pysam.AlignmentFile(str(path), "rb") as src:
                        for read in src:
                            out.write(read)
                with pysam.AlignmentFile(str(out_path), "rb") as copied:
                    return sum(1 for _ in copied)

        return run

    if backend == "samtools":
        def run() -> int:
            raise RuntimeError("samtools write roundtrip is not part of this benchmark")

        return run

    raise ValueError(f"unsupported backend '{backend}'")


def task_region_fetch_bulk(path: Path, backend: str, region: str = REGION) -> Callable[[], int]:
    if backend == "bamboo":
        import bamboo as bm

        def run() -> int:
            with bm.AlignmentFile(str(path)) as bam:
                reads = bam.fetch_bulk(region=region)
            rows = [_materialize_fields(read) for read in reads]
            return len(rows)

        return run

    if backend == "pysam":
        import pysam

        contig, interval = region.split(":", 1)
        start_s, end_s = interval.split("-", 1)
        start = max(int(start_s) - 1, 0)
        end = int(end_s)

        def run() -> int:
            rows: list[tuple[str | None, str | None, int | None, int | None, str]] = []
            with pysam.AlignmentFile(str(path), "rb") as bam:
                for read in bam.fetch(contig, start, end):
                    rows.append(_materialize_fields(read))
            return len(rows)

        return run

    if backend == "samtools":
        def run() -> int:
            completed = subprocess.run(
                ["samtools", "view", str(path), region],
                check=True,
                capture_output=True,
                text=True,
            )
            return len(completed.stdout.splitlines())

        return run

    raise ValueError(f"unsupported backend '{backend}'")


def task_region_fetch_arrow(path: Path, backend: str, region: str = REGION) -> Callable[[], int]:
    if backend == "bamboo":
        import bamboo as bm

        def run() -> int:
            with bm.AlignmentFile(str(path)) as bam:
                table = bam.fetch_arrow(columns=ARROW_COLUMNS, region=region)
            if not isinstance(table, pa.Table):
                raise TypeError("expected pyarrow.Table from bamboo.fetch_arrow")
            return table.num_rows

        return run

    if backend == "pysam":
        import pysam

        contig, interval = region.split(":", 1)
        start_s, end_s = interval.split("-", 1)
        start = max(int(start_s) - 1, 0)
        end = int(end_s)

        def run() -> int:
            qname: list[str | None] = []
            rname: list[str | None] = []
            pos: list[int | None] = []
            mapq: list[int | None] = []
            cigar: list[str] = []
            with pysam.AlignmentFile(str(path), "rb") as bam:
                for read in bam.fetch(contig, start, end):
                    qname.append(read.query_name)
                    rname.append(read.reference_name)
                    pos.append(read.reference_start)
                    mapq.append(read.mapping_quality)
                    cigar.append(read.cigarstring)
            table = pa.table(
                {
                    "qname": qname,
                    "rname": rname,
                    "pos": pos,
                    "mapq": mapq,
                    "cigar": cigar,
                }
            )
            return table.num_rows

        return run

    if backend == "samtools":
        def run() -> int:
            raise RuntimeError("samtools has no native region Arrow export")

        return run

    raise ValueError(f"unsupported backend '{backend}'")


def task_region_columnar(path: Path, backend: str, region: str = REGION) -> Callable[[], int]:
    if backend == "bamboo":
        import bamboo as bm

        def run() -> int:
            table = bm.read_columns(
                str(path),
                columns=ARROW_COLUMNS,
                region=region,
            )
            return table.num_rows

        return run

    if backend == "pysam":
        import pysam

        contig, interval = region.split(":", 1)
        start_s, end_s = interval.split("-", 1)
        start = max(int(start_s) - 1, 0)
        end = int(end_s)

        def run() -> int:
            qname: list[str | None] = []
            rname: list[str | None] = []
            pos: list[int | None] = []
            mapq: list[int | None] = []
            cigar: list[str] = []
            with pysam.AlignmentFile(str(path), "rb") as bam:
                for read in bam.fetch(contig, start, end):
                    qname.append(read.query_name)
                    rname.append(read.reference_name)
                    pos.append(read.reference_start)
                    mapq.append(read.mapping_quality)
                    cigar.append(read.cigarstring)
            table = pa.table(
                {
                    "qname": qname,
                    "rname": rname,
                    "pos": pos,
                    "mapq": mapq,
                    "cigar": cigar,
                }
            )
            return table.num_rows

        return run

    if backend == "samtools":
        def run() -> int:
            raise RuntimeError("samtools has no native region columnar export")

        return run

    raise ValueError(f"unsupported backend '{backend}'")


def task_columnar_materialize(path: Path, backend: str) -> Callable[[], int]:
    if backend == "bamboo":
        import bamboo as bm

        def run() -> int:
            table = bm.read_columns(
                str(path),
                columns=["qname", "rname", "pos", "mapq", "cigar"],
            )
            return table.num_rows

        return run

    if backend == "pysam":
        import pysam

        def run() -> int:
            qname: list[str | None] = []
            rname: list[str | None] = []
            pos: list[int | None] = []
            mapq: list[int | None] = []
            cigar: list[str] = []
            with pysam.AlignmentFile(str(path), "rb") as bam:
                for read in bam:
                    qname.append(read.query_name)
                    rname.append(read.reference_name)
                    pos.append(read.reference_start)
                    mapq.append(read.mapping_quality)
                    cigar.append(read.cigarstring)
            table = pa.table(
                {
                    "qname": qname,
                    "rname": rname,
                    "pos": pos,
                    "mapq": mapq,
                    "cigar": cigar,
                }
            )
            return table.num_rows

        return run

    if backend == "samtools":
        def run() -> int:
            raise RuntimeError("samtools has no native columnar export")

        return run

    raise ValueError(f"unsupported backend '{backend}'")


TASK_SPECS: list[tuple[str, Callable[[Path, str], Callable[[], int]]]] = [
    ("count_records", task_count_records),
    ("iterate_materialize", task_iterate_materialize),
    ("columnar_materialize", task_columnar_materialize),
    ("region_columnar", task_region_columnar),
    ("region_fetch_arrow", task_region_fetch_arrow),
    ("region_fetch_bulk", task_region_fetch_bulk),
    ("region_fetch", task_region_fetch),
    ("arrow_export", task_arrow_export),
    ("write_roundtrip", task_write_roundtrip),
]