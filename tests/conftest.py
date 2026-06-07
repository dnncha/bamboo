"""Shared pytest fixtures."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

DATA_DIR = Path(__file__).parent / "data"
BENCH_DATA_DIR = Path(__file__).resolve().parents[1] / "benchmarks" / "data"
TINY_BAM = DATA_DIR / "tiny.bam"
TINY_BAI = DATA_DIR / "tiny.bam.bai"
BENCH_RECORD_COUNT = 50_000


def _ensure_bench_assets(record_count: int = BENCH_RECORD_COUNT) -> Path:
    """Generate synthetic benchmark BAM/CRAM if missing (realistic parity gate)."""
    bam_path = BENCH_DATA_DIR / f"bench_{record_count}.bam"
    bai_path = bam_path.with_suffix(".bam.bai")
    cram_path = BENCH_DATA_DIR / f"bench_{record_count}.cram"
    crai_path = Path(f"{cram_path}.crai")
    fasta_path = BENCH_DATA_DIR / f"bench_{record_count}.fasta"
    fai_path = Path(f"{fasta_path}.fai")

    if all(
        path.exists()
        for path in (bam_path, bai_path, cram_path, crai_path, fasta_path, fai_path)
    ):
        return bam_path

    repo_root = Path(__file__).resolve().parents[1]
    subprocess.run(
        [
            "cargo",
            "run",
            "-p",
            "bamboo-noodles",
            "--example",
            "generate_bench_data",
            "--",
            str(BENCH_DATA_DIR),
            str(record_count),
        ],
        cwd=repo_root,
        check=True,
    )
    return bam_path


@pytest.fixture(scope="session")
def tiny_bam_path() -> Path:
    if not TINY_BAM.exists():
        pytest.skip(
            "missing tests/data/tiny.bam — run "
            "`cargo run -p bamboo-noodles --example generate_fixtures`"
        )
    return TINY_BAM


@pytest.fixture(scope="session")
def tiny_bam_with_index(tiny_bam_path: Path) -> Path:
    if not TINY_BAI.exists():
        pytest.skip(
            "missing tests/data/tiny.bam.bai — run "
            "`cargo run -p bamboo-noodles --example generate_fixtures`"
        )
    return tiny_bam_path


@pytest.fixture(scope="session")
def bench_bam_with_index() -> Path:
    """50k-read synthetic BAM + BAI for real-world pysam parity tests."""
    return _ensure_bench_assets(BENCH_RECORD_COUNT)


@pytest.fixture(scope="session")
def bench_cram_with_index(bench_bam_with_index: Path) -> Path:
    _ = bench_bam_with_index
    return BENCH_DATA_DIR / f"bench_{BENCH_RECORD_COUNT}.cram"


@pytest.fixture(scope="session")
def bench_fasta_with_index(bench_bam_with_index: Path) -> Path:
    _ = bench_bam_with_index
    return BENCH_DATA_DIR / f"bench_{BENCH_RECORD_COUNT}.fasta"