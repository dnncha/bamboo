"""Shared pytest fixtures."""

from __future__ import annotations

from pathlib import Path

import pytest

DATA_DIR = Path(__file__).parent / "data"
TINY_BAM = DATA_DIR / "tiny.bam"
TINY_BAI = DATA_DIR / "tiny.bam.bai"


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