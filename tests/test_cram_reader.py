"""Integration tests for the Bamboo CRAM reader."""

from __future__ import annotations

from pathlib import Path

import pytest

bamboo = pytest.importorskip("bamboo")


@pytest.fixture(scope="session")
def tiny_cram_path() -> Path:
    path = Path(__file__).parent / "data" / "tiny.cram"
    if not path.exists():
        pytest.skip(
            "missing tests/data/tiny.cram — run "
            "`cargo run -p bamboo-noodles --example generate_fixtures`"
        )
    return path


def test_cram_file_count_and_references(tiny_cram_path: Path) -> None:
    with bamboo.CramFile(str(tiny_cram_path)) as cram:
        assert cram.count() == 2
        assert cram.references() == ["chr1", "chr2"]
        assert cram.reference_lengths == [1000, 1000]
        assert cram.header() == {"chr1": 1000, "chr2": 1000}