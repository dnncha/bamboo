"""Smoke tests for the benchmark harness."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from benchmarks.harness import run_benchmark
from benchmarks.tasks import TASK_SPECS, available_backends, task_count_records


def test_available_backends_includes_bamboo() -> None:
    assert "bamboo" in available_backends()


def test_count_records_tiny_fixture(tiny_bam_path: Path) -> None:
    pytest.importorskip("bamboo")
    result = run_benchmark(
        name="count_records",
        backend="bamboo",
        fn=task_count_records(tiny_bam_path, "bamboo"),
        warmup=0,
        rounds=1,
    )
    assert result.record_count == 2
    assert result.timing.median_s >= 0


def test_task_registry_is_non_empty() -> None:
    assert TASK_SPECS
    names = {name for name, _ in TASK_SPECS}
    assert "count_records" in names
    assert "arrow_export" in names


def test_visualize_accepts_minimal_payload(tmp_path: Path) -> None:
    pytest.importorskip("matplotlib")
    from benchmarks.visualize import render_report

    payload = {
        "generated_at": "2026-06-06T00:00:00+00:00",
        "platform": {"system": "test", "machine": "test", "release": "", "processor": ""},
        "config": {"records": 2},
        "results": [
            {
                "name": "count_records",
                "backend": "bamboo",
                "record_count": 2,
                "timing": {
                    "median_s": 0.01,
                    "mean_s": 0.01,
                    "stdev_s": 0.0,
                    "min_s": 0.01,
                    "max_s": 0.01,
                    "rounds": 1,
                },
                "throughput_records_per_s": 200.0,
                "metadata": {},
            },
            {
                "name": "count_records",
                "backend": "pysam",
                "record_count": 2,
                "timing": {
                    "median_s": 0.02,
                    "mean_s": 0.02,
                    "stdev_s": 0.0,
                    "min_s": 0.02,
                    "max_s": 0.02,
                    "rounds": 1,
                },
                "throughput_records_per_s": 100.0,
                "metadata": {},
            },
        ],
    }
    source = tmp_path / "sample.json"
    source.write_text(json.dumps(payload))
    report_dir = render_report(payload, source)
    assert (report_dir / "index.html").exists()
    assert (report_dir / "dashboard.png").exists()