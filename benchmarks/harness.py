"""Timing utilities for repeatable benchmark runs."""

from __future__ import annotations

import gc
import statistics
import time
from dataclasses import asdict, dataclass
from typing import Any, Callable


@dataclass(frozen=True)
class TimingStats:
    median_s: float
    mean_s: float
    stdev_s: float
    min_s: float
    max_s: float
    rounds: int

    def to_dict(self) -> dict[str, float | int]:
        return asdict(self)


@dataclass(frozen=True)
class BenchmarkResult:
    name: str
    backend: str
    record_count: int
    timing: TimingStats
    throughput_records_per_s: float
    metadata: dict[str, Any]

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "backend": self.backend,
            "record_count": self.record_count,
            "timing": self.timing.to_dict(),
            "throughput_records_per_s": self.throughput_records_per_s,
            "metadata": self.metadata,
        }


def time_callable(
    fn: Callable[[], int],
    *,
    warmup: int = 2,
    rounds: int = 5,
) -> tuple[TimingStats, int]:
    """Run `fn` with warmup and return timing stats plus the last record count."""
    gc.collect()

    last_count = 0
    for _ in range(warmup):
        last_count = fn()
        gc.collect()

    samples: list[float] = []
    for _ in range(rounds):
        gc.collect()
        start = time.perf_counter()
        last_count = fn()
        samples.append(time.perf_counter() - start)

    stats = TimingStats(
        median_s=statistics.median(samples),
        mean_s=statistics.mean(samples),
        stdev_s=statistics.pstdev(samples) if len(samples) > 1 else 0.0,
        min_s=min(samples),
        max_s=max(samples),
        rounds=rounds,
    )
    return stats, last_count


def run_benchmark(
    *,
    name: str,
    backend: str,
    fn: Callable[[], int],
    warmup: int = 2,
    rounds: int = 5,
    metadata: dict[str, Any] | None = None,
) -> BenchmarkResult:
    timing, record_count = time_callable(fn, warmup=warmup, rounds=rounds)
    throughput = record_count / timing.median_s if timing.median_s > 0 else 0.0
    return BenchmarkResult(
        name=name,
        backend=backend,
        record_count=record_count,
        timing=timing,
        throughput_records_per_s=throughput,
        metadata=metadata or {},
    )