#!/usr/bin/env python3
"""Run Bamboo vs competitor benchmarks and write JSON results."""

from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from benchmarks.harness import run_benchmark
from benchmarks.tasks import (
    BAM_TASK_SPECS,
    CRAM_TASK_SPECS,
    TASK_SPECS,
    available_backends,
)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = ROOT / "benchmarks" / "data"
DEFAULT_RESULTS_DIR = ROOT / "benchmarks" / "results"


def ensure_bench_data(data_dir: Path, record_count: int) -> tuple[Path, Path, Path]:
    bam_path = data_dir / f"bench_{record_count}.bam"
    bai_path = bam_path.with_suffix(".bam.bai")
    cram_path = data_dir / f"bench_{record_count}.cram"
    crai_path = Path(f"{cram_path}.crai")
    fasta_path = data_dir / f"bench_{record_count}.fasta"
    if (
        bam_path.exists()
        and bai_path.exists()
        and cram_path.exists()
        and crai_path.exists()
        and fasta_path.exists()
    ):
        return bam_path, cram_path, fasta_path

    data_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "cargo",
        "run",
        "-p",
        "bamboo-noodles",
        "--example",
        "generate_bench_data",
        "--",
        str(data_dir),
        str(record_count),
    ]
    print(f"Generating benchmark BAM/CRAM ({record_count:,} records)...")
    subprocess.run(cmd, cwd=ROOT, check=True)
    return bam_path, cram_path, fasta_path


def package_versions() -> dict[str, str | None]:
    versions: dict[str, str | None] = {"python": platform.python_version()}
    for package in ("bamboo", "pysam", "pyarrow"):
        try:
            module = __import__(package)
            versions[package] = getattr(module, "__version__", "unknown")
        except ImportError:
            versions[package] = None

    versions["samtools"] = None
    try:
        completed = subprocess.run(
            ["samtools", "--version"],
            check=True,
            capture_output=True,
            text=True,
        )
        versions["samtools"] = completed.stdout.splitlines()[0]
    except (FileNotFoundError, subprocess.CalledProcessError):
        pass
    return versions


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--records",
        type=int,
        default=100_000,
        help="number of synthetic reads in the benchmark BAM (default: 100000)",
    )
    parser.add_argument(
        "--rounds",
        type=int,
        default=5,
        help="timed iterations per task/backend (default: 5)",
    )
    parser.add_argument(
        "--warmup",
        type=int,
        default=2,
        help="warmup iterations per task/backend (default: 2)",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=DEFAULT_DATA_DIR,
        help="directory for generated benchmark BAMs",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="JSON output path (default: benchmarks/results/<timestamp>.json)",
    )
    parser.add_argument(
        "--backends",
        nargs="*",
        default=None,
        help="subset of backends to run (default: all available)",
    )
    parser.add_argument(
        "--tasks",
        nargs="*",
        default=None,
        help="subset of task names to run (default: all)",
    )
    parser.add_argument(
        "--visualize",
        action="store_true",
        help="render charts and HTML report after the run",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    backends = args.backends or available_backends()
    if not backends:
        print("No backends available. Install bamboo and at least one competitor.", file=sys.stderr)
        return 1

    wanted_tasks = set(args.tasks) if args.tasks else None
    bam_task_specs = BAM_TASK_SPECS
    cram_task_specs = CRAM_TASK_SPECS
    if wanted_tasks:
        bam_task_specs = [(name, fn) for name, fn in BAM_TASK_SPECS if name in wanted_tasks]
        cram_task_specs = [(name, fn) for name, fn in CRAM_TASK_SPECS if name in wanted_tasks]
        if not bam_task_specs and not cram_task_specs:
            print(f"No matching tasks for: {', '.join(sorted(wanted_tasks))}", file=sys.stderr)
            return 1

    bam_path, cram_path, fasta_path = ensure_bench_data(args.data_dir, args.records)
    results_dir = DEFAULT_RESULTS_DIR
    results_dir.mkdir(parents=True, exist_ok=True)
    output_path = args.output or results_dir / f"benchmark_{args.records}_{datetime.now(timezone.utc):%Y%m%dT%H%M%SZ}.json"

    selected_tasks = [name for name, _ in bam_task_specs] + [name for name, _ in cram_task_specs]
    print(f"Benchmark BAM: {bam_path}")
    print(f"Benchmark CRAM: {cram_path}")
    print(f"Backends: {', '.join(backends)}")
    print(f"Tasks: {', '.join(selected_tasks)}")
    print()

    run_results: list[dict] = []
    expected_counts: dict[str, int] = {}

    for task_name, task_fn in bam_task_specs:
        print(f"== {task_name} ==")
        for backend in backends:
            if task_name in {"arrow_export", "columnar_materialize"} and backend == "samtools":
                print(f"  - samtools: skipped (no {task_name})")
                continue
            if task_name == "write_roundtrip" and backend == "samtools":
                print("  - samtools: skipped (write not benchmarked)")
                continue

            try:
                fn = task_fn(bam_path, backend)
                result = run_benchmark(
                    name=task_name,
                    backend=backend,
                    fn=fn,
                    warmup=args.warmup,
                    rounds=args.rounds,
                    metadata={"bam_path": str(bam_path), "region": "chr1:1000000-5000000"},
                )
            except Exception as exc:  # noqa: BLE001 - benchmark harness reports failures
                print(f"  - {backend}: FAILED ({exc})")
                continue

            if task_name in expected_counts:
                if result.record_count != expected_counts[task_name]:
                    print(
                        f"  - {backend}: WARNING count mismatch "
                        f"({result.record_count} != {expected_counts[task_name]})"
                    )
            else:
                expected_counts[task_name] = result.record_count

            print(
                f"  - {backend}: {result.record_count:,} records, "
                f"{result.throughput_records_per_s:,.0f} rec/s "
                f"(median {result.timing.median_s:.3f}s)"
            )
            run_results.append(result.to_dict())
        print()

    for task_name, task_fn in cram_task_specs:
        print(f"== {task_name} ==")
        for backend in backends:
            if backend == "samtools":
                print(f"  - samtools: skipped (no {task_name})")
                continue

            try:
                fn = task_fn(
                    cram_path,
                    backend,
                    reference_path=fasta_path,
                )
                result = run_benchmark(
                    name=task_name,
                    backend=backend,
                    fn=fn,
                    warmup=args.warmup,
                    rounds=args.rounds,
                    metadata={
                        "cram_path": str(cram_path),
                        "fasta_path": str(fasta_path),
                        "region": "chr1:1000000-5000000",
                    },
                )
            except Exception as exc:  # noqa: BLE001 - benchmark harness reports failures
                print(f"  - {backend}: FAILED ({exc})")
                continue

            if task_name in expected_counts:
                if result.record_count != expected_counts[task_name]:
                    print(
                        f"  - {backend}: WARNING count mismatch "
                        f"({result.record_count} != {expected_counts[task_name]})"
                    )
            else:
                expected_counts[task_name] = result.record_count

            print(
                f"  - {backend}: {result.record_count:,} records, "
                f"{result.throughput_records_per_s:,.0f} rec/s "
                f"(median {result.timing.median_s:.3f}s)"
            )
            run_results.append(result.to_dict())
        print()

    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
            "processor": platform.processor(),
        },
        "versions": package_versions(),
        "config": {
            "records": args.records,
            "rounds": args.rounds,
            "warmup": args.warmup,
            "bam_path": str(bam_path),
            "backends": backends,
            "tasks": selected_tasks,
        },
        "results": run_results,
    }

    output_path.write_text(json.dumps(payload, indent=2))
    print(f"Wrote results to {output_path}")

    if args.visualize:
        from benchmarks.visualize import render_report

        report_dir = render_report(payload, output_path)
        print(f"Wrote report to {report_dir}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())