#!/usr/bin/env python3
"""Render benchmark charts and an HTML report from JSON results."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patheffects as pe
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RESULTS_DIR = ROOT / "benchmarks" / "results"

PALETTE = {
    "bamboo": "#1B9E77",
    "pysam": "#7570B3",
    "samtools": "#D95F02",
}
TASK_LABELS = {
    "count_records": "Count records",
    "iterate_materialize": "Iterate + materialize",
    "columnar_materialize": "Columnar scan (recommended)",
    "region_fetch_bulk": "Indexed region fetch (bulk)",
    "region_fetch": "Indexed region fetch",
    "arrow_export": "Arrow table export",
    "write_roundtrip": "Write round-trip",
}
BACKEND_LABELS = {
    "bamboo": "Bamboo",
    "pysam": "pysam",
    "samtools": "samtools",
}


def _speedup_vs_baseline(results: list[dict[str, Any]], baseline: str = "pysam") -> dict[tuple[str, str], float]:
    baseline_medians: dict[str, float] = {}
    for row in results:
        if row["backend"] == baseline:
            baseline_medians[row["name"]] = row["timing"]["median_s"]

    speedups: dict[tuple[str, str], float] = {}
    for row in results:
        base = baseline_medians.get(row["name"])
        if not base or base <= 0:
            continue
        speedups[(row["name"], row["backend"])] = base / row["timing"]["median_s"]
    return speedups


def _latest_result_file(results_dir: Path) -> Path:
    candidates = sorted(results_dir.glob("benchmark_*.json"))
    if not candidates:
        raise FileNotFoundError(f"No benchmark JSON files found in {results_dir}")
    return candidates[-1]


def render_report(payload: dict[str, Any], source_json: Path) -> Path:
    results = payload["results"]
    if not results:
        raise ValueError("No benchmark results to visualize")

    report_dir = source_json.parent / source_json.stem
    report_dir.mkdir(parents=True, exist_ok=True)

    tasks = list(dict.fromkeys(row["name"] for row in results))
    backends = list(dict.fromkeys(row["backend"] for row in results))
    speedups = _speedup_vs_baseline(results)

    _render_throughput_chart(results, tasks, backends, report_dir / "throughput.png")
    _render_latency_chart(results, tasks, backends, report_dir / "latency.png")
    _render_speedup_chart(speedups, tasks, backends, report_dir / "speedup_vs_pysam.png")
    _render_dashboard(results, tasks, backends, speedups, payload, report_dir / "dashboard.png")
    _write_html_report(payload, report_dir, source_json.name)

    return report_dir


def _style_axes(ax: plt.Axes, title: str, xlabel: str, ylabel: str) -> None:
    ax.set_title(title, fontsize=16, fontweight="bold", pad=14, color="#102A43")
    ax.set_xlabel(xlabel, fontsize=11, color="#334E68")
    ax.set_ylabel(ylabel, fontsize=11, color="#334E68")
    ax.grid(axis="x", color="#D9E2EC", linewidth=0.8, alpha=0.9)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_color("#D9E2EC")


def _render_throughput_chart(
    results: list[dict[str, Any]],
    tasks: list[str],
    backends: list[str],
    output_path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(12, 6.5), facecolor="#F8FAFC")
    ax.set_facecolor("#FFFFFF")

    y = np.arange(len(tasks))
    height = 0.8 / max(len(backends), 1)

    for index, backend in enumerate(backends):
        values = []
        for task in tasks:
            match = next(
                (row for row in results if row["name"] == task and row["backend"] == backend),
                None,
            )
            values.append(match["throughput_records_per_s"] if match else 0.0)
        offset = (index - (len(backends) - 1) / 2) * height
        bars = ax.barh(
            y + offset,
            values,
            height=height,
            label=BACKEND_LABELS.get(backend, backend),
            color=PALETTE.get(backend, "#627D98"),
            edgecolor="white",
            linewidth=0.8,
        )
        for bar, value in zip(bars, values, strict=True):
            if value <= 0:
                continue
            ax.text(
                bar.get_width() * 1.01,
                bar.get_y() + bar.get_height() / 2,
                f"{value:,.0f}",
                va="center",
                ha="left",
                fontsize=9,
                color="#243B53",
            )

    ax.set_yticks(y)
    ax.set_yticklabels([TASK_LABELS.get(task, task) for task in tasks])
    _style_axes(ax, "Throughput (records / second)", "Records per second", "")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _render_latency_chart(
    results: list[dict[str, Any]],
    tasks: list[str],
    backends: list[str],
    output_path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(12, 6.5), facecolor="#F8FAFC")
    ax.set_facecolor("#FFFFFF")

    y = np.arange(len(tasks))
    height = 0.8 / max(len(backends), 1)

    for index, backend in enumerate(backends):
        values = []
        for task in tasks:
            match = next(
                (row for row in results if row["name"] == task and row["backend"] == backend),
                None,
            )
            values.append(match["timing"]["median_s"] if match else 0.0)
        offset = (index - (len(backends) - 1) / 2) * height
        ax.barh(
            y + offset,
            values,
            height=height,
            label=BACKEND_LABELS.get(backend, backend),
            color=PALETTE.get(backend, "#627D98"),
            edgecolor="white",
            linewidth=0.8,
        )

    ax.set_yticks(y)
    ax.set_yticklabels([TASK_LABELS.get(task, task) for task in tasks])
    _style_axes(ax, "Median latency (seconds)", "Seconds (lower is better)", "")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _render_speedup_chart(
    speedups: dict[tuple[str, str], float],
    tasks: list[str],
    backends: list[str],
    output_path: Path,
) -> None:
    if not speedups:
        fig, ax = plt.subplots(figsize=(12, 3), facecolor="#F8FAFC")
        ax.axis("off")
        ax.text(
            0.5,
            0.5,
            "Install pysam to compute speedup vs baseline",
            ha="center",
            va="center",
            fontsize=14,
            color="#627D98",
        )
        fig.savefig(output_path, dpi=180, bbox_inches="tight")
        plt.close(fig)
        return

    fig, ax = plt.subplots(figsize=(12, 6.5), facecolor="#F8FAFC")
    ax.set_facecolor("#FFFFFF")

    y = np.arange(len(tasks))
    height = 0.8 / max(len(backends), 1)

    for index, backend in enumerate(backends):
        values = [speedups.get((task, backend), 0.0) for task in tasks]
        offset = (index - (len(backends) - 1) / 2) * height
        bars = ax.barh(
            y + offset,
            values,
            height=height,
            label=BACKEND_LABELS.get(backend, backend),
            color=PALETTE.get(backend, "#627D98"),
            edgecolor="white",
            linewidth=0.8,
        )
        for bar, value in zip(bars, values, strict=True):
            if value <= 0:
                continue
            ax.text(
                bar.get_width() + 0.03,
                bar.get_y() + bar.get_height() / 2,
                f"{value:.2f}x",
                va="center",
                ha="left",
                fontsize=9,
                color="#243B53",
            )

    ax.axvline(1.0, color="#9FB3C8", linestyle="--", linewidth=1.2)
    ax.set_yticks(y)
    ax.set_yticklabels([TASK_LABELS.get(task, task) for task in tasks])
    _style_axes(ax, "Speedup vs pysam (median latency)", "× faster than pysam", "")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _render_dashboard(
    results: list[dict[str, Any]],
    tasks: list[str],
    backends: list[str],
    speedups: dict[tuple[str, str], float],
    payload: dict[str, Any],
    output_path: Path,
) -> None:
    fig = plt.figure(figsize=(16, 10), facecolor="#0B1F33")
    gs = fig.add_gridspec(2, 2, width_ratios=[1.2, 1.0], height_ratios=[0.9, 1.1], wspace=0.25, hspace=0.3)

    ax_title = fig.add_subplot(gs[0, :])
    ax_title.axis("off")
    title = "Bamboo BAM Reader Benchmarks"
    subtitle = (
        f"{payload['config']['records']:,} synthetic reads · "
        f"{payload['platform']['system']} {payload['platform']['machine']} · "
        f"generated {payload['generated_at'][:19]} UTC"
    )
    ax_title.text(0.03, 0.72, title, fontsize=28, fontweight="bold", color="#F0F4F8")
    ax_title.text(0.03, 0.34, subtitle, fontsize=13, color="#9FB3C8")

    bamboo_wins = 0
    comparisons = 0
    for task in tasks:
        bamboo = next((row for row in results if row["name"] == task and row["backend"] == "bamboo"), None)
        pysam_row = next((row for row in results if row["name"] == task and row["backend"] == "pysam"), None)
        if bamboo and pysam_row:
            comparisons += 1
            if bamboo["timing"]["median_s"] < pysam_row["timing"]["median_s"]:
                bamboo_wins += 1

    card = FancyBboxPatch(
        (0.62, 0.18),
        0.34,
        0.62,
        boxstyle="round,pad=0.02,rounding_size=0.02",
        linewidth=0,
        facecolor="#16324F",
    )
    ax_title.add_patch(card)
    ax_title.text(0.64, 0.66, "Headline", fontsize=12, color="#9FB3C8")
    if comparisons:
        ax_title.text(
            0.64,
            0.48,
            f"Bamboo faster on\n{bamboo_wins}/{comparisons} tasks",
            fontsize=22,
            fontweight="bold",
            color="#63E6BE",
        )
    else:
        ax_title.text(0.64, 0.48, "Install pysam\nfor comparisons", fontsize=18, color="#F0F4F8")

    ax_throughput = fig.add_subplot(gs[1, 0])
    ax_throughput.set_facecolor("#102A43")
    y = np.arange(len(tasks))
    for index, backend in enumerate(backends):
        values = []
        for task in tasks:
            match = next(
                (row for row in results if row["name"] == task and row["backend"] == backend),
                None,
            )
            values.append(match["throughput_records_per_s"] if match else 0.0)
        ax_throughput.plot(
            values,
            y + index * 0.08,
            marker="o",
            linewidth=2.4,
            markersize=7,
            label=BACKEND_LABELS.get(backend, backend),
            color=PALETTE.get(backend, "#627D98"),
        )
    ax_throughput.set_yticks(y)
    ax_throughput.set_yticklabels([TASK_LABELS.get(task, task) for task in tasks], color="#D9E2EC")
    ax_throughput.set_xlabel("Records / second", color="#D9E2EC")
    ax_throughput.set_title("Throughput profile", color="#F0F4F8", fontsize=14, fontweight="bold")
    ax_throughput.tick_params(colors="#D9E2EC")
    ax_throughput.grid(color="#243B53", alpha=0.6)
    for spine in ax_throughput.spines.values():
        spine.set_color("#243B53")
    ax_throughput.legend(frameon=False, loc="lower right", labelcolor="#D9E2EC")

    ax_heatmap = fig.add_subplot(gs[1, 1])
    heatmap_backends = [b for b in backends if any((task, b) in speedups for task in tasks)]
    if heatmap_backends:
        matrix = np.array(
            [[speedups.get((task, backend), np.nan) for backend in heatmap_backends] for task in tasks],
            dtype=float,
        )
        im = ax_heatmap.imshow(matrix, aspect="auto", cmap="viridis", vmin=0.5, vmax=max(2.0, np.nanmax(matrix)))
        ax_heatmap.set_xticks(np.arange(len(heatmap_backends)))
        ax_heatmap.set_xticklabels([BACKEND_LABELS.get(b, b) for b in heatmap_backends], color="#D9E2EC")
        ax_heatmap.set_yticks(np.arange(len(tasks)))
        ax_heatmap.set_yticklabels([TASK_LABELS.get(task, task) for task in tasks], color="#D9E2EC")
        ax_heatmap.set_title("Speedup vs pysam", color="#F0F4F8", fontsize=14, fontweight="bold")
        for i, task in enumerate(tasks):
            for j, backend in enumerate(heatmap_backends):
                value = speedups.get((task, backend))
                if value is None:
                    continue
                text = ax_heatmap.text(
                    j,
                    i,
                    f"{value:.2f}x",
                    ha="center",
                    va="center",
                    color="white" if value < 1.4 else "#102A43",
                    fontsize=10,
                    fontweight="bold",
                )
                text.set_path_effects([pe.withStroke(linewidth=2, foreground="#102A43")])
        cbar = fig.colorbar(im, ax=ax_heatmap, fraction=0.046, pad=0.04)
        cbar.ax.tick_params(colors="#D9E2EC")
    else:
        ax_heatmap.axis("off")
        ax_heatmap.text(0.5, 0.5, "No pysam baseline available", ha="center", va="center", color="#D9E2EC")

    fig.savefig(output_path, dpi=180, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


def _write_html_report(payload: dict[str, Any], report_dir: Path, source_name: str) -> None:
    rows = payload["results"]
    table_rows = []
    for row in rows:
        table_rows.append(
            "<tr>"
            f"<td>{TASK_LABELS.get(row['name'], row['name'])}</td>"
            f"<td>{BACKEND_LABELS.get(row['backend'], row['backend'])}</td>"
            f"<td>{row['record_count']:,}</td>"
            f"<td>{row['timing']['median_s']:.4f}</td>"
            f"<td>{row['throughput_records_per_s']:,.0f}</td>"
            "</tr>"
        )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>Bamboo Benchmark Report</title>
  <style>
    :root {{
      --bg: #f8fafc;
      --card: #ffffff;
      --ink: #102a43;
      --muted: #627d98;
      --accent: #1b9e77;
      --border: #d9e2ec;
    }}
    body {{
      margin: 0;
      font-family: "SF Pro Text", "Segoe UI", sans-serif;
      background: linear-gradient(180deg, #f8fafc 0%, #e6eff7 100%);
      color: var(--ink);
    }}
    .wrap {{
      max-width: 1200px;
      margin: 0 auto;
      padding: 32px 24px 64px;
    }}
    h1 {{
      margin: 0 0 8px;
      font-size: 2.2rem;
    }}
    .meta {{
      color: var(--muted);
      margin-bottom: 24px;
    }}
    .grid {{
      display: grid;
      grid-template-columns: 1fr;
      gap: 20px;
    }}
    .card {{
      background: var(--card);
      border: 1px solid var(--border);
      border-radius: 16px;
      padding: 18px;
      box-shadow: 0 10px 30px rgba(16, 42, 67, 0.06);
    }}
    img {{
      width: 100%;
      border-radius: 12px;
      display: block;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.95rem;
    }}
    th, td {{
      text-align: left;
      padding: 10px 12px;
      border-bottom: 1px solid var(--border);
    }}
    th {{
      color: var(--muted);
      font-weight: 600;
      font-size: 0.8rem;
      text-transform: uppercase;
      letter-spacing: 0.04em;
    }}
    .pill {{
      display: inline-block;
      background: rgba(27, 158, 119, 0.12);
      color: var(--accent);
      padding: 4px 10px;
      border-radius: 999px;
      font-size: 0.85rem;
      font-weight: 600;
    }}
  </style>
</head>
<body>
  <div class="wrap">
    <span class="pill">Bamboo benchmark harness</span>
    <h1>BAM reader competitor report</h1>
    <p class="meta">
      {payload["config"]["records"]:,} reads ·
      {payload["platform"]["system"]} {payload["platform"]["machine"]} ·
      source <code>{source_name}</code>
    </p>
    <div class="grid">
      <div class="card">
        <img src="dashboard.png" alt="Benchmark dashboard" />
      </div>
      <div class="card">
        <img src="throughput.png" alt="Throughput chart" />
      </div>
      <div class="card">
        <img src="latency.png" alt="Latency chart" />
      </div>
      <div class="card">
        <img src="speedup_vs_pysam.png" alt="Speedup vs pysam chart" />
      </div>
      <div class="card">
        <table>
          <thead>
            <tr>
              <th>Task</th>
              <th>Backend</th>
              <th>Records</th>
              <th>Median (s)</th>
              <th>Throughput</th>
            </tr>
          </thead>
          <tbody>
            {"".join(table_rows)}
          </tbody>
        </table>
      </div>
    </div>
  </div>
</body>
</html>
"""
    (report_dir / "index.html").write_text(html)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "json_path",
        nargs="?",
        type=Path,
        default=None,
        help="benchmark JSON file (default: latest in benchmarks/results)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    json_path = args.json_path or _latest_result_file(DEFAULT_RESULTS_DIR)
    payload = json.loads(json_path.read_text())
    report_dir = render_report(payload, json_path)
    print(f"Wrote report to {report_dir / 'index.html'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())