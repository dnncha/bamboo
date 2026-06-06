# Bamboo

**High-performance, modern Python access to high-throughput sequencing data.**

Bamboo is the spiritual successor to `pysam`. It provides a fast, ergonomic, and data-science-native Python interface to the core HTS formats (BAM, SAM, CRAM, VCF, BCF, and related index formats).

## Why Bamboo?

- **pysam is showing its age**: The Cython wrapper around htslib is powerful but has ergonomic, performance, and integration limitations in 2026-era Python data workflows.
- **The world moved on**: We now have excellent pure-Rust parsers (`noodles`), Arrow as the lingua franca of data, Polars, DuckDB, cloud object stores as first-class citizens, and a strong desire for zero-copy and streaming access.
- **Real production needs**: Large cohorts, long-read data, cloud-native pipelines, single-cell scale, and tight integration with the rest of the scientific Python stack.

Bamboo aims to be the library you actually want to use when writing modern genomics code in Python.

## Goals

- **Blazing fast** — Rust core (built on or alongside `noodles` where it makes sense).
- **Pythonic & delightful** — Modern API design, great error messages, context managers, iterators that feel native.
- **Data-native** — First-class Arrow support. `bam.records()` should be able to give you a Polars DataFrame or pyarrow Table with tags as columns with minimal friction.
- **Cloud first** — Excellent support for reading directly from S3, GCS, Azure Blob, with smart caching and range requests.
- **Safe & correct** — Memory safety from Rust + extensive testing against real-world data and edge cases.
- **Interoperable** — Play nicely with existing ecosystems (pysam compatibility shims where helpful, but not at the cost of a better design).
- **Minimal dependencies** for the core path.

## Current Status

**MVP milestone: BAM reading works.**

Implemented today:
- Rust workspace (`bamboo-core`, `bamboo-noodles`, `bamboo-py`)
- `bamboo.AlignmentFile` with iteration, region fetch, and indexed fetch (`.bai`)
- Columnar scan to PyArrow via `read_bam_table()` and `AlignmentFile.to_arrow()`
- Polars helper via `bamboo.to_polars()`
- Test fixtures and examples

Still planned for full MVP:
- BAM/CRAM writing and CRAM decoding (via `bamboo-htslib`)
- Basic VCF/BCF support
- Cloud object store support (S3 + GCS at minimum)
- PyPI / Bioconda packaging
- Comprehensive test suite against real-world data (long reads, weird tags)

Later:
- Full pysam-like API surface + compatibility layer
- Writing support for more formats
- Indexed random access (bai/crai/csi/tabix)
- Pileup / depth engines
- Variant annotation helpers
- WASM / browser support (stretch)

## Installation (future)

```bash
pip install bamboo
# or
conda install -c bioconda bamboo
```

For development:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
pip install maturin pyarrow pytest
maturin develop
pytest
```

Generate test fixtures:

```bash
cargo run -p bamboo-noodles --example generate_fixtures
```

## Quick Start

```python
import bamboo as bm

with bm.AlignmentFile("aligned.bam") as bam:
    print(bam.header())          # {'chr1': 248956422, ...}
    print(bam.references())      # ['chr1', 'chr2', ...]

    for read in bam.fetch(region="chr1:1000000-1001000"):
        print(read.query_name, read.reference_start, read.cigarstring)

    table = bam.to_arrow(columns=["qname", "rname", "pos", "mapq"], region="chr1:1-1000000")
    df = bm.to_polars(table)     # requires polars
```

See `examples/read_bam.py` for a command-line demo.

## Quick Vision for the API

```python
import bamboo as bm

# Open a BAM (local today; s3:// planned)
with bm.AlignmentFile("aligned.bam") as bam:
    for read in bam.fetch(region="chr1:1000-2000"):
        print(read.query_name, read.reference_start)

    df = bm.to_polars(bam.to_arrow())

# VCF support is planned
# with bm.VariantFile("cohort.vcf.gz") as vcf:
#     variants = vcf.to_polars()
```

The exact API will be refined with user feedback. The north star is: "it should feel like it was designed in 2025 for people who live in Polars/Jupyter/cloud environments", not "a thin wrapper over C structs".

## Contributing

We are at the very beginning. Feedback on API design, performance targets, and must-have features is extremely welcome.

## License

MIT

## Name

"Bamboo" — fast-growing, strong yet flexible, and a nice break from the `pyhts`/`pysam` naming crowd. Also: "Bamboo for your BAMs".
