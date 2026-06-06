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

**Early stage / proof of concept.**

This repository is the initial skeleton for the Bamboo project.

### Planned MVP (first usable release)

- High-quality BAM/CRAM reader + writer (via Rust)
- Clean Python API: `bamboo.AlignmentFile`, record iteration, header access, etc.
- Native export to Arrow (zero-copy friendly path to Polars / pandas / DuckDB)
- Basic VCF/BCF support
- Cloud object store support (S3 + GCS at minimum)
- Good documentation and examples
- Bioconda + PyPI packages
- Comprehensive test suite against real data (including CRAM, long reads, weird tags)

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

For development (once built):

```bash
pip install maturin
maturin develop
```

## Quick Vision for the API

```python
import bamboo as bm
import polars as pl

# Open a BAM (local or s3://...)
with bm.AlignmentFile("s3://my-bucket/aligned.bam") as bam:
    # Iterate like pysam, but better
    for read in bam.fetch("chr1", 1000, 2000):
        print(read.query_name, read.reference_start)

    # Or get a proper DataFrame directly
    df = bam.to_arrow().to_polars()   # or similar
    # tags become real columns, sequences as large_string, etc.

# Same for VCF
with bm.VariantFile("cohort.vcf.gz") as vcf:
    variants = vcf.to_polars()
```

The exact API will be refined with user feedback. The north star is: "it should feel like it was designed in 2025 for people who live in Polars/Jupyter/cloud environments", not "a thin wrapper over C structs".

## Contributing

We are at the very beginning. Feedback on API design, performance targets, and must-have features is extremely welcome.

## License

MIT

## Name

"Bamboo" — fast-growing, strong yet flexible, and a nice break from the `pyhts`/`pysam` naming crowd. Also: "Bamboo for your BAMs".
