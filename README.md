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

**MVP: BAM + CRAM reading, columnar Arrow export, pysam parity CI, pileup (htslib).**

Implemented today:
- Rust workspace (`bamboo-core`, `bamboo-io`, `bamboo-noodles`, `bamboo-htslib`, `bamboo-py`)
- `bamboo.AlignmentFile` — iteration, region fetch, indexed fetch, BAM write
- `bamboo.CramFile` — CRAM decode with external reference, columnar + pileup
- VCF/BCF readers with Arrow export
- Columnar scan to PyArrow via `read_columns()` / `to_arrow()` (BAM + CRAM)
- pysam parity tests on tiny fixtures **and** 50k-read synthetic cohort data
- `from bamboo.compat import pysam` drop-in shim — see **[MIGRATION.md](MIGRATION.md)**
- Cloud I/O: `s3://`, `gs://`, `https://`, `file://`
- Pileup via htslib (optional build feature)

Still planned:
- PyPI publish (`bamboo-hts`) + Bioconda merge
- Unified `AlignmentFile(..., "rc")` for CRAM
- SAM/CRAM writing, CSI/tabix, coverage APIs
- Parity on messy production files (long reads, exotic tags)

## Installation

```bash
pip install bamboo-hts
```

```python
import bamboo  # Python module name (not bamboo-hts)
```

> **Note:** PyPI package is `bamboo-hts` because [`bamboo`](https://pypi.org/project/bamboo/) is an unrelated imaging library.

Optional extras:

```bash
pip install bamboo-hts[polars]   # Polars adapter
pip install bamboo-hts[pandas]   # Pandas adapter
```

**Conda** (after [Bioconda recipe](conda/bioconda/meta.yaml) is merged):

```bash
conda install -c bioconda -c conda-forge bamboo-hts
```

Wheels bundle htslib — no system `libhts` or maturin required. See [PACKAGING.md](PACKAGING.md) for release and conda-build details.

### Development

```bash
python3.12 -m venv .venv
source .venv/bin/activate
pip install maturin pyarrow pytest 'pysam>=0.22'
maturin develop --release --features htslib
pytest
```

Smoke-test a release wheel locally: `./scripts/verify_wheel.sh`

Generate test fixtures:

```bash
cargo run -p bamboo-noodles --example generate_fixtures
```

## Quick Start

**Migrating from pysam?** Start with [MIGRATION.md](MIGRATION.md) — one-line import swap, API table, validation checklist.

**Killer workflow** (indexed region → Arrow → Polars QC):

```python
import bamboo as bm

table = bm.read_columns(
    "cohort.bam",
    columns=["qname", "rname", "pos", "mapq", "flag"],
    region="chr1:1000000-5000000",
    min_mapq=30,
)
df = bm.to_polars(table)  # requires polars
print(df.group_by("rname").len())
```

Runnable demo: `python examples/cohort_region_qc.py tests/data/tiny.bam --region chr1:100-500`

**Record iteration** (pysam-familiar):

```python
import bamboo as bm

with bm.AlignmentFile("aligned.bam") as bam:
    for read in bam.fetch(region="chr1:1000000-1001000"):
        print(read.query_name, read.reference_start, read.cigarstring)
```

See also `examples/read_bam.py`.

## Cloud and remote paths

Bamboo reads BAMs (and sidecar `.bai` indexes when present) from local paths and cloud URIs through the same API:

```python
import bamboo as bm

# Local path or file:// URI
with bm.AlignmentFile("aligned.bam") as bam:
    ...

# S3 (uses default AWS credential chain: env vars, ~/.aws, IAM role, etc.)
with bm.AlignmentFile("s3://my-bucket/cohort/sample.bam") as bam:
    ...

# GCS (uses Application Default Credentials / GOOGLE_APPLICATION_CREDENTIALS)
with bm.AlignmentFile("gs://my-bucket/cohort/sample.bam") as bam:
    ...

# HTTPS (public or pre-signed URLs)
with bm.AlignmentFile("https://example.com/public/sample.bam") as bam:
    ...
```

Index discovery tries `sample.bam.bai` then `sample.bai` next to the BAM URI. For indexed `fetch()`, the `.bai` must be reachable at one of those locations.

## Writing BAMs

BAM writing uses local paths today (`wb` / `w` mode). Copy reads from an existing file with a pysam-style template header:

```python
import bamboo as bm

with bm.AlignmentFile("input.bam") as src:
    with bm.AlignmentFile("output.bam", "wb", template=src) as out:
        for read in src:
            out.write(read)
```

Or supply a reference dictionary when creating a new file:

```python
with bm.AlignmentFile("output.bam", "wb", header={"chr1": 248956422}) as out:
    out.write(read)
```

## Quick Vision for the API

```python
import bamboo as bm

# Open a BAM from local disk or a cloud URI
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
