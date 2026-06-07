# Migrating from pysam to Bamboo

Bamboo is designed as a **drop-in replacement for common pysam read paths**, with a faster columnar route for analytics workloads. Start with a one-line import change, validate parity on your data, then adopt columnar APIs where you iterate millions of records.

## Install

```bash
pip install bamboo-hts
import bamboo  # module name unchanged
# conda (after bioconda merge): conda install -c bioconda -c conda-forge bamboo-hts
```

## Quick swap (record iteration)

```python
# before
import pysam

with pysam.AlignmentFile("sample.bam", "rb") as bam:
    for read in bam.fetch("chr1", 1_000_000, 5_000_000):
        print(read.query_name, read.reference_start, read.flag)

# after — compat shim, same call patterns
from bamboo.compat import pysam

with pysam.AlignmentFile("sample.bam") as bam:
    for read in bam.fetch("chr1", 1_000_000, 5_000_000):
        print(read.query_name, read.reference_start, read.flag)
```

Or use Bamboo directly (recommended once validated):

```python
import bamboo as bm

with bm.AlignmentFile("sample.bam") as bam:
    for read in bam.fetch(region="chr1:1000001-5000000"):
        print(read.query_name, read.reference_start, read.flag)
```

## API mapping

| pysam | Bamboo | Notes |
|-------|--------|-------|
| `pysam.AlignmentFile(path, "rb")` | `bamboo.AlignmentFile(path)` | Default mode is read-binary |
| `pysam.AlignmentFile(path, "wb", header=...)` | `bamboo.AlignmentFile(path, "wb", template=...)` | Local paths only for write |
| `pysam.AlignmentFile(path, "rc", reference_filename=...)` | `bamboo.CramFile(path, reference_filename=...)` | CRAM uses a dedicated class today |
| `AlignedSegment` / `AlignedRead` | `bamboo.AlignedSegment` | Same core getters |
| `read.query_name` | `read.query_name` | |
| `read.reference_name` | `read.reference_name` | |
| `read.reference_start` | `read.reference_start` | 0-based, same as pysam |
| `read.mapping_quality` | `read.mapping_quality` | |
| `read.cigarstring` | `read.cigarstring` | |
| `read.flag` | `read.flag` | |
| `read.is_paired`, `is_unmapped`, `is_reverse` | same | |
| `read.get_tag("NM")` | `read.get_tag("NM")` | Tags load when requested via columnar scan; see below |
| `bam.fetch(contig, start, end)` | `bam.fetch(contig, start, end)` or `bam.fetch(region="chr1:1-100")` | Region string is samtools-style (1-based start in string) |
| `bam.count()` | `bam.count()` | |
| `bam.pileup(...)` | `bam.pileup(..., reads=False)` for depth-only | Requires htslib build; `reads=False` skips per-read materialization |
| `bam.header` (dict-like) | `bam.header()` → `dict[str, int]` | Reference name → length |
| `bam.references` | `bam.references()` | |
| N/A | `bam.to_arrow(...)`, `bamboo.read_columns(...)` | **Preferred for analytics** |

## When to keep iterating vs go columnar

**Keep record iteration** when you:
- Process a small number of alignments
- Need full read objects with complex per-record logic
- Are porting existing pysam code verbatim (use `bamboo.compat.pysam`)

**Switch to columnar** when you:
- Filter or aggregate over millions of reads (QC, coverage summaries, cohort stats)
- Want Polars/Pandas/Arrow without Python per-record overhead
- Only need a subset of columns (`qname`, `pos`, `mapq`, `flag`, tags)

```python
import bamboo as bm

# One-shot: Rust scan → PyArrow (fastest for large regions)
table = bm.read_columns(
    "cohort.bam",
    columns=["qname", "rname", "pos", "mapq", "flag"],
    tags=["NM"],
    region="chr1:1000000-5000000",
    min_mapq=30,
)
df = bm.to_polars(table)
```

This is the workflow Bamboo is optimized for. See `examples/cohort_region_qc.py`.

## CRAM

```python
# pysam
import pysam
with pysam.AlignmentFile("sample.cram", "rc", reference_filename="ref.fasta") as cram:
    ...

# bamboo
import bamboo as bm
with bm.CramFile("sample.cram", reference_filename="ref.fasta") as cram:
    ...

# columnar CRAM (recommended)
table = bm.read_cram_columns(
    "sample.cram",
    columns=["qname", "rname", "pos", "mapq"],
    region="chr1:1000000-5000000",
    reference_filename="ref.fasta",
)
```

## Cloud paths

Bamboo reads `s3://`, `gs://`, `https://`, and `file://` URIs with the same API. Ensure the `.bai` / `.crai` index is co-located (e.g. `sample.bam.bai` next to `sample.bam`).

## Validation checklist

Before switching production pipelines:

1. Run parity on your files:
   ```bash
   pytest tests/test_pysam_parity.py tests/test_pysam_realworld.py -q
   ```
2. Spot-check a region:
   ```python
   import bamboo as bm, pysam

   region = "chr1:1000000-5000000"
   cols = ["qname", "pos", "mapq", "flag", "cigar"]

   bm_table = bm.read_columns("your.bam", columns=cols, region=region)
   # compare against pysam loop — counts and values should match
   ```
3. Confirm index sidecars are present for `fetch()` / columnar region scans.

## Not yet supported (stick with pysam)

- SAM text mode (`"r"`), CRAM via unified `AlignmentFile("rc")`
- `check_index()`, CSI, tabix on BAM
- `count_coverage()`, `pileup_coverage()`
- Record mutation (`set_tag`, building `AlignedSegment` from scratch)
- SAM/CRAM writing (BAM write is supported)
- Full `Header` object API (`header.to_dict()`, `@PG` records, etc.)

## Recommended migration path

1. **Day 1:** `from bamboo.compat import pysam` — zero code changes, run existing tests
2. **Day 2:** Replace hot loops with `read_columns()` / `to_arrow()` on indexed regions
3. **Day 3:** Adopt `examples/cohort_region_qc.py` pattern for QC notebooks
4. **Later:** Direct `import bamboo as bm`, drop compat shim where columnar covers your use case