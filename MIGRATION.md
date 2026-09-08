# Migrating selected pysam read paths to Bamboo

Bamboo's compatibility module re-exports a subset of alignment classes and
helpers. An import change can be a useful first experiment, but differences in
headers, methods, supported modes, and record mutation often require code changes.

## Install and check scope

Follow [installation and platform support](README.md#install), then list the
pysam methods and fields your workflow actually uses. Compare those against the
mapping below before editing an analysis.

## Compare a record-iteration path

```python
# before
import pysam

with pysam.AlignmentFile("sample.bam", "rb") as bam:
    for read in bam.fetch("chr1", 1_000_000, 5_000_000):
        print(read.query_name, read.reference_start, read.flag)

# candidate — compare outputs before adopting
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
| `read.get_tag("NM")` | `read.get_tag("NM")` | Record tag access; select exported tags separately in columnar scans |
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
- Need a supported record interface while porting existing code

**Switch to columnar** when you:
- Filter or aggregate over millions of reads (QC, coverage summaries, cohort stats)
- Want Polars/Pandas/Arrow without Python per-record overhead
- Only need a subset of columns (`qname`, `pos`, `mapq`, `flag`, tags)

```python
import bamboo as bm

# Export selected columns and tags to PyArrow
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

1. Run the repository parity fixtures (these do not read your own files):
   ```bash
   pytest tests/test_pysam_parity.py tests/test_pysam_realworld.py -q
   ```
2. Build a separate comparison on your own representative files. For example,
   obtain Bamboo columns for a region and compare every selected value against
   a pysam loop with the same coordinates and filters:
   ```python
   import bamboo as bm, pysam

   region = "chr1:1000000-5000000"
   cols = ["qname", "pos", "mapq", "flag", "cigar"]

   bm_table = bm.read_columns("your.bam", columns=cols, region=region)
   # compare against pysam loop — counts and values should match
   ```
3. Confirm index sidecars are present for `fetch()` / columnar region scans.

## Unsupported interfaces

- SAM text mode (`"r"`), CRAM via unified `AlignmentFile("rc")`
- `check_index()`, CSI, tabix on BAM
- `count_coverage()`, `pileup_coverage()`
- Record mutation (`set_tag`, building `AlignedSegment` from scratch)
- SAM/CRAM writing (BAM write is supported)
- Full `Header` object API (`header.to_dict()`, `@PG` records, etc.)

## Adopt one verified path at a time

1. Run the original workflow and retain its outputs and package versions.
2. Port one supported read path into a separate comparison script.
3. Compare records, coordinates, tags, ordering, filtering, and failure behavior
   used downstream. Include empty regions and reads at region boundaries.
4. Measure runtime and memory with the same inputs only after the required
   outputs agree. Keep the original route for unsupported interfaces.

For related investigations and reproducible examples, see
[Cheerful Duck Research](https://cheerfulduck.com/research).
