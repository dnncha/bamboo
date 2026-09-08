# Bamboo

Bamboo provides Python access to sequencing alignments and variants through a
Rust extension. It supports record iteration and Arrow tables for workflows that
select, filter, or aggregate alignment columns in Python, pandas, or Polars.

This is an alpha library. Its compatibility module exposes selected pysam-style
entry points; it does not implement the full pysam API. Start with one read path
and compare the fields your analysis consumes.

[Migration guide](MIGRATION.md) · [Packaging and platforms](PACKAGING.md) ·
[Examples](examples/) · [Cheerful Duck Research](https://cheerfulduck.com/research)

## Install

The distribution is [`bamboo-hts`](https://pypi.org/project/bamboo-hts/0.1.0/);
the Python import is `bamboo`. Release 0.1.0 provides Python 3.12 wheels for
macOS Intel/Apple Silicon and Linux x86_64/ARM64, plus a source distribution.
Linux wheels require glibc 2.34 or newer. Windows wheels are not published for
this release. Other Python versions may require a source build despite the
package's Python >=3.10 metadata.

In a Python 3.12 virtual environment:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install 'bamboo-hts==0.1.0'
```

Check the compiled extension as well as the import: the package can import even
if its extension could not load.

```python
import bamboo as bm

assert bm.AlignmentFile is not None, "Bamboo's compiled extension did not load"
print(bm.__version__)
```

For Polars or pandas adapters, install `'bamboo-hts[polars]==0.1.0'` or
`'bamboo-hts[pandas]==0.1.0'`. See [packaging](PACKAGING.md) for source builds
and the Bioconda recipe; a checked-in recipe is not evidence of a channel release.

## Read an indexed BAM region into Arrow

Provide an indexed BAM and a reference name present in its header. Region
strings use a 1-based inclusive start; numeric `fetch(contig, start, end)` uses
0-based, half-open coordinates.

```python
import bamboo as bm

table = bm.read_columns(
    "cohort.bam",
    columns=["qname", "rname", "pos", "mapq", "flag"],
    region="chr1:1000000-5000000",
    min_mapq=30,
)
print(table.num_rows)
print(table.schema)
```

The returned `pos` column uses 0-based alignment starts. Convert the table with
`bm.to_polars(table)` or `bm.to_pandas(table)` when the corresponding dependency
is installed. [The region-QC example](examples/cohort_region_qc.py) demonstrates
filtering and grouping; Arrow export alone does not establish a speedup or
zero-copy behavior for an entire analysis.

## Current scope

| Task | Entry point | Boundary |
| --- | --- | --- |
| Read BAM records | `AlignmentFile` | Iteration and region fetch; indexed fetch needs a reachable BAI |
| Export BAM columns | `read_columns`, `AlignmentFile.to_arrow` | Select columns and tags explicitly |
| Read CRAM | `CramFile`, `read_cram_columns` | Dedicated CRAM API; supply the required reference |
| Read variants | `VariantFile`, `read_vcf_table`, `read_bcf_table` | VCF/BCF readers and table export are implemented |
| Write BAM | `AlignmentFile(..., "wb")` | Local output paths; see the example below |
| Inspect pileup | `pileup` | Requires an htslib-enabled extension |
| Port selected pysam calls | `bamboo.compat.pysam` | Aliases and helpers, not complete API compatibility |

SAM/CRAM writing, the full pysam header and mutation APIs, and unified CRAM
opening through `AlignmentFile(..., "rc")` are outside the documented scope.
See [migration limits](MIGRATION.md#unsupported-interfaces) before porting code.

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


## Evaluate and contribute

The parity tests compare selected fields on small fixtures and a generated
50,000-read dataset. The file named `test_pysam_realworld.py` uses synthetic
cohort data; it does not establish parity on arbitrary production files.
Keep comparisons for your own references, indexes, tags, and filtering rules.

From a source checkout, with Rust and the native build prerequisites installed:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install 'maturin>=1.5,<2' pyarrow pytest 'pysam>=0.22'
maturin develop --release --features htslib
python -m pytest
```

Use [the migration guide](MIGRATION.md) to structure a comparison and
[the wheel verifier](scripts/verify_wheel.sh) to check packaging. Bug reports
should include versions, the command or Python call, expected and observed
fields, and a small shareable fixture.

[Cheerful Duck Research](https://cheerfulduck.com/research) collects our
bioinformatics investigations and reproducible corrections. Bamboo's test suite
and your own workflow comparisons remain the evidence for a migration decision.

## License

MIT. See [LICENSE](LICENSE).
