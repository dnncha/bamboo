# Packaging and releases

Release 0.1.0 is [published on PyPI](https://pypi.org/project/bamboo-hts/0.1.0/)
as `bamboo-hts`; users import `bamboo`. The published files include CPython 3.12
wheels for macOS x86_64/ARM64 and manylinux_2_34 x86_64/ARM64, plus an sdist.
There is no Windows wheel in this release. A source distribution's presence
does not establish that it builds on every platform listed in package metadata.

See [the README](README.md#install) for an isolated installation and extension
load check. Matching wheels contain the compiled extension; a source build
requires Rust, maturin, and native build dependencies. The release workflow
uses the `htslib` feature. Check `bamboo.pileup_available()` before relying on it.

The repository contains [Bioconda metadata](conda/bioconda/meta.yaml). Verify a
package in the channel index before documenting `conda install` as available.

## Maintainer release checklist

### 1. Verify a candidate, then bump the version

Sync version in:
- `pyproject.toml` (`project.version`)
- `Cargo.toml` (`workspace.package.version`)
- `conda/recipe/meta.yaml`
- `conda/bioconda/meta.yaml`

### 2. Tag and publish to PyPI

```bash
git tag v0.1.0
git push origin v0.1.0
```

The [Release workflow](.github/workflows/release.yml) builds:
- manylinux wheels (`x86_64`, `aarch64`)
- macOS wheels (`x86_64`, `aarch64`)
- Windows wheel attempt (allowed to fail; not a release prerequisite)
- source distribution (`bamboo-hts-0.1.0.tar.gz`)

and uploads to PyPI when `PYPI_API_TOKEN` is set in GitHub repository secrets.

Create a PyPI API token at https://pypi.org/manage/account/token/ with scope for the `bamboo-hts` project (or entire account for first upload), then:

```bash
gh secret set PYPI_API_TOKEN --repo dnncha/bamboo
```

Re-run the failed workflow for the same reviewed tag after fixing configuration.
Do not move a published tag or infer publication from successful build jobs;
check the PyPI file list and install the artifact on its supported platform.

### 3. Verify locally before tagging

```bash
./scripts/verify_wheel.sh
```

### 4. Submit / update Bioconda

After PyPI publish completes:

```bash
# Get sdist checksum
curl -L -o /tmp/bamboo-hts-0.1.0.tar.gz \
  https://pypi.org/packages/source/b/bamboo-hts/bamboo-hts-0.1.0.tar.gz
shasum -a 256 /tmp/bamboo-hts-0.1.0.tar.gz
```

Update `conda/bioconda/meta.yaml`:
- `version`
- `sha256`
- increment `build.number` if re-building same version

Open a PR to https://github.com/bioconda/bioconda-recipes adding `recipes/bamboo-hts/meta.yaml` (copy from `conda/bioconda/meta.yaml`).

### 5. Local conda build (dev)

```bash
conda install -c conda-forge conda-build maturin rust
conda build conda/recipe --output-folder conda/dist
conda install --use-local -c file://$PWD/conda/dist bamboo-hts
```

## CI requirements

- **CI** (`.github/workflows/ci.yml`): `libhts-dev` for faster native builds in dev; wheels bundle htslib independently.
- **Release**: uses `PyO3/maturin-action` with `--features htslib` (pileup + CRAM htslib paths enabled in published wheels).

## Troubleshooting

| Problem | Fix |
|---------|-----|
| `import bamboo` fails after pip install | Ensure Python >= 3.10; wheel must match platform (no pure-Python fallback) |
| `pileup_available()` is False | Check the installed artifact and build feature; verify the compiled extension loaded |
| Conda build can't find `maturin` | Add `maturin >=1.5,<2` to host requirements |
| htslib compile errors on Linux | Ensure `zlib`, `bzip2`, `xz`, `libdeflate`, `libcurl` host deps (see recipe) |
| Windows wheel build fails in CI | `hts-sys` runs `version.sh` at compile time — Windows wheels are optional until upstream fixes land |