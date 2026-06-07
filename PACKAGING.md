# Packaging and releases

Bamboo ships as **PyPI wheels** (primary) and a **Bioconda recipe** (for conda/mamba environments).

## User install

```bash
# PyPI (recommended)
pip install bamboo

# Optional extras
pip install bamboo[polars]
pip install bamboo[pandas]

# Conda (after bioconda recipe is merged)
conda install -c bioconda -c conda-forge bamboo
```

No maturin or Rust toolchain required for end users — wheels bundle the compiled extension and bundled htslib (via `rust-htslib` / `hts-sys`).

## Maintainer release checklist

### 1. Bump version

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
- Windows wheels
- source distribution (`bamboo-0.1.0.tar.gz`)

and uploads to PyPI when `PYPI_API_TOKEN` is set in GitHub repository secrets.

Create a PyPI API token at https://pypi.org/manage/account/token/ with scope for the `bamboo` project (or entire account for first upload).

### 3. Verify locally before tagging (optional)

```bash
./scripts/verify_wheel.sh
```

### 4. Submit / update Bioconda

After PyPI publish completes:

```bash
# Get sdist checksum
curl -L -o /tmp/bamboo-0.1.0.tar.gz \
  https://pypi.org/packages/source/b/bamboo/bamboo-0.1.0.tar.gz
shasum -a 256 /tmp/bamboo-0.1.0.tar.gz
```

Update `conda/bioconda/meta.yaml`:
- `version`
- `sha256`
- increment `build.number` if re-building same version

Open a PR to https://github.com/bioconda/bioconda-recipes adding `recipes/bamboo/meta.yaml` (copy from `conda/bioconda/meta.yaml`).

### 5. Local conda build (dev)

```bash
conda install -c conda-forge conda-build maturin rust
conda build conda/recipe --output-folder conda/dist
conda install --use-local -c file://$PWD/conda/dist bamboo
```

## CI requirements

- **CI** (`.github/workflows/ci.yml`): `libhts-dev` for faster native builds in dev; wheels bundle htslib independently.
- **Release**: uses `PyO3/maturin-action` with `--features htslib` (pileup + CRAM htslib paths enabled in published wheels).

## Troubleshooting

| Problem | Fix |
|---------|-----|
| `import bamboo` fails after pip install | Ensure Python >= 3.10; wheel must match platform (no pure-Python fallback) |
| `pileup_available()` is False | Reinstall from PyPI wheel built with htslib feature (not a noodles-only dev build) |
| Conda build can't find `maturin` | Add `maturin >=1.5,<2` to host requirements |
| htslib compile errors on Linux | Ensure `zlib`, `bzip2`, `xz`, `libdeflate`, `libcurl` host deps (see recipe) |