#!/usr/bin/env bash
# Build a release wheel and smoke-test install in a clean venv.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${ROOT}/dist"
VENV="${ROOT}/.wheel-test-venv"

cd "${ROOT}"

rm -rf "${OUT}"
mkdir -p "${OUT}"

echo "== Building wheel =="
python -m pip install -q maturin pyarrow
maturin build --release --strip --features htslib -o "${OUT}"

echo "== Building sdist =="
maturin sdist -o "${OUT}"

echo "== Smoke install =="
rm -rf "${VENV}"
python -m venv "${VENV}"
# shellcheck disable=SC1091
source "${VENV}/bin/activate"
python -m pip install -q pyarrow
python -m pip install "${OUT}"/bamboo-*.whl

python - <<'PY'
import bamboo as bm

assert bm.__version__
assert bm.read_columns is not None
assert bm.pileup_available(), "expected htslib-enabled wheel"
print("OK:", bm.__version__, "pileup=", bm.pileup_available())
PY

echo "== Artifacts =="
ls -lh "${OUT}"