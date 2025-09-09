#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)"
LOG_DIR="$PROJECT_ROOT/logs"
mkdir -p "$LOG_DIR"

echo "Running integration and annotation pipeline..."

PY="python"
SM="snakemake"
# Detect micromamba (prefer project-local, fallback to PATH)
MAMBA_BIN="$PROJECT_ROOT/.micromamba/bin/micromamba"
if [[ -x "$MAMBA_BIN" ]]; then
  MAMBA="$MAMBA_BIN"
elif command -v micromamba >/dev/null 2>&1; then
  MAMBA="$(command -v micromamba)"
fi

if [[ -n "${MAMBA:-}" ]]; then
  ROOT_PREFIX="$PROJECT_ROOT/.micromamba"
  if [[ -d "$ROOT_PREFIX" ]]; then
    export MAMBA_ROOT_PREFIX="$ROOT_PREFIX"
    PY="$MAMBA -r $ROOT_PREFIX run -n immune-atlas python"
    SM="$MAMBA -r $ROOT_PREFIX run -n immune-atlas snakemake"
  else
    PY="$MAMBA run -n immune-atlas python"
    SM="$MAMBA run -n immune-atlas snakemake"
  fi
  export PYTHONNOUSERSITE=1
fi

# Use snakemake if available; otherwise call modules directly
if command -v $SM >/dev/null 2>&1 || [[ "$SM" == *micromamba* ]]; then
  $SM -j 6 | tee -a "$LOG_DIR/pipeline.log"
else
  echo "snakemake not found in PATH, running modules directly" | tee -a "$LOG_DIR/pipeline.log"
  $PY -m src.atlas.qc --run --config "$PROJECT_ROOT/config/atlas.yaml" | tee -a "$LOG_DIR/qc.log"
  $PY -m src.atlas.doublets --run --config "$PROJECT_ROOT/config/atlas.yaml" | tee -a "$LOG_DIR/doublets.log"
  $PY -m src.atlas.integration --config "$PROJECT_ROOT/config/atlas.yaml" | tee -a "$LOG_DIR/integration.log"
  $PY -m src.atlas.annotate --config "$PROJECT_ROOT/config/atlas.yaml" | tee -a "$LOG_DIR/annotation.log"
  $PY -m src.atlas.export --cellxgene --config "$PROJECT_ROOT/config/atlas.yaml" | tee -a "$LOG_DIR/export.log"
  $PY -m src.atlas.viz --all --config "$PROJECT_ROOT/config/atlas.yaml" | tee -a "$LOG_DIR/viz.log"
fi

echo "Pipeline completed."

