#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)"
LOG_DIR="$PROJECT_ROOT/logs"
LOG="$LOG_DIR/orchestrator.log"
ROOT_PREFIX="$PROJECT_ROOT/.micromamba"
MAMBA_BIN="$ROOT_PREFIX/bin/micromamba"
CLI="scimmuneatlas"

if [[ -x "$MAMBA_BIN" ]]; then
  CLI="$MAMBA_BIN -r $ROOT_PREFIX run -n immune-atlas scimmuneatlas"
elif command -v micromamba >/dev/null 2>&1; then
  CLI="$(command -v micromamba) run -n immune-atlas scimmuneatlas"
fi

mkdir -p "$LOG_DIR"

echo "[$(date +%F' '%T)] Watcher start" >>"$LOG"

# Wait for all fetch jobs to finish
while pgrep -fa "micromamba .* run -n immune-atlas python -m src.atlas.fetch_census" >/dev/null; do
  sleep 60
done

echo "[$(date +%F' '%T)] Fetch complete. Upgrading numpy/numba..." >>"$LOG"

# Attempt to upgrade numpy and numba before pipeline (to satisfy numba requirements)
if command -v micromamba >/dev/null 2>&1; then
  micromamba -r "$ROOT_PREFIX" install -y -n immune-atlas "numpy>=1.24" numba >>"$LOG" 2>&1 || true
fi

echo "[$(date +%F' '%T)] Starting pipeline..." >>"$LOG"
$CLI pipeline --config "$PROJECT_ROOT/config/atlas.yaml" --jobs 8 >>"$LOG" 2>&1 || true
echo "[$(date +%F' '%T)] Watcher done" >>"$LOG"

