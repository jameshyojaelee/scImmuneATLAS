#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)"
LOG_DIR="$PROJECT_ROOT/logs"
LOG="$LOG_DIR/orchestrator.log"
ROOT_PREFIX="$PROJECT_ROOT/.micromamba"

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
bash "$PROJECT_ROOT/scripts/run_pipeline_after_fetch.sh" >>"$LOG" 2>&1 || true
echo "[$(date +%F' '%T)] Watcher done" >>"$LOG"


