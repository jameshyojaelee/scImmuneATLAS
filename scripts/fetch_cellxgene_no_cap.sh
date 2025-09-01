#!/usr/bin/env bash
set -euo pipefail

# Requires micromamba env prepared via scripts/bootstrap_micromamba.sh

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)"
# Prefer executing via micromamba to ensure correct Python/site-packages
PY="python"
if [[ -x "$PROJECT_ROOT/.micromamba/bin/micromamba" ]]; then
  PY="$PROJECT_ROOT/.micromamba/bin/micromamba run -n immune-atlas python"
  export PYTHONNOUSERSITE=1
fi
LOG_DIR="$PROJECT_ROOT/logs/fetch"
OUT_DIR="$PROJECT_ROOT/data/raw"
mkdir -p "$LOG_DIR" "$OUT_DIR"

# Datasets to fetch (dataset_id cancer_type output_basename disease_terms...)
declare -a JOBS=(
  "MELANOMA_CENSUS melanoma melanoma_census melanoma,melanoma progression"
  "NSCLC_CENSUS NSCLC nsclc_census non-small cell lung carcinoma,lung adenocarcinoma,lung squamous cell carcinoma"
  "BREAST_CENSUS Breast breast_census breast carcinoma,breast cancer"
  "RCC_NIVOLUMAB_CENSUS RCC rcc_nivo_census renal cell carcinoma,clear cell renal cell carcinoma"
  "UROTHELIAL_ATEZOLIZUMAB_CENSUS Urothelial urothelial_atezo_census urothelial carcinoma,bladder cancer"
  "PDAC_ICI_CENSUS PDAC pdac_ici_census pancreatic ductal adenocarcinoma,pancreatic cancer"
)

echo "Launching ${#JOBS[@]} CELLxGENE fetch jobs (no cell cap)"

PIDS=()
N=0
for entry in "${JOBS[@]}"; do
  dataset_id=$(echo "$entry" | awk '{print $1}')
  cancer_type=$(echo "$entry" | awk '{print $2}')
  out_base=$(echo "$entry" | awk '{print $3}')
  terms_csv=$(echo "$entry" | cut -d' ' -f4-)

  out="${OUT_DIR}/${out_base}.h5ad"
  log="${LOG_DIR}/${dataset_id}.log"

  cmd=($PY -m src.atlas.fetch_census \
      --dataset-id "$dataset_id" \
      --cancer-type "$cancer_type" \
      --out "$out" \
      --all-cells)

  # Expand comma-separated disease terms into multiple --disease-term args
  IFS=',' read -ra TERMS <<< "$terms_csv"
  for t in "${TERMS[@]}"; do
    term_trimmed=$(echo "$t" | sed 's/^ *//;s/ *$//')
    cmd+=(--disease-term "$term_trimmed")
  done

  echo "[$(date +%F' '%T)] Starting $dataset_id -> $out" | tee -a "$log"
  ("${cmd[@]}" >>"$log" 2>&1 && echo "[$(date +%F' '%T)] DONE $dataset_id" >>"$log") &
  PIDS+=("$!")
  N=$((N+1))
done

echo "Launched $N jobs. Waiting..."

FAIL=0
for pid in "${PIDS[@]}"; do
  if ! wait "$pid"; then
    FAIL=1
  fi
done

if [[ "$FAIL" -ne 0 ]]; then
  echo "One or more fetch jobs failed. Check logs in $LOG_DIR" >&2
  exit 1
fi

echo "All fetch jobs completed successfully."

