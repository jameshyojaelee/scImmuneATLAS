#!/usr/bin/env bash
set -euo pipefail

# Bootstrap micromamba locally in project dir and create/activate env

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)"
MAMBA_DIR="$PROJECT_ROOT/.micromamba"
ENV_NAME="immune-atlas"
ENV_FILE_CPU="$PROJECT_ROOT/env.cpu.yml"

mkdir -p "$MAMBA_DIR"

if [[ ! -x "$MAMBA_DIR/bin/micromamba" ]]; then
  echo "Installing micromamba to $MAMBA_DIR..."
  curl -Ls https://micro.mamba.pm/install.sh | bash -s - -b -p "$MAMBA_DIR"
else
  echo "micromamba already installed at $MAMBA_DIR"
fi

export PATH="$MAMBA_DIR/bin:$PATH"
export MAMBA_ROOT_PREFIX="$MAMBA_DIR"

echo "Creating/updating env $ENV_NAME from $ENV_FILE_CPU..."
micromamba create -y -n "$ENV_NAME" -f "$ENV_FILE_CPU" || true

eval "$(micromamba shell hook -s bash)"
micromamba activate "$ENV_NAME"
export PYTHONNOUSERSITE=1

python -V
python -c "import anndata; import numpy as np; print('anndata', anndata.__version__, 'numpy', np.__version__)"

echo "Environment $ENV_NAME is ready."

