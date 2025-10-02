.PHONY: setup lint test demo all app env-bootstrap fetch-no-cap run-after-fetch validate-data

ENV_NAME=immune-atlas

setup:
	mamba env create -f env.yml || true
	@echo "Activate with: mamba activate $(ENV_NAME)"
	pre-commit install || true

lint:
	black src tests app --check
	isort src tests app --check-only
	ruff src tests app
	flake8 src tests app

test:
	pytest -q

demo:
	python -c "from src.atlas.utils import run_demo; run_demo('config/atlas.yaml')"

all:
	snakemake -j 8

report:
	python -m atlas.cli report

validate-data:
	scimmuneatlas validate-data --config config/atlas.yaml

app:
	streamlit run app/streamlit_app.py

env-bootstrap:
	bash scripts/bootstrap_micromamba.sh

fetch-no-cap:
	bash scripts/fetch_cellxgene_no_cap.sh

run-after-fetch:
	bash scripts/run_pipeline_after_fetch.sh
