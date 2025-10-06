.PHONY: setup lint test demo all app env-bootstrap fetch-no-cap run-after-fetch validate-data pipeline receptor tcr report

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
	scimmuneatlas pipeline --config config/atlas.yaml --jobs 8

pipeline:
	scimmuneatlas pipeline --config config/atlas.yaml

receptor:
	scimmuneatlas receptor --config config/atlas.yaml --stage all

tcr:
	scimmuneatlas tcr --config config/atlas.yaml

report:
	scimmuneatlas report --config config/atlas.yaml

validate-data:
	scimmuneatlas validate-data --config config/atlas.yaml

app:
	streamlit run app/streamlit_app.py

env-bootstrap:
	bash scripts/bootstrap_micromamba.sh

fetch-no-cap:
	bash scripts/fetch_cellxgene_no_cap.sh

run-after-fetch: pipeline
