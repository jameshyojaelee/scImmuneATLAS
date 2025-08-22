.PHONY: setup lint test demo all app

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

app:
	streamlit run app/streamlit_app.py
