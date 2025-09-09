"""Fetch scRNA-seq subsets from CELLxGENE Census and save as H5AD.

This utility lets us programmatically materialize disease-specific immune subsets
as local `.h5ad` files compatible with the pipeline (adds `dataset_id` and
`cancer_type` to `obs`).
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import List

import anndata as ad
import numpy as np
import pandas as pd


def setup_logging(level: str = "INFO") -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def build_obs_filter(
    disease_terms: List[str],
    immune_only: bool,
    tissue_terms: List[str] | None = None,
    tissue_general_terms: List[str] | None = None,
) -> str:
    """Build a robust obs_value_filter using disease and optional tissue filters.

    Uses only supported ops (==, !=, in, and, or). Avoids functions like lower()/contains().
    """
    clauses: list[str] = []

    if disease_terms:
        diseases_list = ", ".join([f"'{d}'" for d in disease_terms])
        clauses.append(f"disease in [{diseases_list}]")

    if tissue_terms:
        tissues_list = ", ".join([f"'{t}'" for t in tissue_terms])
        clauses.append(f"tissue in [{tissues_list}]")

    if tissue_general_terms:
        tgs_list = ", ".join([f"'{t}'" for t in tissue_general_terms])
        clauses.append(f"tissue_general in [{tgs_list}]")

    # Combine disease/tissue clauses with OR
    core = " or ".join(f"({c})" for c in clauses) if clauses else "True"

    # Optionally require primary data
    if immune_only:
        return f"({core}) and is_primary_data == True"
    return f"({core})"


def fetch_subset(
    dataset_id: str,
    cancer_type: str,
    disease_terms: List[str],
    out_path: Path,
    immune_only: bool = True,
    organism: str = "homo_sapiens",
    max_cells: int | None = 50000,
    tissue_terms: List[str] | None = None,
    tissue_general_terms: List[str] | None = None,
) -> Path:
    import cellxgene_census

    setup_logging()
    ensure_dir(out_path.parent)

    logging.info(
        f"Opening CELLxGENE Census for organism={organism}, dataset_id={dataset_id}, immune_only={immune_only}"
    )

    obs_filter = build_obs_filter(
        disease_terms,
        immune_only,
        tissue_terms=tissue_terms,
        tissue_general_terms=tissue_general_terms,
    )

    with cellxgene_census.open_soma() as census:
        logging.info(f"Querying obs with filter: {obs_filter}")
        # Request minimal, schema-robust fields
        adata: ad.AnnData = cellxgene_census.get_anndata(
            census,
            organism=organism,
            obs_value_filter=obs_filter,
            obs_column_names=[
                "cell_type",
                "tissue",
                "tissue_general",
                "disease",
                "assay",
                "is_primary_data",
            ],
        )

    if adata.n_obs == 0:
        # Fallback: relax primary filter and add common kidney tissue hints for RCC
        logging.warning("No cells matched the query; retrying with relaxed filter")
        alt_filter = build_obs_filter(
            disease_terms,
            immune_only=False,
            tissue_terms=(tissue_terms or []) + ["kidney"],
            tissue_general_terms=(tissue_general_terms or []) + ["kidney"],
        )
        with cellxgene_census.open_soma() as census:
            logging.info(f"Retrying obs with filter: {alt_filter}")
            adata = cellxgene_census.get_anndata(
                census,
                organism=organism,
                obs_value_filter=alt_filter,
                obs_column_names=[
                    "cell_type",
                    "tissue",
                    "tissue_general",
                    "disease",
                    "assay",
                    "is_primary_data",
                ],
            )
        if adata.n_obs == 0:
            raise RuntimeError("No cells matched after relaxed filter; refine terms.")

    # Subsample only if a positive cap was requested
    if max_cells is not None and max_cells > 0 and adata.n_obs > max_cells:
        logging.info(f"Subsampling cells: {adata.n_obs} -> {max_cells}")
        rng = np.random.default_rng(0)
        idx = rng.choice(adata.n_obs, size=max_cells, replace=False)
        adata = adata[idx, :].copy()

    # Ensure names are unique
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    # Add pipeline-required metadata
    adata.obs["dataset_id"] = dataset_id
    adata.obs["cancer_type"] = cancer_type

    logging.info(f"Writing H5AD to {out_path} ({adata.n_obs} cells, {adata.n_vars} genes)")
    adata.write_h5ad(out_path)
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch subsets from CELLxGENE Census")
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--cancer-type", required=True)
    parser.add_argument(
        "--disease-term",
        action="append",
        required=True,
        help="One or more disease search terms (append multiple times)",
    )
    parser.add_argument("--out", required=True, help="Output .h5ad path")
    parser.add_argument(
        "--tissue-term",
        action="append",
        help="Optional tissue terms (append multiple times)",
    )
    parser.add_argument(
        "--tissue-general-term",
        action="append",
        help="Optional tissue_general terms (append multiple times)",
    )
    parser.add_argument("--no-immune-only", action="store_true", help="Do not attempt immune-only filter")
    parser.add_argument(
        "--max-cells",
        type=int,
        default=50000,
        help="Maximum number of cells to keep. Use 0 or negative to keep all.",
    )
    parser.add_argument(
        "--all-cells",
        action="store_true",
        help="Ignore any max cell cap and keep all matching cells",
    )

    args = parser.parse_args()
    immune_only = not args.no_immune_only
    out_path = Path(args.out)
    # Interpret unlimited request
    max_cells: int | None
    if args.all_cells or (args.max_cells is not None and args.max_cells <= 0):
        max_cells = None
    else:
        max_cells = args.max_cells

    fetch_subset(
        dataset_id=args.dataset_id,
        cancer_type=args.cancer_type,
        disease_terms=args.disease_term,
        out_path=out_path,
        immune_only=immune_only,
        max_cells=max_cells,
        tissue_terms=args.tissue_term or None,
        tissue_general_terms=args.tissue_general_term or None,
    )


if __name__ == "__main__":
    main()

