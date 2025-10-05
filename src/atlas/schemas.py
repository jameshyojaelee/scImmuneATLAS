"""Pandera schemas for AnnData validation."""

from __future__ import annotations

import logging
from typing import List, Optional

import anndata as ad
import pandera as pa
import pandas as pd


# Define expected obs (cell metadata) schema
obs_schema = pa.DataFrameSchema(
    {
        "dataset_id": pa.Column(str, required=True, nullable=False),
        "cancer_type": pa.Column(str, required=True, nullable=False),
        "platform": pa.Column(str, required=True, nullable=False),
    },
    strict=False,  # Allow additional columns beyond required ones
    coerce=True,
)

# Define expected var (gene metadata) schema
# Gene names should be strings (index), allowing any additional columns
var_schema = pa.DataFrameSchema(
    {},  # No required columns in var for now
    strict=False,
    coerce=True,
)


receptor_table_schema = pa.DataFrameSchema(
    {
        "dataset_id": pa.Column(str, required=True, nullable=False),
        "cell_id": pa.Column(str, required=True, nullable=False),
        "chain": pa.Column(str, required=True, nullable=False),
        "locus": pa.Column(str, required=False, nullable=True),
        "clonotype_id": pa.Column(str, required=False, nullable=True),
        "productive": pa.Column(bool, required=False, nullable=True),
        "cdr3_nt": pa.Column(str, required=False, nullable=True),
        "cdr3_aa": pa.Column(str, required=False, nullable=True),
        "v_gene": pa.Column(str, required=False, nullable=True),
        "d_gene": pa.Column(str, required=False, nullable=True),
        "j_gene": pa.Column(str, required=False, nullable=True),
    },
    strict=False,
    coerce=True,
)


def validate_receptor_table(
    df: pd.DataFrame,
    *,
    chain_sets: Optional[List[List[str]]] = None,
    raise_on_error: bool = True,
) -> tuple[bool, Optional[str]]:
    """Validate receptor tables and optionally enforce chain pairing completeness."""

    try:
        receptor_table_schema.validate(df, lazy=True)

        if chain_sets:
            required_sets = [set(group) for group in chain_sets if group]
            grouped = df.groupby("cell_id")["chain"].agg(lambda values: set(v for v in values if v))
            incomplete = [
                cell_id
                for cell_id, chains in grouped.items()
                if required_sets and not any(req.issubset(chains) for req in required_sets)
            ]
            if incomplete:
                raise ValueError(
                    "Receptor pairing incomplete for cells: "
                    + ", ".join(incomplete[:10])
                    + ("..." if len(incomplete) > 10 else "")
                )

        return True, None

    except (pa.errors.SchemaError, ValueError) as exc:
        message = f"Receptor table validation failed: {exc}"
        if raise_on_error:
            raise
        return False, message


def validate_anndata(
    adata: ad.AnnData,
    dataset_id: Optional[str] = None,
    raise_on_error: bool = True,
) -> tuple[bool, Optional[str]]:
    """
    Validate an AnnData object against expected schemas.

    Parameters
    ----------
    adata : ad.AnnData
        The AnnData object to validate
    dataset_id : Optional[str]
        Dataset ID for logging context
    raise_on_error : bool
        If True, raise exception on validation failure. If False, return (False, error_msg)

    Returns
    -------
    tuple[bool, Optional[str]]
        (is_valid, error_message). error_message is None if valid.

    Raises
    ------
    pa.errors.SchemaError
        If validation fails and raise_on_error is True
    """
    context = f" for dataset {dataset_id}" if dataset_id else ""

    try:
        # Validate obs schema
        obs_schema.validate(adata.obs, lazy=True)
        logging.debug(f"obs schema validation passed{context}")

        # Validate var schema
        var_schema.validate(adata.var, lazy=True)
        logging.debug(f"var schema validation passed{context}")

        # Additional integrity checks
        if adata.n_obs == 0:
            raise ValueError("AnnData has 0 observations (cells)")

        if adata.n_vars == 0:
            raise ValueError("AnnData has 0 variables (genes)")

        # Check for unique indices
        if len(adata.obs_names) != len(set(adata.obs_names)):
            raise ValueError("AnnData obs_names are not unique")

        if len(adata.var_names) != len(set(adata.var_names)):
            raise ValueError("AnnData var_names are not unique")

        logging.info(f"AnnData validation passed{context}: {adata.n_obs} cells, {adata.n_vars} genes")
        return True, None

    except (pa.errors.SchemaError, ValueError) as e:
        error_msg = f"Validation failed{context}: {str(e)}"
        logging.error(error_msg)

        if raise_on_error:
            raise

        return False, error_msg


def validate_anndata_list(
    adata_list: list[ad.AnnData],
    dataset_ids: Optional[list[str]] = None,
    raise_on_error: bool = True,
) -> list[tuple[bool, Optional[str]]]:
    """
    Validate a list of AnnData objects.

    Parameters
    ----------
    adata_list : list[ad.AnnData]
        List of AnnData objects to validate
    dataset_ids : Optional[list[str]]
        Optional list of dataset IDs for logging
    raise_on_error : bool
        If True, raise exception on first validation failure

    Returns
    -------
    list[tuple[bool, Optional[str]]]
        List of (is_valid, error_message) tuples for each dataset
    """
    if dataset_ids is None:
        dataset_ids = [None] * len(adata_list)

    results = []
    for adata, dataset_id in zip(adata_list, dataset_ids):
        result = validate_anndata(adata, dataset_id=dataset_id, raise_on_error=raise_on_error)
        results.append(result)

    return results
