"""I/O functions for the Single-cell Immune Atlas."""

import hashlib
import logging
import tempfile
from pathlib import Path
from typing import Dict, List, Optional
from urllib.parse import urlparse

import anndata as ad
import numpy as np
import pandas as pd
import requests
import scanpy as sc
from scipy import io as scipy_io

from .schemas import validate_anndata
from .utils import ensure_dir, load_config, setup_logging, timer


def _compute_sha256(filepath: Path) -> str:
    """Compute SHA256 hash of a file."""
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        # Read file in chunks to handle large files
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()


def download_if_needed(entry: Dict, outdir: Path) -> Dict:
    """
    Download dataset files if they are URLs, otherwise use local paths.

    Supports optional SHA256 verification via 'sha256' field in entry dict.
    When a hash is present, downloads to temp file, verifies checksum,
    then atomically moves to final destination.

    Parameters
    ----------
    entry : Dict
        Dataset entry with url/genes/barcodes/metadata paths and optional sha256 hashes
    outdir : Path
        Output directory for downloaded files

    Returns
    -------
    Dict
        Entry dict with URLs replaced by local paths
    """
    result = entry.copy()

    for key in ["url", "genes", "barcodes", "metadata"]:
        if key not in entry:
            continue

        path_or_url = entry[key]
        parsed = urlparse(path_or_url)

        if parsed.scheme in ["http", "https"]:
            # It's a URL - download it
            filename = Path(parsed.path).name
            if not filename:
                filename = f"{entry['id']}_{key}.txt"

            outpath = outdir / filename

            # Check for optional SHA256 hash for this key
            hash_key = f"{key}_sha256"
            expected_sha256 = entry.get(hash_key)

            if not outpath.exists():
                logging.info(f"Downloading {path_or_url} to {outpath}")
                response = requests.get(path_or_url)
                response.raise_for_status()

                if expected_sha256:
                    # Download to temp file, verify hash, then move atomically
                    with tempfile.NamedTemporaryFile(
                        mode="wb", delete=False, dir=outdir
                    ) as tmp_file:
                        tmp_path = Path(tmp_file.name)
                        tmp_file.write(response.content)

                    # Verify SHA256
                    actual_sha256 = _compute_sha256(tmp_path)
                    if actual_sha256 != expected_sha256:
                        tmp_path.unlink()  # Clean up temp file
                        raise ValueError(
                            f"SHA256 mismatch for {key} in dataset {entry['id']}: "
                            f"expected {expected_sha256}, got {actual_sha256}"
                        )

                    # Atomic move to final destination
                    tmp_path.replace(outpath)
                    logging.info(f"SHA256 verification passed for {outpath}")
                else:
                    # No hash verification - write directly
                    with open(outpath, "wb") as f:
                        f.write(response.content)
            else:
                # File already exists - optionally verify hash if provided
                if expected_sha256:
                    actual_sha256 = _compute_sha256(outpath)
                    if actual_sha256 != expected_sha256:
                        logging.warning(
                            f"Existing file {outpath} has incorrect SHA256. "
                            f"Expected {expected_sha256}, got {actual_sha256}. "
                            f"Consider deleting and re-downloading."
                        )

            result[key] = str(outpath)
        else:
            # It's a local path - use as is
            result[key] = path_or_url

    return result


def load_matrix(entry: Dict) -> ad.AnnData:
    """Load count matrix from various formats and ensure proper metadata."""
    dataset_id = entry["id"]
    cancer_type = entry["cancer_type"]
    platform = entry.get("platform", "unknown")
    
    if "url" in entry and entry["url"].endswith(".h5ad"):
        # Load H5AD directly
        adata = ad.read_h5ad(entry["url"])
    elif "url" in entry and entry["url"].endswith(".mtx.gz"):
        # Load MTX format
        X = scipy_io.mmread(entry["url"]).T.tocsr()  # Transpose to cells x genes
        adata = ad.AnnData(X)
        
        # Load gene names if provided
        if "genes" in entry:
            genes_df = pd.read_csv(entry["genes"], sep="\t", header=None)
            if genes_df.shape[1] >= 2:
                adata.var_names = genes_df.iloc[:, 1].values  # Use gene symbols
            else:
                adata.var_names = genes_df.iloc[:, 0].values
        
        # Load cell barcodes if provided
        if "barcodes" in entry:
            barcodes_df = pd.read_csv(entry["barcodes"], sep="\t", header=None)
            adata.obs_names = barcodes_df.iloc[:, 0].values
    else:
        raise ValueError(f"Unsupported format for dataset {dataset_id}")
    
    # Ensure required metadata
    adata.obs["dataset_id"] = dataset_id
    adata.obs["cancer_type"] = cancer_type
    adata.obs["platform"] = platform
    
    # Load additional metadata if provided
    if "metadata" in entry:
        meta_df = pd.read_csv(entry["metadata"], sep="\t", index_col=0)
        # Join with existing obs
        adata.obs = adata.obs.join(meta_df, how="left")
    
    # Make names unique
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    # Validate AnnData against schema
    validate_anndata(adata, dataset_id=dataset_id, raise_on_error=True)

    return adata


def download_datasets(config: Dict) -> List[str]:
    """Download all datasets specified in config and return list of dataset IDs."""
    raw_dir = Path("data/raw")
    ensure_dir(raw_dir)
    
    dataset_ids = []
    
    for dataset_info in config["datasets"]:
        dataset_id = dataset_info["id"]
        logging.info(f"Processing dataset: {dataset_id}")
        
        # Download files if needed
        downloaded_info = download_if_needed(dataset_info, raw_dir)
        
        # Create a flag file to indicate download completion
        flag_file = raw_dir / f"{dataset_id}_downloaded.flag"
        with open(flag_file, "w") as f:
            f.write(f"Downloaded at {pd.Timestamp.now()}\n")
            f.write(f"Config: {downloaded_info}\n")
        
        dataset_ids.append(dataset_id)
    
    return dataset_ids


def write_metadata_table(config: Dict, outpath: str = "data/external/metadata.tsv") -> None:
    """Write a metadata table summarizing all datasets."""
    ensure_dir(Path(outpath).parent)
    
    rows = []
    for dataset_info in config["datasets"]:
        rows.append({
            "dataset_id": dataset_info["id"],
            "modality": dataset_info.get("modality", "scRNA"),
            "platform": dataset_info.get("platform", "unknown"),
            "cancer_type": dataset_info["cancer_type"],
            "url": dataset_info.get("url", ""),
        })
    
    df = pd.DataFrame(rows)
    df.to_csv(outpath, sep="\t", index=False)
    logging.info(f"Wrote metadata table to {outpath}")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Download and prepare datasets")
    parser.add_argument("--download", action="store_true", help="Download datasets")
    parser.add_argument("--config", default="config/atlas.yaml", help="Config file path")
    
    args = parser.parse_args()
    
    setup_logging()
    config = load_config(args.config)
    
    if args.download:
        with timer("Downloading datasets"):
            dataset_ids = download_datasets(config)
            write_metadata_table(config)
            logging.info(f"Downloaded {len(dataset_ids)} datasets: {dataset_ids}")
