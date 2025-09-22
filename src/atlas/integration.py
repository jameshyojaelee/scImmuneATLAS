"""Integration functions for the Single-cell Immune Atlas."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional

import anndata as ad
import numpy as np
import scanpy as sc
from scipy import sparse as sp

from .utils import ensure_dir, load_config, set_seed, setup_logging, timer


def _resolve_use_gpu(flag: object) -> bool:
    """Determine whether GPU should be used based on config flag."""
    if isinstance(flag, str):
        flag_lower = flag.lower()
        if flag_lower == "auto":
            try:
                import torch

                return torch.cuda.is_available()
            except ImportError:
                return False
        return flag_lower in {"true", "1", "yes"}
    return bool(flag)


def _compute_shared_hvgs(
    adata: ad.AnnData,
    batch_key: str,
    n_hvg: int,
) -> np.ndarray:
    """Compute highly variable genes shared across batches."""
    logging.info("Computing %d shared HVGs using Seurat v3 flavor", n_hvg)
    tmp = adata.copy()
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)
    sc.pp.highly_variable_genes(
        tmp,
        n_top_genes=n_hvg,
        flavor="seurat_v3",
        batch_key=batch_key,
    )
    hvg_mask = tmp.var["highly_variable"].values
    logging.info("Identified %d HVGs", int(hvg_mask.sum()))
    return hvg_mask


def _apply_neighbors_umap(
    adata: ad.AnnData,
    config: Dict,
    representation_key: Optional[str],
) -> None:
    """Compute neighbors and UMAP on the given representation."""
    neighbors_cfg = config.get("neighbors", {})
    umap_cfg = config.get("umap", {})

    n_neighbors = neighbors_cfg.get("n_neighbors", 15)
    n_pcs = neighbors_cfg.get("n_pcs", 50)
    random_state = config.get("seed", 0)

    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=None if representation_key else n_pcs,
        use_rep=representation_key,
        random_state=random_state,
    )
    sc.tl.umap(
        adata,
        min_dist=umap_cfg.get("min_dist", 0.4),
        spread=umap_cfg.get("spread", 1.0),
        random_state=random_state,
    )


def integrate_scvi(
    adata: ad.AnnData,
    config: Dict,
    hvg_mask: Optional[np.ndarray] = None,
) -> ad.AnnData:
    """Integrate datasets using scVI."""
    try:
        import scvi
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError("scvi-tools is required for scVI integration") from exc

    integration_cfg = config["integration"]
    batch_key = integration_cfg["batch_key"]
    latent_dim = integration_cfg.get("latent_dim", 30)
    max_epochs = integration_cfg.get("max_epochs", 200)
    use_gpu = _resolve_use_gpu(integration_cfg.get("use_gpu", "auto"))

    logging.info("Starting scVI integration (use_gpu=%s)", use_gpu)

    adata_scvi = adata.copy()
    adata_scvi.layers["counts"] = adata_scvi.X.copy()

    if hvg_mask is not None:
        adata_scvi.var["highly_variable"] = hvg_mask

    scvi.settings.seed = config.get("seed", 0)
    scvi.settings.dl_gpu_training = use_gpu

    scvi.model.SCVI.setup_anndata(
        adata_scvi,
        layer="counts",
        batch_key=batch_key,
    )

    model = scvi.model.SCVI(adata_scvi, n_latent=latent_dim)
    model.train(max_epochs=max_epochs, early_stopping=True)

    adata_scvi.obsm["X_scvi"] = model.get_latent_representation()
    _apply_neighbors_umap(adata_scvi, config, representation_key="X_scvi")

    logging.info(
        "scVI integration completed: %d cells, %d genes",
        adata_scvi.n_obs,
        adata_scvi.n_vars,
    )
    return adata_scvi


def integrate_harmony(
    adata: ad.AnnData,
    config: Dict,
    hvg_mask: Optional[np.ndarray] = None,
) -> ad.AnnData:
    """Integrate datasets using Harmony."""
    try:
        import harmonypy as hm
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError("harmonypy is required for Harmony integration") from exc

    integration_cfg = config["integration"]
    batch_key = integration_cfg["batch_key"]

    logging.info("Starting Harmony integration")

    # Work on a copy to avoid mutating the raw counts matrix
    adata_norm = adata.copy()

    def _normalize_total_no_numba(a: ad.AnnData, target_sum: float = 1e4) -> None:
        X = a.X
        if sp.issparse(X):
            counts = np.asarray(X.sum(axis=1)).ravel()
            counts[counts == 0] = 1.0
            scaling_factors = target_sum / counts
            a.X = sp.diags(scaling_factors) @ X
        else:
            counts = X.sum(axis=1)
            counts[counts == 0] = 1.0
            scaling_factors = (target_sum / counts).astype(X.dtype)
            a.X = (X.T * scaling_factors).T

    _normalize_total_no_numba(adata_norm)
    sc.pp.log1p(adata_norm)

    if hvg_mask is not None:
        adata_norm = adata_norm[:, hvg_mask].copy()
    else:
        sc.pp.highly_variable_genes(adata_norm, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata_norm = adata_norm[:, adata_norm.var.highly_variable]

    sc.pp.scale(adata_norm, max_value=10)
    sc.tl.pca(adata_norm, svd_solver="arpack", random_state=config.get("seed", 0))

    harmony_out = hm.run_harmony(
        adata_norm.obsm["X_pca"],
        adata.obs,
        batch_key,
        random_state=config.get("seed", 0),
    )

    adata_out = adata.copy()
    adata_out.obsm["X_harmony"] = harmony_out.Z_corr.T
    _apply_neighbors_umap(adata_out, config, representation_key="X_harmony")

    logging.info(
        "Harmony integration completed: %d cells, %d genes",
        adata_out.n_obs,
        adata_out.n_vars,
    )
    return adata_out


def _sample_embedding(embedding: np.ndarray, max_cells: int = 5000) -> np.ndarray:
    if embedding.shape[0] <= max_cells:
        return embedding
    rng = np.random.default_rng(0)
    idx = rng.choice(embedding.shape[0], size=max_cells, replace=False)
    return embedding[idx]


def compute_integration_metrics(
    adata: ad.AnnData,
    config: Dict,
    embedding_key: str,
) -> Dict[str, float]:
    """Compute user-specified integration diagnostics."""
    metrics_requested = config["integration"].get("metrics", [])
    if not metrics_requested:
        return {}

    batch_key = config["integration"]["batch_key"]
    metrics: Dict[str, float] = {}
    embedding = adata.obsm.get(embedding_key)
    if embedding is None:
        logging.warning("Embedding %s not present; skipping integration metrics", embedding_key)
        return metrics

    embedding_sampled = _sample_embedding(embedding)
    batch_labels = adata.obs[batch_key].to_numpy()

    if "silhouette" in metrics_requested:
        try:
            from sklearn.metrics import silhouette_score

            batch_subset = batch_labels[: embedding_sampled.shape[0]]
            if np.unique(batch_subset).size > 1:
                metrics["silhouette_batch"] = float(
                    silhouette_score(embedding_sampled, batch_subset)
                )

            if "cell_type" in adata.obs:
                cell_labels = adata.obs["cell_type"].astype(str).to_numpy()
                cell_subset = cell_labels[: embedding_sampled.shape[0]]
                if np.unique(cell_subset).size > 1:
                    metrics["silhouette_cell_type"] = float(
                        silhouette_score(embedding_sampled, cell_subset)
                    )
        except Exception as exc:  # pragma: no cover - best effort metric
            logging.warning("Failed to compute silhouette metrics: %s", exc)

    if "lisi" in metrics_requested:
        try:  # pragma: no cover - optional dependency
            from scib.metrics import lisi_graph

            lisi_res = lisi_graph(
                adata,
                batch_key=batch_key,
                label_key="cell_type" if "cell_type" in adata.obs else batch_key,
                use_rep=embedding_key,
                n_neighbors=config.get("neighbors", {}).get("n_neighbors", 15),
            )
            metrics["lisi_batch_mean"] = float(np.mean(lisi_res["lisi_batch"]))
            metrics["lisi_label_mean"] = float(np.mean(lisi_res["lisi_label"]))
        except ImportError:
            logging.warning("scIB not available; skipping LISI metric")
        except Exception as exc:
            logging.warning("Failed to compute LISI metric: %s", exc)

    return metrics


def _persist_integration_metrics(metrics: Dict[str, float], config: Dict) -> None:
    metrics_dir = Path(config["outputs"].get("metrics_dir", "processed/metrics"))
    ensure_dir(metrics_dir)
    metrics_path = metrics_dir / "integration_metrics.json"
    with open(metrics_path, "w") as f:
        json.dump(metrics, f, indent=2)
    logging.info("Saved integration metrics to %s", metrics_path)


def run_integration(config: Dict) -> None:
    """Run integration based on configuration settings."""
    adatas: List[ad.AnnData] = []
    for dataset_info in config["datasets"]:
        dataset_id = dataset_info["id"]
        input_path = Path("data/interim") / f"{dataset_id}.doublet_filtered.h5ad"
        logging.info("Loading %s", input_path)
        adatas.append(ad.read_h5ad(input_path))

    if not adatas:
        raise ValueError("No datasets found for integration")

    set_seed(config["seed"])

    adata_concat = ad.concat(adatas, join="outer", index_unique="-")
    adata_concat.var_names_make_unique()

    integration_cfg = config["integration"]
    hvg_mask: Optional[np.ndarray] = None
    if integration_cfg.get("shared_hvg", True):
        hvg_mask = _compute_shared_hvgs(
            adata_concat,
            integration_cfg["batch_key"],
            integration_cfg.get("n_hvg", 3000),
        )
        adata_concat.var["highly_variable"] = hvg_mask

    method = integration_cfg["method"].lower()
    if method == "scvi":
        adata_integrated = integrate_scvi(adata_concat, config, hvg_mask=hvg_mask)
        embedding_key = "X_scvi"
    elif method == "harmony":
        adata_integrated = integrate_harmony(adata_concat, config, hvg_mask=hvg_mask)
        embedding_key = "X_harmony"
    else:
        raise ValueError(f"Unknown integration method: {integration_cfg['method']}")

    metrics = compute_integration_metrics(adata_integrated, config, embedding_key)
    _persist_integration_metrics(metrics, config)

    output_path = Path("processed") / "integrated_atlas.h5ad"
    ensure_dir(output_path.parent)
    adata_integrated.write(output_path)
    logging.info("Saved integrated atlas to %s", output_path)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run dataset integration")
    parser.add_argument("--method", choices=["scvi", "harmony"], help="Integration method")
    parser.add_argument("--config", default="config/atlas.yaml", help="Config file path")

    args = parser.parse_args()

    setup_logging()
    cfg = load_config(args.config)

    if args.method:
        cfg["integration"]["method"] = args.method

    with timer(f"Integration using {cfg['integration']['method']}"):
        run_integration(cfg)
