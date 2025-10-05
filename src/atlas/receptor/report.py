"""Generate report snippets for receptor analytics."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, List

from .analytics import AnalyticsResult
from .config import global_receptor_config

REPORT_FRAGMENT = "report_section.md"


def _rel_path(path: Path, base: Path) -> str:
    return os.path.relpath(path, base)


def write_report_fragments(config: Dict, result: AnalyticsResult) -> None:
    summary = result.summary_payload
    global_summary = summary.get("global", {})
    if not global_summary:
        # Nothing to report
        fragment_path = result.metrics_dir / REPORT_FRAGMENT
        fragment_path.write_text("\n")
        return

    receptor_cfg = global_receptor_config(config)
    report_cfg = receptor_cfg.get("report", {})
    include_figures: List[str] = report_cfg.get(
        "include_figures",
        [
            "clonotype_frequency",
            "umap_clonal_expansion",
            "cdr3_spectratype",
        ],
    )

    lines: List[str] = []
    lines.append("## Immune Receptor Repertoire")
    lines.append("")
    lines.append(
        f"- Cells with receptor calls: {global_summary.get('n_cells', 0):,}"
    )
    lines.append(
        f"- Unique clonotypes detected: {global_summary.get('n_clonotypes', 0):,}"
    )
    diversity = global_summary.get("diversity", {})
    if diversity:
        lines.append(
            "- Diversity indices (global): "
            f"Shannon={diversity.get('shannon', 0.0):.2f}, "
            f"Simpson={diversity.get('simpson', 0.0):.2f}, "
            f"Gini={diversity.get('gini', 0.0):.2f}"
        )

    per_dataset = summary.get("datasets", {})
    if per_dataset:
        lines.append("")
        lines.append("### Dataset Highlights")
        lines.append("")
        for dataset_id, payload in sorted(per_dataset.items()):
            lines.append(
                f"- **{dataset_id}**: {payload.get('n_clonotypes', 0):,} clonotypes, "
                f"{payload.get('n_cells', 0):,} cells with receptor calls"
            )
            top = payload.get("top_clonotypes", [])[:3]
            if top:
                formatted = ", ".join(
                    f"{item['clonotype_id']} ({item['cell_count']} cells)" for item in top
                )
                lines.append(f"  - Top clonotypes: {formatted}")

    report_root = Path("processed")
    figures = {
        "clonotype_frequency": result.figures_dir / "clonotype_frequency.png",
        "umap_clonal_expansion": result.figures_dir / "umap_clonal_expansion.png",
        "cdr3_spectratype": result.figures_dir / "cdr3_spectratype.png",
        "vj_usage_heatmap": result.figures_dir / "vj_usage_heatmap.png",
    }

    lines.append("")
    for figure_key in include_figures:
        fig_path = figures.get(figure_key)
        if fig_path and fig_path.exists():
            rel = _rel_path(fig_path, report_root)
            title = figure_key.replace("_", " ").title()
            lines.append(f"![{title}]({rel})")
    lines.append("")

    fragment_path = result.metrics_dir / REPORT_FRAGMENT
    fragment_path.write_text("\n".join(lines))
