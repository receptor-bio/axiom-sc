# License: Apache 2.0
"""
Input format detection for AXIOM-SC.

Detects input type, modalities, and preprocessing needed in <2s.
Result is shown in the UI's InspectionPanel before the annotation job starts.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Optional


class InputType(str, Enum):
    CLUSTERED_H5AD    = "clustered_h5ad"    # has cluster obs column + normalized X
    UNCLUSTERED_H5AD  = "unclustered_h5ad"  # normalized X, no clusters
    RAW_COUNT_H5AD    = "raw_count_h5ad"    # integer counts, needs full preprocessing
    MARKER_CSV        = "marker_csv"        # cluster DE marker table → skip to Tier 2
    CELLRANGER_H5     = "cellranger_h5"     # 10x CellRanger output


@dataclass
class InputInspection:
    input_type: str
    n_cells: int = 0
    n_genes: int = 0
    n_clusters: int = 0
    cluster_key: Optional[str] = None
    has_velocity: bool = False
    has_atac: bool = False
    has_spatial: bool = False
    preprocessing_needed: list[str] = field(default_factory=list)
    estimated_minutes: float = 5.0
    active_streams: int = 4
    warnings: list[str] = field(default_factory=list)


# Cluster obs column candidates (checked in order)
_CLUSTER_KEYS = ["leiden", "louvain", "seurat_clusters", "cell_type",
                 "celltype", "cluster", "cluster_id", "CellType"]


def inspect_input(path: str) -> InputInspection:
    """
    Auto-detects input type and modalities in <2 seconds (backed=r for h5ad).
    Called before queuing the annotation job so the UI can show InspectionPanel.
    """
    p = Path(path)

    if p.suffix == ".csv":
        return _inspect_csv(p)
    if p.suffix == ".h5":
        return _inspect_cellranger_h5(p)
    if p.suffix in (".h5ad",):
        return _inspect_h5ad(p)

    # Fallback — treat as raw h5ad
    return InputInspection(
        input_type=InputType.RAW_COUNT_H5AD,
        warnings=[f"Unknown extension '{p.suffix}', treating as raw h5ad"],
    )


def _inspect_csv(p: Path) -> InputInspection:
    import csv
    with open(p, newline="") as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames or []
        rows = list(reader)

    n_clusters = len({r.get("cluster_id", r.get("cluster", "")) for r in rows}) if rows else 0
    return InputInspection(
        input_type=InputType.MARKER_CSV,
        n_clusters=n_clusters,
        preprocessing_needed=[],
        estimated_minutes=2.0,
        active_streams=0,
        warnings=[] if any(h in headers for h in ["logfoldchanges", "avg_log2FC", "markers"])
                 else ["CSV missing expected column: logfoldchanges / avg_log2FC / markers"],
    )


def _inspect_cellranger_h5(p: Path) -> InputInspection:
    try:
        import h5py
        with h5py.File(p, "r") as f:
            n_cells = f["matrix/barcodes"].shape[0] if "matrix/barcodes" in f else 0
            n_genes = f["matrix/features/id"].shape[0] if "matrix/features/id" in f else 0
    except Exception:
        n_cells = n_genes = 0

    est = max(5.0, n_cells / 10_000)
    return InputInspection(
        input_type=InputType.CELLRANGER_H5,
        n_cells=n_cells,
        n_genes=n_genes,
        preprocessing_needed=["10k_normalization", "log1p", "hvg_selection",
                              "pca", "knn_graph", "leiden_clustering"],
        estimated_minutes=est,
        active_streams=4,
        warnings=[],
    )


def _inspect_h5ad(p: Path) -> InputInspection:
    try:
        import anndata
        import numpy as np
        adata = anndata.read_h5ad(str(p), backed="r")
    except ImportError:
        return InputInspection(
            input_type=InputType.RAW_COUNT_H5AD,
            warnings=["anndata not installed — cannot inspect h5ad fully"],
        )
    except Exception as e:
        return InputInspection(
            input_type=InputType.RAW_COUNT_H5AD,
            warnings=[f"Could not read h5ad: {e}"],
        )

    n_cells, n_genes = adata.n_obs, adata.n_vars

    # Detect cluster key
    cluster_key = next((k for k in _CLUSTER_KEYS if k in adata.obs.columns), None)
    n_clusters = adata.obs[cluster_key].nunique() if cluster_key else 0

    # Detect modalities
    has_velocity = any(k in adata.layers for k in ("velocity", "spliced", "unspliced"))
    has_atac = any(k in adata.obsm for k in ("X_atac", "gene_activities"))
    has_spatial = "spatial" in adata.obsm

    # Detect normalization state (sample first 100 cells to avoid loading full matrix)
    try:
        import scipy.sparse as sp
        X_sample = adata.X[:100]
        if sp.issparse(X_sample):
            X_sample = X_sample.toarray()
        import numpy as np
        max_val = float(np.max(X_sample))
        is_raw = max_val > 50  # raw counts are integers, max usually in hundreds-thousands
    except Exception:
        is_raw = False
        max_val = 0.0

    # Count active streams based on available modalities
    active_streams = 4  # baseline RNA-only (communication, spatial_niche, cross_species, SCENIC)
    if has_velocity:
        active_streams += 1
    if has_atac:
        active_streams += 1

    # Determine input type and preprocessing needed
    preprocessing_needed: list[str] = []
    if is_raw:
        input_type = InputType.RAW_COUNT_H5AD
        preprocessing_needed = ["10k_normalization", "log1p", "hvg_selection", "pca",
                                 "knn_graph", "leiden_clustering"]
        est = max(10.0, n_cells / 5_000)
    elif cluster_key:
        input_type = InputType.CLUSTERED_H5AD
        preprocessing_needed = ["de_markers_per_cluster"]
        est = max(5.0, n_cells / 8_000)
    else:
        input_type = InputType.UNCLUSTERED_H5AD
        preprocessing_needed = ["leiden_clustering", "de_markers_per_cluster"]
        est = max(7.0, n_cells / 7_000)

    warnings: list[str] = []
    if has_velocity:
        warnings.append("Velocity layers detected — scVelo stream will run (Tier 3)")
    if has_spatial:
        warnings.append("Spatial coordinates detected — spatial niche stream enhanced")
    if is_raw:
        warnings.append(f"Raw counts detected (max={max_val:.0f}) — full preprocessing will run")

    return InputInspection(
        input_type=input_type,
        n_cells=n_cells,
        n_genes=n_genes,
        n_clusters=n_clusters,
        cluster_key=cluster_key,
        has_velocity=has_velocity,
        has_atac=has_atac,
        has_spatial=has_spatial,
        preprocessing_needed=preprocessing_needed,
        estimated_minutes=est,
        active_streams=active_streams,
        warnings=warnings,
    )
