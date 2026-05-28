# License: Apache 2.0
"""
Preprocessing pipeline for AXIOM-SC.

Prepares AnnData for annotation based on InputInspection results.
Runs only the steps marked in inspection.preprocessing_needed.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from axiom_sc.input.detector import InputInspection


def prepare_for_annotation(adata, inspection: "InputInspection"):
    """
    Applies preprocessing steps listed in inspection.preprocessing_needed.
    Returns adata ready for tier annotation (normalized, log1p, HVGs selected).

    Input adata is modified in-place and returned.
    """
    import numpy as np

    steps = set(inspection.preprocessing_needed)

    if "10k_normalization" in steps:
        _normalize_10k(adata)

    if "log1p" in steps:
        import scanpy as sc
        sc.pp.log1p(adata)

    if "hvg_selection" in steps:
        _select_hvg(adata)

    if "pca" in steps:
        import scanpy as sc
        sc.tl.pca(adata, svd_solver="arpack", n_comps=min(50, adata.n_vars - 1))

    if "knn_graph" in steps:
        import scanpy as sc
        sc.pp.neighbors(adata, n_neighbors=15)

    if "leiden_clustering" in steps:
        import scanpy as sc
        sc.tl.leiden(adata, resolution=0.5, key_added="leiden")

    if "de_markers_per_cluster" in steps:
        _compute_de_markers(adata, inspection.cluster_key or "leiden")

    return adata


def _normalize_10k(adata) -> None:
    import scipy.sparse as sp
    import numpy as np

    X = adata.X
    if sp.issparse(X):
        totals = np.asarray(X.sum(axis=1)).flatten()
    else:
        totals = X.sum(axis=1)

    scale = 10_000 / np.where(totals == 0, 1, totals)

    if sp.issparse(X):
        from scipy.sparse import diags
        adata.X = diags(scale) @ X
    else:
        adata.X = X * scale[:, None]


def _select_hvg(adata, n_top: int = 2000) -> None:
    import scanpy as sc
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top, subset=True)


def _compute_de_markers(adata, cluster_key: str) -> None:
    """Runs scanpy rank_genes_groups for Tier 2 marker input."""
    import scanpy as sc
    if cluster_key not in adata.obs.columns:
        return
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method="wilcoxon",
                             key_added="rank_genes_groups")
