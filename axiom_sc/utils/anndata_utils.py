"""
AnnData utility functions shared across tiers.

License: Apache 2.0
"""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import anndata


def get_cluster_marker_zscores(
    adata: "anndata.AnnData",
    cluster_key: str = "leiden",
    n_genes: int = 50,
) -> dict[str, dict[str, float]]:
    """
    Extract per-cluster marker gene z-scores from an AnnData object.

    Expects adata.uns[f"rank_genes_groups_{cluster_key}"] to exist
    (run sc.tl.rank_genes_groups first).

    Returns
    -------
    dict[cluster_id -> dict[gene -> z_score]]
    """
    try:
        import scanpy as sc
        import numpy as np
    except ImportError:
        raise ImportError("scanpy is required for get_cluster_marker_zscores. "
                          "Install with: pip install 'axiom-sc[science]'")

    key = "rank_genes_groups"
    if key not in adata.uns:
        raise ValueError(
            f"adata.uns['{key}'] not found. "
            "Run sc.tl.rank_genes_groups(adata) first."
        )

    result: dict[str, dict[str, float]] = {}
    rgg = adata.uns[key]
    groups = rgg["names"].dtype.names

    for group in groups:
        genes = list(rgg["names"][group][:n_genes])
        scores = list(rgg["scores"][group][:n_genes])
        result[str(group)] = dict(zip(genes, [float(s) for s in scores]))

    return result


def validate_adata_normalized(adata: "anndata.AnnData") -> None:
    """
    Verify the AnnData is 10k-normalized and log1p-transformed.
    Raises ValueError if suspicious.
    """
    import numpy as np

    X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()

    sample = X[:100, :100]
    max_val = float(sample.max())

    if max_val > 20:
        raise ValueError(
            f"adata.X max value {max_val:.1f} looks like raw counts, not log1p-normalized. "
            "Run sc.pp.normalize_total(adata, target_sum=1e4) then sc.pp.log1p(adata)."
        )
