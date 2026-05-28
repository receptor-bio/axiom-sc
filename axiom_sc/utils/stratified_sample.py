"""
Stratified sampling of AnnData cells ensuring minimum coverage per cell type.

Phase 1 lesson: 20k cells total is insufficient to recover rare TF regulons
(FOXP3, AIRE) because the rare cell types contribute too few cells for GRN
training. Minimum 300 cells per type ensures signal in grnboost2.

Target: 50k total cells (recommended_total in SCENIC_DEFAULTS), with at least
300 per type. Rare types (<300 cells) are included in full.

License: Apache 2.0
"""
from __future__ import annotations

import numpy as np
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import anndata


def stratified_sample(
    adata: "anndata.AnnData",
    cell_type_key: str,
    target_total: int = 50_000,
    min_per_type: int = 300,
    seed: int = 42,
) -> "anndata.AnnData":
    """
    Stratified subsample of an AnnData object.

    Ensures at least `min_per_type` cells per cell type (or all available
    cells if a type has fewer). Total cell count ≤ `target_total`.

    Parameters
    ----------
    adata : AnnData
        Input AnnData with cell type labels.
    cell_type_key : str
        Key in adata.obs with cell type labels.
    target_total : int
        Target total number of cells after sampling.
    min_per_type : int
        Minimum cells to include per cell type.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    AnnData
        Subsampled AnnData with ≤ target_total cells.
    """
    try:
        import anndata as ad
    except ImportError:
        raise ImportError("anndata required. Install with: pip install 'axiom-sc[science]'")

    if cell_type_key not in adata.obs.columns:
        raise ValueError(f"Cell type key {cell_type_key!r} not found in adata.obs")

    rng = np.random.default_rng(seed)
    cell_types = adata.obs[cell_type_key].unique()
    selected_indices: list[int] = []

    # First pass: guarantee minimums
    type_budgets: dict[str, int] = {}
    guaranteed_total = 0
    for ct in cell_types:
        ct_mask = adata.obs[cell_type_key] == ct
        ct_count = int(ct_mask.sum())
        guaranteed = min(ct_count, min_per_type)
        type_budgets[ct] = guaranteed
        guaranteed_total += guaranteed

    # Second pass: distribute remaining budget proportionally
    remaining = max(0, target_total - guaranteed_total)
    if remaining > 0:
        for ct in cell_types:
            ct_mask = adata.obs[cell_type_key] == ct
            ct_count = int(ct_mask.sum())
            available = ct_count - type_budgets[ct]
            if available > 0:
                # Proportional allocation based on type size
                extra = int(remaining * available / max(1, sum(
                    max(0, int((adata.obs[cell_type_key] == c).sum()) - type_budgets[c])
                    for c in cell_types
                )))
                type_budgets[ct] += min(extra, available)

    # Sample cells per type
    for ct in cell_types:
        ct_indices = np.where(adata.obs[cell_type_key] == ct)[0]
        n_sample = min(type_budgets[ct], len(ct_indices))
        chosen = rng.choice(ct_indices, size=n_sample, replace=False)
        selected_indices.extend(chosen.tolist())

    selected_indices.sort()
    return adata[selected_indices].copy()


def get_type_counts(adata: "anndata.AnnData", cell_type_key: str) -> dict[str, int]:
    """Return per-cell-type cell counts."""
    return dict(adata.obs[cell_type_key].value_counts())


def check_min_cells(
    adata: "anndata.AnnData",
    cell_type_key: str,
    min_per_type: int = 300,
) -> dict[str, int]:
    """
    Return cell types that have fewer than min_per_type cells.
    Empty dict means all types meet the minimum.
    """
    counts = get_type_counts(adata, cell_type_key)
    return {ct: n for ct, n in counts.items() if n < min_per_type}
