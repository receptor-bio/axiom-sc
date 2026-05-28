"""
scType scoring algorithm reimplemented in Python.

The original scType R code is GPL v3 — NOT imported here.
This file reimplements the scoring algorithm in Python (~200 lines).

Reference:
  Ianevski et al. (2022) Nature Communications 13:1246
  doi: 10.1038/s41467-022-28803-w

License: Apache 2.0 (own reimplementation, original R code not used)
"""
from __future__ import annotations

import numpy as np
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd


def score_cell_types(
    expr_matrix: "np.ndarray",
    gene_names: list[str],
    cell_type_markers: dict[str, list[str]],
    cell_type_negative: dict[str, list[str]] | None = None,
) -> "pd.DataFrame":
    """
    Score each cell for each cell type using positive and negative marker genes.

    Reimplements the scType scoring formula:
      score(c, t) = sum(log1p(expr[g, c])) for g in positive_markers[t]
                  - sum(log1p(expr[g, c])) for g in negative_markers[t]

    Parameters
    ----------
    expr_matrix : np.ndarray
        Cells × genes expression matrix (log1p-normalized).
    gene_names : list[str]
        Gene names corresponding to matrix columns.
    cell_type_markers : dict[str, list[str]]
        Positive marker genes per cell type.
    cell_type_negative : dict[str, list[str]] | None
        Negative marker genes per cell type (optional).

    Returns
    -------
    pd.DataFrame
        Cells × cell_types score matrix.
    """
    try:
        import pandas as pd
    except ImportError:
        raise ImportError("pandas required. Install with: pip install 'axiom-sc[science]'")

    gene_index = {g: i for i, g in enumerate(gene_names)}
    cell_type_negative = cell_type_negative or {}
    results = {}

    for ct, pos_genes in cell_type_markers.items():
        neg_genes = cell_type_negative.get(ct, [])

        pos_idx = [gene_index[g] for g in pos_genes if g in gene_index]
        neg_idx = [gene_index[g] for g in neg_genes if g in gene_index]

        pos_score = expr_matrix[:, pos_idx].sum(axis=1) if pos_idx else np.zeros(expr_matrix.shape[0])
        neg_score = expr_matrix[:, neg_idx].sum(axis=1) if neg_idx else np.zeros(expr_matrix.shape[0])

        results[ct] = pos_score - neg_score

    return pd.DataFrame(results)
