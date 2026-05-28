"""
Per-regulon z-score normalization of pySCENIC AUC values.

Phase 1 lesson: raw AUC values are not comparable across runs or datasets.
Thymus runs produced AUC values in range 0–12; lung in 0.01–0.35.
Z-scoring per regulon across all cells in a single run produces
comparable activity scores that are threshold-able at a fixed value.

Phase 1 validated thresholds (hardcoded — do not change):
  z ≥ 1.5  → regulon active
  z ≤ -0.5 → regulon inactive

License: Apache 2.0
"""
from __future__ import annotations

import numpy as np


def zscore_auc_matrix(
    auc_matrix: np.ndarray,
    regulon_names: list[str],
    ddof: int = 1,
) -> dict[str, np.ndarray]:
    """
    Z-score an AUC matrix (cells × regulons) per regulon.

    Parameters
    ----------
    auc_matrix : np.ndarray
        Shape (n_cells, n_regulons). Raw AUC values from pySCENIC AUCell.
    regulon_names : list[str]
        Regulon names corresponding to matrix columns.
    ddof : int
        Delta degrees of freedom for std (default 1 for sample std).

    Returns
    -------
    dict[regulon_name -> np.ndarray of z-scores per cell]
    """
    if auc_matrix.ndim != 2:
        raise ValueError(f"auc_matrix must be 2D, got shape {auc_matrix.shape}")
    if auc_matrix.shape[1] != len(regulon_names):
        raise ValueError(
            f"auc_matrix has {auc_matrix.shape[1]} columns "
            f"but {len(regulon_names)} regulon names provided"
        )

    result: dict[str, np.ndarray] = {}
    for i, regulon in enumerate(regulon_names):
        col = auc_matrix[:, i].astype(float)
        mean = col.mean()
        std = col.std(ddof=ddof)
        if std < 1e-10:
            # Constant regulon (no variation) — all cells get z=0
            result[regulon] = np.zeros_like(col)
        else:
            result[regulon] = (col - mean) / std
    return result


def zscore_auc_dict(
    auc_dict: dict[str, list[float]],
    ddof: int = 1,
) -> dict[str, np.ndarray]:
    """
    Z-score a dict[regulon -> list[auc_per_cell]] per regulon.

    Convenience wrapper around zscore_auc_matrix for the dict format
    returned by scenic_worker.py.
    """
    if not auc_dict:
        return {}
    regulon_names = list(auc_dict.keys())
    matrix = np.array([auc_dict[r] for r in regulon_names], dtype=float).T  # (n_cells, n_regulons)
    return zscore_auc_matrix(matrix, regulon_names, ddof=ddof)


def cluster_mean_zscore(
    zscore_dict: dict[str, np.ndarray],
    cell_indices: list[int],
) -> dict[str, float]:
    """
    Compute mean z-score for each regulon across a set of cells (a cluster).

    Parameters
    ----------
    zscore_dict : dict[regulon -> np.ndarray of per-cell z-scores]
        Output of zscore_auc_dict or zscore_auc_matrix.
    cell_indices : list[int]
        Indices of cells belonging to this cluster.

    Returns
    -------
    dict[regulon -> mean_z_score] — suitable for EvidenceBundle.regulons
    """
    result: dict[str, float] = {}
    for regulon, zscores in zscore_dict.items():
        vals = zscores[cell_indices]
        result[regulon] = float(vals.mean())
    return result
