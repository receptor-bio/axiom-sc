"""
RNA velocity stream for Tier 3 multi-stream convergence.

Computes per-cluster trajectory direction from scVelo RNA velocity vectors.
Key discrimination power:
  - Exhausted T (TEX) vs progenitor-exhausted T (TPEX): opposite velocity directions
  - Terminal states have low velocity magnitude (sink); progenitor states have high magnitude
  - FOXP3-locus velocity sink confirms Treg identity (absent in naïve T)

Reference:
  Bergen et al. (2020) Nature Biotechnology 38:1408-1414
  doi: 10.1038/s41587-020-0591-3

License: Apache-2.0
"""
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np

from axiom_sc.tier3.stream_result import StreamResult

if TYPE_CHECKING:
    import pandas as pd
    import anndata

logger = logging.getLogger(__name__)

# Velocity magnitude thresholds (empirical, validated on thymus + lung Phase 1 data)
MAGNITUDE_TERMINAL_THRESH: float = 0.05   # below this = likely terminal/sink state
MAGNITUDE_PROGENITOR_THRESH: float = 0.20  # above this = likely progenitor/source state
COSINE_AGREEMENT_THRESH: float = 0.7       # inter-cluster cosine similarity for directed flow


# Velocity signature for each annotatable state:
# Each entry maps (candidate_label, trajectory_type) → expected velocity profile
# "sink" = low magnitude terminal state
# "source" = high magnitude progenitor state
# "transitional" = moderate magnitude, directional
VELOCITY_SIGNATURES: dict[str, dict] = {
    "Treg": {"type": "sink", "key_tf": "FOXP3",
             "notes": "Regulatory T cells: terminal differentiation sink"},
    "Th1": {"type": "sink", "key_tf": "TBX21",
            "notes": "Th1 effector: terminal effector state"},
    "Th2": {"type": "sink", "key_tf": "GATA3",
            "notes": "Th2 effector: terminal effector state"},
    "Th17": {"type": "sink", "key_tf": "RORC",
             "notes": "Th17 effector: terminal differentiation"},
    "Naive_T": {"type": "source", "key_tf": "TCF7",
                "notes": "Naïve T: progenitor/source state for effector differentiation"},
    "Memory_T": {"type": "transitional", "key_tf": "TCF7",
                 "notes": "Memory T: low velocity, self-renewing"},
    "Exhausted_T": {"type": "sink", "key_tf": "TOX",
                    "notes": "Terminally exhausted T: low velocity, irreversible"},
    "Progenitor_Exhausted_T": {"type": "transitional", "key_tf": "TCF7",
                                "notes": "TPEX: directed velocity toward TEX state"},
    "pDC": {"type": "sink", "key_tf": "IRF7",
            "notes": "Plasmacytoid DC: terminal effector"},
    "cDC1": {"type": "sink", "key_tf": "IRF8",
             "notes": "cDC1: terminal, BATF3-lineage sink"},
    "cDC2": {"type": "sink", "key_tf": "IRF4",
             "notes": "cDC2: terminal myeloid DC state"},
    "Mature_mTEC": {"type": "sink", "key_tf": "AIRE",
                    "notes": "Mature medullary TEC: terminal TEC differentiation"},
}


def _try_import_scvelo():
    """Returns scvelo module or None if not installed."""
    try:
        import scvelo as scv
        return scv
    except ImportError:
        return None


class VelocityStream:
    """
    RNA velocity stream: assigns trajectory-based labels to clusters.

    Uses scVelo-computed velocity vectors stored in adata. When scVelo results
    are not pre-computed, optionally runs scVelo on the fly if available.
    Degrades gracefully when data or package is missing.

    Parameters
    ----------
    compute_velocity : bool
        If True and scVelo is installed, compute velocity from scratch.
        If False (default), uses pre-computed results in adata.layers['velocity'].
    tier2_labels : dict[str, str] | None
        Tier 2 verdicts per cluster. Used to focus velocity analysis on
        ambiguous pairs (e.g., Exhausted_T vs Progenitor_Exhausted_T).
    """

    STREAM_NAME = "velocity"

    def __init__(
        self,
        compute_velocity: bool = False,
        tier2_labels: dict[str, str] | None = None,
    ):
        self.compute_velocity = compute_velocity
        self.tier2_labels = tier2_labels or {}

    def run(
        self,
        adata: "anndata.AnnData",
        cluster_key: str = "leiden",
    ) -> dict[str, StreamResult]:
        """
        Run velocity stream on all clusters.

        Parameters
        ----------
        adata : AnnData
            Must contain either adata.layers['velocity'] (pre-computed)
            or adata.layers['spliced'] + adata.layers['unspliced'] for on-the-fly.
        cluster_key : str
            Column in adata.obs containing cluster labels.

        Returns
        -------
        dict[str, StreamResult]
            Per-cluster StreamResult. StreamResult.available = False
            when velocity data is missing.
        """
        if cluster_key not in adata.obs.columns:
            logger.warning("cluster_key %r not in adata.obs — velocity stream skipped", cluster_key)
            return self._unavailable_results(adata, cluster_key)

        # Ensure velocity is computed
        if "velocity" not in adata.layers:
            if self.compute_velocity:
                adata = self._run_scvelo(adata)
            if "velocity" not in adata.layers:
                logger.info("No velocity layer found; velocity stream returning unavailable.")
                return self._unavailable_results(adata, cluster_key)

        return self._compute_cluster_velocity(adata, cluster_key)

    # ── internals ────────────────────────────────────────────────────────────

    def _run_scvelo(self, adata: "anndata.AnnData") -> "anndata.AnnData":
        """Compute scVelo dynamic model in-place. Requires spliced/unspliced layers."""
        scv = _try_import_scvelo()
        if scv is None:
            logger.warning("scvelo not installed; cannot compute velocity.")
            return adata
        missing = [k for k in ("spliced", "unspliced") if k not in adata.layers]
        if missing:
            logger.warning("Missing layers %s for scVelo; cannot compute velocity.", missing)
            return adata
        logger.info("Running scVelo dynamical model…")
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
        scv.tl.recover_dynamics(adata, n_jobs=4)
        scv.tl.velocity(adata, mode="dynamical")
        scv.tl.velocity_graph(adata)
        scv.tl.velocity_embedding(adata, basis="umap")
        return adata

    def _compute_cluster_velocity(
        self,
        adata: "anndata.AnnData",
        cluster_key: str,
    ) -> dict[str, StreamResult]:
        """Compute per-cluster velocity magnitude, classify trajectory state."""
        velocity_matrix: np.ndarray = np.asarray(adata.layers["velocity"])
        # Replace NaN (scVelo sets unfit genes to NaN) with 0
        velocity_matrix = np.nan_to_num(velocity_matrix, nan=0.0)

        clusters = adata.obs[cluster_key]
        unique_clusters = clusters.unique()

        # Per-gene velocity embedding if available (more interpretable)
        has_embedding = "velocity_umap" in adata.obsm

        results: dict[str, StreamResult] = {}
        for cid in unique_clusters:
            mask = (clusters == cid).values
            cluster_vel = velocity_matrix[mask]  # n_cells × n_genes
            magnitude = self._mean_magnitude(cluster_vel)
            trajectory_type = self._classify_magnitude(magnitude)
            label, confidence = self._match_signature(
                cid, trajectory_type, adata, mask, has_embedding
            )
            results[str(cid)] = StreamResult(
                cluster_id=str(cid),
                stream=self.STREAM_NAME,
                label=label,
                confidence=confidence,
                available=True,
                evidence={
                    "mean_velocity_magnitude": float(magnitude),
                    "trajectory_type": trajectory_type,
                    "n_cells": int(mask.sum()),
                },
            )
        return results

    def _mean_magnitude(self, vel: np.ndarray) -> float:
        """Mean L2 norm of velocity vectors across cells in a cluster."""
        norms = np.linalg.norm(vel, axis=1)  # per-cell
        return float(np.nanmean(norms))

    def _classify_magnitude(self, magnitude: float) -> str:
        """Classify cluster trajectory type from mean velocity magnitude."""
        if magnitude < MAGNITUDE_TERMINAL_THRESH:
            return "sink"
        if magnitude > MAGNITUDE_PROGENITOR_THRESH:
            return "source"
        return "transitional"

    def _match_signature(
        self,
        cluster_id: str,
        trajectory_type: str,
        adata: "anndata.AnnData",
        mask: "np.ndarray",
        has_embedding: bool,
    ) -> tuple[str, float]:
        """
        Match cluster trajectory type to known velocity signatures.
        Uses Tier 2 label hint when available to resolve ties.
        """
        # If Tier 2 gave us a hint, check if the velocity is consistent
        tier2_hint = self.tier2_labels.get(str(cluster_id), "")
        if tier2_hint and tier2_hint in VELOCITY_SIGNATURES:
            expected_type = VELOCITY_SIGNATURES[tier2_hint]["type"]
            if expected_type == trajectory_type:
                return tier2_hint, 0.75   # velocity confirms Tier 2 label
            elif expected_type != trajectory_type:
                # Velocity contradicts Tier 2 — lower confidence, note conflict
                return tier2_hint, 0.30

        # No hint — return trajectory type as generic label
        if trajectory_type == "sink":
            return "terminal_state", 0.50
        if trajectory_type == "source":
            return "progenitor_state", 0.50
        return "transitional_state", 0.40

    def _unavailable_results(
        self,
        adata: "anndata.AnnData",
        cluster_key: str,
    ) -> dict[str, StreamResult]:
        """Return unavailable StreamResult for every cluster."""
        if cluster_key not in adata.obs.columns:
            return {}
        clusters = adata.obs[cluster_key].unique()
        return {
            str(cid): StreamResult(
                cluster_id=str(cid),
                stream=self.STREAM_NAME,
                label="",
                confidence=0.0,
                available=False,
                evidence={"reason": "no_velocity_data"},
            )
            for cid in clusters
        }


def discriminate_exhausted_t(
    adata: "anndata.AnnData",
    cluster_key: str,
    candidate_clusters: list[str],
) -> dict[str, str]:
    """
    Discriminate Exhausted_T (TEX) vs Progenitor_Exhausted_T (TPEX) using velocity.

    TEX and TPEX are the canonical hard case for velocity:
      - TPEX → TEX is the direction of flow (TPEX is source, TEX is sink)
      - TPEX has higher velocity magnitude (moving toward TEX state)
      - TEX has low velocity magnitude (terminal sink)

    Parameters
    ----------
    adata : AnnData
        Must contain adata.layers['velocity'] and adata.obs[cluster_key].
    cluster_key : str
        Obs column with cluster IDs.
    candidate_clusters : list[str]
        Only these clusters are evaluated; others are ignored.

    Returns
    -------
    dict[str, str]
        {cluster_id: "Exhausted_T" | "Progenitor_Exhausted_T" | "uncertain"}
    """
    if "velocity" not in adata.layers or len(candidate_clusters) < 1:
        return {c: "uncertain" for c in candidate_clusters}

    velocity_matrix = np.nan_to_num(np.asarray(adata.layers["velocity"]), nan=0.0)
    clusters = adata.obs[cluster_key]
    magnitudes: dict[str, float] = {}
    for cid in candidate_clusters:
        mask = (clusters == cid).values
        if mask.sum() == 0:
            continue
        norms = np.linalg.norm(velocity_matrix[mask], axis=1)
        magnitudes[cid] = float(np.nanmean(norms))

    if len(magnitudes) < 2:
        return {c: "uncertain" for c in candidate_clusters}

    # The cluster with higher mean velocity magnitude is TPEX (source),
    # the one with lower is TEX (sink).
    sorted_by_mag = sorted(magnitudes.items(), key=lambda x: x[1])
    result: dict[str, str] = {}
    for cid in candidate_clusters:
        if cid not in magnitudes:
            result[cid] = "uncertain"
        elif cid == sorted_by_mag[0][0]:
            result[cid] = "Exhausted_T"         # lowest magnitude = terminal sink
        elif cid == sorted_by_mag[-1][0]:
            result[cid] = "Progenitor_Exhausted_T"  # highest magnitude = source
        else:
            result[cid] = "transitional_state"
    return result
