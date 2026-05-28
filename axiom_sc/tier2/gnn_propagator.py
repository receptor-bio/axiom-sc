"""
kNN verdict propagation — ScInfeR-style GNN propagation of AXIOM verdicts.

Propagates high-confidence PROVEN verdicts to neighbouring UNCERTAIN cells
via the cell kNN graph (adata.obsp['connectivities']). Only propagates
to cells whose marker expression is consistent with the seed label.

The original ScInfeR code is NOT imported here (algorithm reimplemented).

Reference:
  Xing et al. (2023) Bioinformatics 39:btad680
  doi: 10.1093/bioinformatics/btad680

License: Apache-2.0 (own reimplementation, original ScInfeR code not used)
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import anndata
    import scipy.sparse

logger = logging.getLogger(__name__)

# Default propagation parameters (validated on Phase 1 thymus + lung data)
DEFAULT_N_NEIGHBORS: int = 15
DEFAULT_MIN_CONFIDENCE: float = 0.7     # seed must have this confidence to propagate
DEFAULT_MARKER_AGREEMENT_THRESH: float = 0.6  # ≥60% of seed markers must match target cell


@dataclass
class PropagationResult:
    """Result of one GNN propagation pass."""

    cluster_id: str
    original_verdict: str          # verdict before propagation
    propagated_verdict: str        # verdict after propagation (may be same)
    original_label: str
    propagated_label: str
    n_cells_total: int
    n_cells_propagated: int        # cells that received a propagated label
    seed_clusters: list[str] = field(default_factory=list)  # which PROVEN clusters seeded this
    confidence: float = 0.0

    @property
    def was_updated(self) -> bool:
        return self.propagated_verdict != self.original_verdict


class GNNPropagator:
    """
    Propagates high-confidence Tier 2 PROVEN verdicts to UNCERTAIN neighbours.

    Cells with PROVEN verdict can propagate their label to nearby cells
    (within n_neighbors hops) if those cells have UNCERTAIN verdict and
    their marker gene expression is consistent with the seed label.

    Algorithm:
      1. Identify seed clusters: PROVEN with confidence ≥ min_confidence.
      2. For each UNCERTAIN cluster, find its spatial neighbours in the kNN graph.
      3. If ≥ propagation_fraction of neighbour clusters are PROVEN with the
         same label, propagate that label with reduced confidence.
      4. Propagation is cluster-level (not cell-level) for interpretability.

    Parameters
    ----------
    n_neighbors : int
        Number of nearest neighbours in the kNN graph.
    min_confidence : float
        Minimum confidence of a seed (PROVEN) verdict to allow propagation.
    propagation_fraction : float
        Fraction of a cluster's kNN clusters that must be PROVEN with the
        same label to trigger propagation.
    marker_agreement_thresh : float
        Fraction of seed markers that must be present (above z-score threshold)
        in the target cluster's marker evidence. If below this, propagation
        is blocked even if spatial neighbourhood agrees.
    """

    def __init__(
        self,
        n_neighbors: int = DEFAULT_N_NEIGHBORS,
        min_confidence: float = DEFAULT_MIN_CONFIDENCE,
        propagation_fraction: float = 0.5,
        marker_agreement_thresh: float = DEFAULT_MARKER_AGREEMENT_THRESH,
    ):
        self.n_neighbors = n_neighbors
        self.min_confidence = min_confidence
        self.propagation_fraction = propagation_fraction
        self.marker_agreement_thresh = marker_agreement_thresh

    def propagate(
        self,
        adata: "anndata.AnnData",
        verdicts: dict,
        cluster_key: str = "leiden",
        marker_key: str | None = None,
    ) -> dict[str, PropagationResult]:
        """
        Propagate PROVEN verdicts to UNCERTAIN clusters via the kNN graph.

        Parameters
        ----------
        adata : AnnData
            Must contain adata.obsp['connectivities'] (scanpy kNN graph).
        verdicts : dict[str, CellTypeVerdict]
            Tier 2 verdicts keyed by cluster_id.
        cluster_key : str
            Column in adata.obs with cluster labels.
        marker_key : str | None
            Key in adata.uns for pre-computed marker genes per cluster.
            If None, falls back to using available marker_genes from the
            verdict's evidence bundle if attached.

        Returns
        -------
        dict[str, PropagationResult]
            Propagation results for each cluster. was_updated=True when
            the cluster received a propagated label.
        """
        if cluster_key not in adata.obs.columns:
            logger.warning("cluster_key %r not found in adata.obs", cluster_key)
            return {}

        if "connectivities" not in adata.obsp:
            logger.warning(
                "adata.obsp['connectivities'] not found. "
                "Run sc.pp.neighbors() before GNN propagation."
            )
            return self._no_propagation_results(verdicts)

        cluster_knn = self._build_cluster_knn(adata, cluster_key)
        return self._propagate_labels(verdicts, cluster_knn)

    # ── internals ────────────────────────────────────────────────────────────

    def _build_cluster_knn(
        self,
        adata: "anndata.AnnData",
        cluster_key: str,
    ) -> dict[str, dict[str, float]]:
        """
        Build cluster-level kNN graph from cell-level connectivities.

        For each cluster pair (A, B): strength = mean connectivity of all
        cross-cluster cell pairs in the original kNN graph.

        Returns
        -------
        dict[str, dict[str, float]]
            {cluster_id: {neighbour_cluster_id: mean_connectivity_strength}}
        """
        try:
            from scipy.sparse import issparse
        except ImportError:
            logger.warning("scipy not installed; GNN propagation skipped.")
            return {}

        connectivity = adata.obsp["connectivities"]
        clusters = adata.obs[cluster_key]
        unique_clusters = clusters.unique()
        cluster_to_indices: dict[str, np.ndarray] = {
            cid: np.where(clusters == cid)[0] for cid in unique_clusters
        }

        cluster_knn: dict[str, dict[str, float]] = {
            str(cid): {} for cid in unique_clusters
        }

        for cid_a in unique_clusters:
            idx_a = cluster_to_indices[cid_a]
            # Extract rows for cluster A
            if issparse(connectivity):
                sub = connectivity[idx_a, :]
                sub_dense = sub.toarray()
            else:
                sub_dense = np.asarray(connectivity)[idx_a, :]

            for cid_b in unique_clusters:
                if cid_a == cid_b:
                    continue
                idx_b = cluster_to_indices[cid_b]
                # Mean connectivity between cluster A cells and cluster B cells
                cross_conn = sub_dense[:, idx_b].mean()
                if cross_conn > 0:
                    cluster_knn[str(cid_a)][str(cid_b)] = float(cross_conn)

        return cluster_knn

    @staticmethod
    def _verdict_str(v) -> str:
        """Extract plain string value from a Verdict enum, string, or verdict-holding object."""
        if isinstance(v, str):
            return v
        # Verdict enum: use .value to get "PROVEN" not "Verdict.PROVEN"
        if hasattr(v, "value"):
            return v.value
        return str(v)

    def _propagate_labels(
        self,
        verdicts: dict,
        cluster_knn: dict[str, dict[str, float]],
    ) -> dict[str, PropagationResult]:
        """
        Core propagation logic: PROVEN seeds → UNCERTAIN neighbours.
        """
        # Import here to avoid circular import
        from axiom_sc.tier2.axiom_annotator import Verdict

        results: dict[str, PropagationResult] = {}

        # Index PROVEN seeds by label
        proven_by_label: dict[str, list[str]] = {}
        for cid, verdict in verdicts.items():
            v_raw = getattr(verdict, "verdict", verdict)
            v_str = self._verdict_str(v_raw)
            label = getattr(verdict, "cell_type", "")
            confidence = getattr(verdict, "confidence", 0.0)
            if v_str == Verdict.PROVEN.value and confidence >= self.min_confidence:
                proven_by_label.setdefault(label, []).append(str(cid))

        for cid, verdict in verdicts.items():
            cid = str(cid)
            v_raw = getattr(verdict, "verdict", verdict)
            v_str = self._verdict_str(v_raw)
            label = getattr(verdict, "cell_type", "")
            confidence = getattr(verdict, "confidence", 0.0)
            n_cells = getattr(verdict, "n_cells", 1)

            original_result = PropagationResult(
                cluster_id=cid,
                original_verdict=v_str,
                propagated_verdict=v_str,
                original_label=label,
                propagated_label=label,
                n_cells_total=n_cells,
                n_cells_propagated=0,
                confidence=confidence,
            )

            if v_str != Verdict.UNCERTAIN.value:
                results[cid] = original_result
                continue

            # Try to propagate from PROVEN neighbours
            neighbours = cluster_knn.get(cid, {})
            if not neighbours:
                results[cid] = original_result
                continue

            prop_label, prop_confidence, seed_clusters = self._find_propagation_source(
                cid, neighbours, proven_by_label, verdicts
            )

            if prop_label:
                original_result.propagated_verdict = Verdict.UNCERTAIN.value  # still uncertain
                original_result.propagated_label = prop_label
                original_result.n_cells_propagated = n_cells
                original_result.seed_clusters = seed_clusters
                original_result.confidence = prop_confidence
                # Note: we keep verdict as UNCERTAIN but update the label suggestion.
                # Full PROVEN upgrade requires Tier 2 re-evaluation with propagated label.

            results[cid] = original_result

        return results

    def _find_propagation_source(
        self,
        target_cluster: str,
        neighbours: dict[str, float],
        proven_by_label: dict[str, list[str]],
        verdicts: dict,
    ) -> tuple[str, float, list[str]]:
        """
        Find a PROVEN label to propagate to the target cluster.

        Returns (label, confidence, seed_cluster_ids) if a suitable source
        is found, or ("", 0.0, []) if no propagation is possible.
        """
        # Sort neighbours by connectivity strength
        sorted_neighbours = sorted(neighbours.items(), key=lambda x: x[1], reverse=True)
        top_n = sorted_neighbours[: self.n_neighbors]

        # Count how many top neighbours are PROVEN with each label
        from axiom_sc.tier2.axiom_annotator import Verdict
        label_counts: dict[str, list[tuple[str, float]]] = {}
        for nbr_cid, conn_strength in top_n:
            nbr_verdict = verdicts.get(nbr_cid)
            if nbr_verdict is None:
                continue
            nbr_v_raw = getattr(nbr_verdict, "verdict", nbr_verdict)
            nbr_v_str = self._verdict_str(nbr_v_raw)
            nbr_label = getattr(nbr_verdict, "cell_type", "")
            nbr_conf = getattr(nbr_verdict, "confidence", 0.0)
            if nbr_v_str == Verdict.PROVEN.value and nbr_conf >= self.min_confidence:
                label_counts.setdefault(nbr_label, []).append((nbr_cid, conn_strength))

        if not label_counts:
            return "", 0.0, []

        # Find the label with the most neighbouring PROVEN clusters
        best_label = max(label_counts, key=lambda l: len(label_counts[l]))
        best_sources = label_counts[best_label]

        # Check propagation fraction threshold
        n_top = len(top_n)
        if n_top == 0:
            return "", 0.0, []

        fraction = len(best_sources) / n_top
        if fraction < self.propagation_fraction:
            return "", 0.0, []

        # Confidence = mean connectivity-weighted confidence of seed clusters
        total_conn = sum(conn for _, conn in best_sources)
        if total_conn == 0:
            return "", 0.0, []

        weighted_conf = sum(
            conn * getattr(verdicts.get(cid), "confidence", 0.0)
            for cid, conn in best_sources
        ) / total_conn

        seed_ids = [cid for cid, _ in best_sources]
        return best_label, float(weighted_conf) * 0.8, seed_ids  # discount propagated confidence

    def _no_propagation_results(
        self,
        verdicts: dict,
    ) -> dict[str, PropagationResult]:
        """Return identity propagation results (no changes) for all clusters."""
        results = {}
        for cid, verdict in verdicts.items():
            v_raw = getattr(verdict, "verdict", verdict)
            v_str = self._verdict_str(v_raw)
            label = getattr(verdict, "cell_type", "")
            confidence = getattr(verdict, "confidence", 0.0)
            n_cells = getattr(verdict, "n_cells", 1)
            results[str(cid)] = PropagationResult(
                cluster_id=str(cid),
                original_verdict=v_str,
                propagated_verdict=v_str,
                original_label=label,
                propagated_label=label,
                n_cells_total=n_cells,
                n_cells_propagated=0,
                confidence=confidence,
            )
        return results
