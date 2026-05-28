"""
Spatial niche encoder stream for Tier 3 multi-stream convergence.

Two operational modes:
  1. CELLama mode: uses pre-computed sentence embeddings from
     adata.obsm['X_cellama'] and matches against reference cell type embeddings.
  2. Niche composition mode: builds kNN neighborhood composition from
     adata.obsp['connectivities'] + adata.obs[cluster_key] and matches against
     known niche signatures (what cell types are expected to co-localise).

Key discriminations (from CLAUDE.md Section 10):
  - pDC: found in T cell zones → T cells are dominant neighbours
  - Treg: found adjacent to effector T cells (immunosuppressive niche)
  - Tissue-resident Mφ: near sinusoidal endothelium (spatial position relative to vessels)
  - GCB: inside germinal centers → FDC-like and Tfh neighbour context
  - SMC: coat blood vessels → endothelial neighbours
  - ILC2: gut epithelial niche → near tuft/goblet cells

Reference:
  Lv et al. (2025) Advanced Science
  doi: 10.1002/advs.202413514

License: Apache-2.0
"""
from __future__ import annotations

import logging
import pickle
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from axiom_sc.tier3.stream_result import StreamResult

if TYPE_CHECKING:
    import anndata
    import pandas as pd

logger = logging.getLogger(__name__)

# ── CELLama reference embeddings (Patch B) ────────────────────────────────────
_CELLAMA_EMBEDDINGS_PATH = (
    Path.home() / ".axiom_sc" / "resources" / "spatial" / "cellama_reference_embeddings.pkl"
)


def _load_cellama_reference() -> "dict | None":
    """Load CELLama reference embeddings when pre-downloaded. Returns None if unavailable."""
    if _CELLAMA_EMBEDDINGS_PATH.exists():
        try:
            with open(_CELLAMA_EMBEDDINGS_PATH, "rb") as f:
                ref = pickle.load(f)
            logger.debug("Loaded CELLama reference embeddings: %d cell types", len(ref))
            return ref
        except Exception as exc:
            logger.warning("Could not load CELLama embeddings from %s: %s", _CELLAMA_EMBEDDINGS_PATH, exc)
    return None


CELLAMA_REFERENCE: "dict | None" = _load_cellama_reference()

# Niche signatures: expected neighbor cell type fractions per annotated type.
# Format: {target_cell_type: {neighbour_label_keyword: min_fraction_expected}}
# The keyword is checked as a substring of neighbour cluster labels (case-insensitive).
# All conditions are treated as soft (score = fraction of conditions met).
NICHE_SIGNATURES: dict[str, list[tuple[str, float]]] = {
    "pDC": [
        ("T", 0.15),           # T cells dominant in T cell zones
        ("DC", 0.05),          # co-localise with other DCs
    ],
    "cDC1": [
        ("T", 0.1),
        ("NK", 0.05),
    ],
    "Treg": [
        ("T", 0.2),            # Tregs suppress other T cells → T neighbours abundant
        ("CD8", 0.1),
    ],
    "Th1": [
        ("T", 0.15),
    ],
    "GCB": [
        ("B", 0.15),           # GCBs inside B cell follicles
        ("Tfh", 0.05),
    ],
    "Plasma_cell": [
        ("B", 0.1),
    ],
    "NK": [
        ("T", 0.15),
        ("Mono", 0.05),
    ],
    "Tissue_Resident_Mac": [
        ("Endo", 0.1),         # near endothelium (sinusoidal)
    ],
    "SMC": [
        ("Endo", 0.25),        # SMCs coat blood vessels
    ],
    "CAF": [
        ("Fibro", 0.1),
        ("Endo", 0.05),
    ],
    "ILC2": [
        ("Epi", 0.1),          # near epithelial cells in gut/lung
    ],
    "AT1": [
        ("AT2", 0.15),         # AT1 and AT2 co-localise in alveoli
    ],
    "AT2": [
        ("AT1", 0.1),
    ],
    "Mature_mTEC": [
        ("T", 0.2),            # mature mTEC in medulla surrounded by T cells
        ("DC", 0.05),
    ],
    "Cortical_TEC": [
        ("T", 0.25),           # cTEC surrounded by thymocytes undergoing positive selection
    ],
}

# Minimum fraction of a cluster's total neighbours that must be from
# an expected neighbour type to count it as "present"
NICHE_FRACTION_THRESH: float = 0.05


class SpatialNicheStream:
    """
    Spatial niche stream: assigns cell type labels based on neighbourhood context.

    Parameters
    ----------
    cellama_key : str
        Key in adata.obsm containing pre-computed CELLama embeddings.
    cell_type_obs_key : str | None
        Column in adata.obs containing existing cell type annotations for neighbours.
        Required for niche composition mode. If None, falls back to cluster_key.
    tier2_labels : dict[str, str] | None
        Tier 2 label hints per cluster_id.
    """

    STREAM_NAME = "spatial_niche"

    def __init__(
        self,
        cellama_key: str = "X_cellama",
        cell_type_obs_key: str | None = None,
        tier2_labels: dict[str, str] | None = None,
    ):
        self.cellama_key = cellama_key
        self.cell_type_obs_key = cell_type_obs_key
        self.tier2_labels = tier2_labels or {}

    def run(
        self,
        adata: "anndata.AnnData",
        cluster_key: str = "leiden",
    ) -> dict[str, StreamResult]:
        """
        Run spatial niche stream on all clusters.

        Parameters
        ----------
        adata : AnnData
            May contain adata.obsm['X_cellama'] (pre-computed CELLama embeddings).
            Fallback: uses adata.obsp['connectivities'] for niche composition.
        cluster_key : str

        Returns
        -------
        dict[str, StreamResult]
        """
        if cluster_key not in adata.obs.columns:
            logger.warning("cluster_key %r not in adata.obs — spatial_niche stream skipped", cluster_key)
            return {}

        # Mode 1: CELLama embeddings
        if self.cellama_key in adata.obsm:
            return self._score_from_cellama(adata, cluster_key)

        # Mode 2: kNN neighbourhood composition
        if "connectivities" in adata.obsp:
            return self._score_from_niche_composition(adata, cluster_key)

        logger.info("No spatial niche data available (no CELLama, no kNN). Returning unavailable.")
        return self._unavailable_results(adata, cluster_key)

    # ── CELLama mode ──────────────────────────────────────────────────────────

    def _score_from_cellama(
        self,
        adata: "anndata.AnnData",
        cluster_key: str,
    ) -> dict[str, StreamResult]:
        """
        Score clusters using pre-computed CELLama sentence embeddings.
        Computes cluster-mean embeddings and finds nearest reference cluster
        from a predefined reference library.
        """
        embeddings = np.asarray(adata.obsm[self.cellama_key])
        clusters = adata.obs[cluster_key]
        unique_clusters = clusters.unique()
        results: dict[str, StreamResult] = {}

        # Build per-cluster mean embedding
        cluster_embeddings: dict[str, np.ndarray] = {}
        for cid in unique_clusters:
            mask = (clusters == cid).values
            cluster_embeddings[str(cid)] = embeddings[mask].mean(axis=0)

        for cid in unique_clusters:
            cid_str = str(cid)
            tier2_hint = self.tier2_labels.get(cid_str, "")
            emb = cluster_embeddings[cid_str]

            if CELLAMA_REFERENCE is not None:
                # Embedding mode: cosine similarity to reference cell type embeddings
                label, confidence = self._match_cellama_reference(emb, tier2_hint)
                evidence = {"mode": "cellama_embedding_reference", "embedding_dim": emb.shape[0]}
            else:
                # No reference library — return tier2 hint with low confidence
                label = tier2_hint if tier2_hint else "unknown"
                confidence = 0.5 if tier2_hint else 0.1
                evidence = {"mode": "cellama_embedding", "embedding_dim": emb.shape[0]}

            results[cid_str] = StreamResult(
                cluster_id=cid_str,
                stream=self.STREAM_NAME,
                label=label,
                confidence=confidence,
                available=True,
                evidence=evidence,
            )
        return results

    def _match_cellama_reference(
        self,
        emb: "np.ndarray",
        tier2_hint: str,
    ) -> "tuple[str, float]":
        """Cosine similarity against CELLama reference cell type embeddings."""
        emb_norm = emb / (np.linalg.norm(emb) + 1e-10)
        best_label, best_score = tier2_hint or "unknown", 0.0
        for cell_type, ref_emb in CELLAMA_REFERENCE.items():  # type: ignore[union-attr]
            ref = np.asarray(ref_emb)
            ref_norm = ref / (np.linalg.norm(ref) + 1e-10)
            cosine = float(np.dot(emb_norm, ref_norm))
            if cosine > best_score:
                best_score = cosine
                best_label = cell_type
        confidence = min(1.0, max(0.0, best_score))
        return best_label, confidence

    # ── niche composition mode ────────────────────────────────────────────────

    def _score_from_niche_composition(
        self,
        adata: "anndata.AnnData",
        cluster_key: str,
    ) -> dict[str, StreamResult]:
        """
        Score clusters using kNN neighbourhood composition.

        For each cluster, computes what fraction of its kNN neighbours
        belong to each cluster label, then matches against NICHE_SIGNATURES.
        """
        try:
            from scipy.sparse import issparse
        except ImportError:
            return self._unavailable_results(adata, cluster_key)

        # Choose obs column for neighbour labels — prefer explicit cell type col
        label_key = self.cell_type_obs_key or cluster_key
        if label_key not in adata.obs.columns:
            label_key = cluster_key

        clusters = adata.obs[cluster_key]
        neighbour_labels = adata.obs[label_key]
        unique_clusters = clusters.unique()

        connectivity = adata.obsp["connectivities"]
        results: dict[str, StreamResult] = {}

        for cid in unique_clusters:
            cid_str = str(cid)
            cell_mask = (clusters == cid).values
            cell_indices = np.where(cell_mask)[0]

            # Aggregate neighbour labels for all cells in this cluster
            neighbour_type_counts: dict[str, float] = {}
            total_connectivity = 0.0

            for cell_idx in cell_indices:
                if issparse(connectivity):
                    row = connectivity[cell_idx].toarray()[0]
                else:
                    row = np.asarray(connectivity[cell_idx]).flatten()

                for nbr_idx, conn_weight in enumerate(row):
                    if conn_weight == 0 or nbr_idx == cell_idx:
                        continue
                    nbr_label = str(neighbour_labels.iloc[nbr_idx])
                    neighbour_type_counts[nbr_label] = (
                        neighbour_type_counts.get(nbr_label, 0.0) + conn_weight
                    )
                    total_connectivity += conn_weight

            if total_connectivity == 0:
                results[cid_str] = StreamResult(
                    cluster_id=cid_str, stream=self.STREAM_NAME,
                    label="", confidence=0.0, available=True,
                    evidence={"reason": "isolated_cluster"},
                )
                continue

            # Normalise to fractions
            niche_fractions = {
                k: v / total_connectivity
                for k, v in neighbour_type_counts.items()
            }

            label, confidence, evidence = self._match_niche_signature(
                niche_fractions, cid_str
            )
            results[cid_str] = StreamResult(
                cluster_id=cid_str,
                stream=self.STREAM_NAME,
                label=label,
                confidence=confidence,
                available=True,
                evidence=evidence,
            )
        return results

    def _match_niche_signature(
        self,
        niche_fractions: dict[str, float],
        cluster_id: str,
    ) -> tuple[str, float, dict]:
        """Match neighbourhood composition against known niche signatures."""
        tier2_hint = self.tier2_labels.get(cluster_id, "")
        candidates = [tier2_hint] + [c for c in NICHE_SIGNATURES if c != tier2_hint] \
                     if tier2_hint else list(NICHE_SIGNATURES.keys())

        best_label, best_score = "", 0.0
        best_evidence: dict = {}

        for ct in candidates:
            conditions = NICHE_SIGNATURES.get(ct)
            if conditions is None:
                continue
            n_conditions = len(conditions)
            n_met = 0
            condition_evidence: dict = {}

            for keyword, min_fraction in conditions:
                # Sum fractions of neighbour clusters whose label contains the keyword
                matched_fraction = sum(
                    frac for label, frac in niche_fractions.items()
                    if keyword.lower() in label.lower()
                )
                met = matched_fraction >= max(min_fraction, NICHE_FRACTION_THRESH)
                condition_evidence[keyword] = {
                    "fraction": round(matched_fraction, 3),
                    "min_required": min_fraction,
                    "met": met,
                }
                if met:
                    n_met += 1

            score = n_met / max(n_conditions, 1)
            if score > best_score:
                best_score = score
                best_label = ct
                best_evidence = {
                    "mode": "niche_composition",
                    "conditions": condition_evidence,
                    "top_neighbours": dict(sorted(niche_fractions.items(), key=lambda x: -x[1])[:5]),
                }

        confidence = min(1.0, best_score)
        if confidence < 0.25:
            # No niche match: use tier2 hint if available, else unknown
            return (tier2_hint or "unknown"), 0.2, best_evidence
        return best_label, confidence, best_evidence

    def _unavailable_results(
        self,
        adata: "anndata.AnnData",
        cluster_key: str,
    ) -> dict[str, StreamResult]:
        if cluster_key not in adata.obs.columns:
            return {}
        return {
            str(cid): StreamResult(
                cluster_id=str(cid), stream=self.STREAM_NAME,
                label="", confidence=0.0, available=False,
                evidence={"reason": "no_spatial_data"},
            )
            for cid in adata.obs[cluster_key].unique()
        }
