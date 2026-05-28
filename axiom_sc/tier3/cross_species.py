"""
Cross-species conservation stream for Tier 3 multi-stream convergence.

Validates cell type annotation by checking whether a cluster's key markers
are conserved across species. Strongly conserved marker programs support
fundamental, well-annotated cell identities.

Two operational modes:
  1. OrthoFinder mode: uses pre-computed OrthoFinder results stored in
     adata.uns['cross_species_conservation'] (researcher-run offline).
  2. Conserved marker mode: scores cluster markers against a curated set of
     known cross-species conserved markers per cell type (from published
     cross-species single-cell atlases).

Key principle: a cell type whose canonical TF/marker program is conserved
across mouse-human (and ideally other vertebrates) is a fundamental, ancient
cell type — our annotation confidence is higher.

Reference:
  Emms & Kelly (2019) Genome Biology 20:238
  doi: 10.1186/s13059-019-1832-y

License: Apache-2.0
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from axiom_sc.tier3.stream_result import StreamResult

if TYPE_CHECKING:
    import anndata
    import pandas as pd

logger = logging.getLogger(__name__)

# ── External conservation index (Patch B) ────────────────────────────────────
_CONSERVATION_PATH = Path.home() / ".axiom_sc" / "resources" / "orthology" / "conservation_index.json"


def _load_conservation_index() -> dict:
    """Load ENSEMBL one-to-one ortholog index when pre-downloaded. Falls back to {}."""
    if _CONSERVATION_PATH.exists():
        try:
            with open(_CONSERVATION_PATH) as f:
                index = json.load(f)
            logger.debug("Loaded conservation index: %d genes", len(index))
            return index
        except Exception as exc:
            logger.warning("Could not load conservation index from %s: %s", _CONSERVATION_PATH, exc)
    return {}


CONSERVATION_INDEX: dict = _load_conservation_index()

# Per cell type: genes that are one-to-one orthologs across human/mouse/rat
# and whose expression in the cell type is conserved (from cross-species atlases).
# Sources: Tabula Sapiens, Tabula Muris, cross-species thymus/lung/brain atlases.
CONSERVED_MARKERS: dict[str, set[str]] = {
    # T cell subtypes
    "Treg": {"FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18", "TIGIT"},
    "Th1": {"TBX21", "IFNG", "CXCR3", "IL12RB2", "HAVCR2"},
    "Th2": {"GATA3", "IL4", "IL5", "IL13", "PTGDR2", "CRTH2"},
    "Th17": {"RORC", "IL17A", "IL17F", "CCR6", "IL23R"},
    "Naive_T": {"CCR7", "TCF7", "SELL", "LEF1", "CD62L"},
    "Memory_T": {"IL7R", "TCF7", "CXCR3"},
    "Exhausted_T": {"TOX", "PDCD1", "LAG3", "HAVCR2", "TIGIT", "ENTPD1"},
    "Progenitor_Exhausted_T": {"TCF7", "TOX", "PDCD1", "SLAMF6"},
    "CD8_Effector": {"GZMB", "PRF1", "NKG7", "GZMH", "GNLY"},
    # NK / ILC
    "NK": {"NCAM1", "KLRD1", "NKG7", "GNLY", "GZMB", "PRF1"},
    "ILC1": {"TBX21", "IFNG", "NCR1"},
    "ILC2": {"GATA3", "IL13", "IL5", "KLRG1"},
    "ILC3": {"RORC", "NCR2", "NCR3", "IL22"},
    # B cell / plasma
    "Naive_B": {"MS4A1", "CD19", "IGHD", "TCL1A", "FCER2"},
    "Memory_B": {"MS4A1", "CD19", "CD27", "AIM2"},
    "GCB": {"BCL6", "AICDA", "CXCR4", "CXCR5", "RGS13"},
    "Plasma_cell": {"JCHAIN", "PRDM1", "XBP1", "SDC1", "IGHG1"},
    "Pre_B": {"VPREB1", "IGLL1", "RAG1", "RAG2"},
    "Pro_B": {"CD34", "MME", "DNTT", "RAG1"},
    # Myeloid
    "Monocyte": {"CD14", "LYZ", "S100A8", "S100A9", "CSF1R"},
    "Classical_Monocyte": {"CD14", "LYZ", "S100A8", "FCGR3A"},
    "Non_Classical_Monocyte": {"FCGR3A", "CX3CR1", "CDKB"},
    "pDC": {"IRF7", "SIGLEC1", "LILRA4", "GZMB", "CLEC4C"},
    "cDC1": {"CLEC9A", "XCR1", "CADM1", "BATF3", "IRF8"},
    "cDC2": {"CD1C", "FCER1A", "CLEC10A", "FCGR2A"},
    "Tissue_Resident_Mac": {"CSF1R", "CD68", "MRC1", "CD163", "LYVE1"},
    "M1_Macrophage": {"IDO1", "NOS2", "CXCL9", "CXCL10", "CD86"},
    "M2_Macrophage": {"MRC1", "CD163", "TREM2", "PPARG"},
    "Mast_cell": {"TPSAB1", "KIT", "FCER1A", "CPA3"},
    "Basophil": {"BASOPH", "FUT7"},
    "Eosinophil": {"EPX", "RNASE2", "SMPD3"},
    "Neutrophil": {"S100A8", "CXCR2", "FCGR3B", "CSF3R"},
    # DC subsets
    "pDC_precursor": {"IRF8", "LTB", "ITM2C"},
    # Thymic epithelial
    "Mature_mTEC": {"AIRE", "FN1", "CCL19", "CCL21"},
    "Cortical_TEC": {"PSMB11", "PRSS16", "CTSV"},
    # Stromal
    "CAF": {"COL1A1", "FAP", "ACTA2", "PDPN", "POSTN"},
    "SMC": {"ACTA2", "MYH11", "TAGLN", "CNN1"},
    "Fibroblast": {"COL1A1", "COL1A2", "PDGFRA", "DCN", "LUM"},
    "Endothelial": {"PECAM1", "CDH5", "CLDN5", "FLT1", "ENG"},
    # Lung epithelial
    "AT1": {"AGER", "PDPN", "CLDN18", "AQP5"},
    "AT2": {"SFTPB", "SFTPC", "ABCA3", "NKX2-1"},
    "Club": {"SCGB1A1", "SCGB3A2", "CYP2F1"},
    "Ciliated": {"FOXJ1", "DNAI1", "TUBA1A", "PIFO"},
}

# Minimum fraction of conserved markers that must be detected (above threshold)
# in the cluster for a match to be called
CONSERVATION_THRESH: float = 0.3

# Z-score threshold for calling a gene "detected" in a cluster
GENE_DETECTED_THRESH: float = 1.0


class CrossSpeciesStream:
    """
    Cross-species conservation stream.

    Scores each cluster's expressed genes against curated conserved marker sets.
    High conservation fraction → annotation is a well-established, ancient cell type.

    Parameters
    ----------
    orthofinder_uns_key : str
        Key in adata.uns for pre-computed OrthoFinder conservation scores.
        If found, used directly; otherwise falls back to conserved marker scoring.
    marker_genes_obs_key : str | None
        Column in adata.obs pointing to pre-computed per-cluster marker gene lists.
        If None, cluster markers are derived from adata.X directly.
    tier2_labels : dict[str, str] | None
        Tier 2 label hints per cluster_id.
    """

    STREAM_NAME = "cross_species"

    def __init__(
        self,
        orthofinder_uns_key: str = "cross_species_conservation",
        tier2_labels: dict[str, str] | None = None,
    ):
        self.orthofinder_uns_key = orthofinder_uns_key
        self.tier2_labels = tier2_labels or {}

    def run(
        self,
        adata: "anndata.AnnData",
        cluster_key: str = "leiden",
    ) -> dict[str, StreamResult]:
        """
        Run cross-species conservation stream on all clusters.

        Parameters
        ----------
        adata : AnnData
        cluster_key : str

        Returns
        -------
        dict[str, StreamResult]
        """
        if cluster_key not in adata.obs.columns:
            logger.warning("cluster_key %r not in adata.obs — cross_species stream skipped", cluster_key)
            return {}

        # Mode 1: pre-computed OrthoFinder results
        if self.orthofinder_uns_key in adata.uns:
            return self._score_from_orthofinder(adata, cluster_key)

        # Mode 2: conserved marker gene expression
        return self._score_from_conserved_markers(adata, cluster_key)

    # ── OrthoFinder mode ──────────────────────────────────────────────────────

    def _score_from_orthofinder(
        self,
        adata: "anndata.AnnData",
        cluster_key: str,
    ) -> dict[str, StreamResult]:
        """
        Score clusters using pre-computed OrthoFinder conservation table.

        Expected format of adata.uns[orthofinder_uns_key]:
          {cluster_id: {"conservation_score": float, "label": str}}
        OR:
          {cluster_id: float}  (plain conservation score)
        """
        conservation_data = adata.uns[self.orthofinder_uns_key]
        clusters = adata.obs[cluster_key].unique()
        results: dict[str, StreamResult] = {}

        for cid in clusters:
            cid_str = str(cid)
            entry = conservation_data.get(cid_str, conservation_data.get(cid))
            if entry is None:
                results[cid_str] = StreamResult(
                    cluster_id=cid_str, stream=self.STREAM_NAME,
                    label="", confidence=0.0, available=True,
                    evidence={"reason": "not_in_orthofinder_table"},
                )
                continue

            if isinstance(entry, dict):
                score = float(entry.get("conservation_score", 0.0))
                label = entry.get("label", self.tier2_labels.get(cid_str, ""))
            else:
                score = float(entry)
                label = self.tier2_labels.get(cid_str, "")

            results[cid_str] = StreamResult(
                cluster_id=cid_str,
                stream=self.STREAM_NAME,
                label=label,
                confidence=min(1.0, score),
                available=True,
                evidence={"mode": "orthofinder", "conservation_score": score},
            )
        return results

    # ── conserved marker mode ─────────────────────────────────────────────────

    def _score_from_conserved_markers(
        self,
        adata: "anndata.AnnData",
        cluster_key: str,
    ) -> dict[str, StreamResult]:
        """
        Score clusters by measuring overlap with curated conserved marker sets.

        For each cluster, computes mean expression z-score of each gene,
        then checks what fraction of each cell type's conserved markers
        are detected (z ≥ GENE_DETECTED_THRESH) in the cluster.
        """
        try:
            import pandas as pd
        except ImportError:
            return self._unavailable_results(adata, cluster_key)

        # Build set of all conserved genes that are in adata.
        # When the ENSEMBL index is loaded (Patch B), every gene in adata can be
        # checked for conservation — not just the ~300 in CONSERVED_MARKERS.
        all_conserved = set()
        if CONSERVATION_INDEX:
            # Use genome-wide conservation index: include all conserved genes present in adata
            all_conserved.update(g for g in adata.var_names if g in CONSERVATION_INDEX)
        for gene_set in CONSERVED_MARKERS.values():
            all_conserved.update(gene_set)
        present_genes = [g for g in all_conserved if g in adata.var_names]

        if not present_genes:
            logger.info(
                "No conserved marker genes found in adata.var_names. "
                "Returning tier2 hints or unknown labels."
            )
            return self._fallback_tier2_results(adata, cluster_key)

        # Extract expression for conserved genes
        gene_idx = [list(adata.var_names).index(g) for g in present_genes]
        X = adata.X
        if hasattr(X, "toarray"):
            X = X.toarray()
        else:
            X = np.asarray(X)
        expr = X[:, gene_idx]

        # Z-score across all cells
        mean = expr.mean(axis=0)
        std = expr.std(axis=0)
        std[std == 0] = 1.0
        zscored = (expr - mean) / std

        gene_to_col = {g: i for i, g in enumerate(present_genes)}

        clusters = adata.obs[cluster_key]
        unique_clusters = clusters.unique()
        results: dict[str, StreamResult] = {}

        for cid in unique_clusters:
            cid_str = str(cid)
            mask = (clusters == cid).values
            cluster_mean_z = zscored[mask].mean(axis=0)  # (n_present_genes,)

            label, confidence, evidence = self._match_conserved_signature(
                cluster_mean_z, gene_to_col, cid_str
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

    def _match_conserved_signature(
        self,
        cluster_mean_z: np.ndarray,
        gene_to_col: dict[str, int],
        cluster_id: str,
    ) -> tuple[str, float, dict]:
        """Find the best conserved marker match for a cluster."""
        tier2_hint = self.tier2_labels.get(cluster_id, "")
        candidates = [tier2_hint] + [c for c in CONSERVED_MARKERS if c != tier2_hint] \
                     if tier2_hint else list(CONSERVED_MARKERS.keys())

        best_label, best_score = "", 0.0
        best_evidence: dict = {}

        for ct in candidates:
            gene_set = CONSERVED_MARKERS.get(ct)
            if gene_set is None:
                continue

            present_in_adata = [g for g in gene_set if g in gene_to_col]
            if not present_in_adata:
                continue

            n_detected = sum(
                1 for g in present_in_adata
                if float(cluster_mean_z[gene_to_col[g]]) >= GENE_DETECTED_THRESH
            )
            fraction = n_detected / len(gene_set)

            if fraction > best_score:
                best_score = fraction
                best_label = ct
                best_evidence = {
                    "mode": "conserved_markers",
                    "cell_type": ct,
                    "n_detected": n_detected,
                    "n_present": len(present_in_adata),
                    "n_in_signature": len(gene_set),
                    "fraction_detected": round(fraction, 3),
                    "detected_genes": [
                        g for g in present_in_adata
                        if float(cluster_mean_z[gene_to_col[g]]) >= GENE_DETECTED_THRESH
                    ],
                }

        confidence = min(1.0, best_score / max(CONSERVATION_THRESH, 0.1))
        confidence = min(1.0, confidence)
        if confidence < 0.2:
            return (tier2_hint or "unknown"), 0.1, best_evidence
        return best_label, confidence, best_evidence

    def _fallback_tier2_results(
        self,
        adata: "anndata.AnnData",
        cluster_key: str,
    ) -> dict[str, StreamResult]:
        """Return tier2 hint labels when no genes present for scoring."""
        return {
            str(cid): StreamResult(
                cluster_id=str(cid), stream=self.STREAM_NAME,
                label=self.tier2_labels.get(str(cid), ""),
                confidence=0.3 if self.tier2_labels.get(str(cid)) else 0.0,
                available=bool(self.tier2_labels.get(str(cid))),
                evidence={"reason": "no_conserved_genes_in_adata"},
            )
            for cid in adata.obs[cluster_key].unique()
        }

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
                evidence={"reason": "no_cross_species_data"},
            )
            for cid in adata.obs[cluster_key].unique()
        }
