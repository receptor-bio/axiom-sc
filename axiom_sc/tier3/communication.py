"""
Ligand-receptor communication stream for Tier 3 multi-stream convergence.

Two operational modes (in priority order):
  1. COMMOT mode: uses pre-computed communication scores from
     adata.obsm['commot-sum-sender-*'] / adata.obsm['commot-sum-receiver-*'].
  2. Expression mode: infers communication capacity directly from ligand and
     receptor gene z-scores in adata. Used when COMMOT has not been run.

Key discriminations (from CLAUDE.md Section 10):
  - CAF vs SMC: CAFs are high TGFβ senders; SMCs are PDGF receivers
  - Monocyte vs tissue-resident Mφ: monocytes send CCL2; Mφ send SPP1/CXCL10
  - pDC: strong IFN-I sender upon stimulation (IFNA1, IFNB1 high)
  - Treg: TGFβ sender + IL-2 receiver
  - GCB: CXCL13 receiver (from FDCs in germinal center)

Reference:
  Cang et al. (2023) Nature Communications 14:7706
  doi: 10.1038/s41467-023-43600-9

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

# ── External L-R database (Patch B) ──────────────────────────────────────────
_LR_INDEX_PATH = Path.home() / ".axiom_sc" / "resources" / "communication" / "lr_gene_index.json"


def _load_lr_database() -> dict:
    """Load CellPhoneDB L-R gene index when pre-downloaded. Falls back to {}."""
    if _LR_INDEX_PATH.exists():
        try:
            with open(_LR_INDEX_PATH) as f:
                db = json.load(f)
            logger.debug("Loaded external L-R database: %d genes", len(db))
            return db
        except Exception as exc:
            logger.warning("Could not load L-R database from %s: %s", _LR_INDEX_PATH, exc)
    return {}


LR_DATABASE: dict = _load_lr_database()

# Ligand/receptor gene sets per cell-type communication role.
# Format: {cell_type: {"sender_genes": [...], "receiver_genes": [...], "threshold": float}}
# threshold = minimum fraction of sender/receiver genes that must be active (z ≥ 1.5)
# for the profile to match.
LR_SIGNATURES: dict[str, dict] = {
    "CAF": {
        "sender_genes": ["TGFB1", "TGFB2", "FGF7", "FGF10"],
        "receiver_genes": [],
        "threshold": 0.4,
        "notes": "CAFs are primary TGFβ/FGF senders in stromal niche",
    },
    "SMC": {
        "sender_genes": ["NOTCH3"],
        "receiver_genes": ["PDGFRB", "ACVRL1"],
        "threshold": 0.4,
        "notes": "SMCs receive PDGF; express Notch3 for vessel communication",
    },
    "Treg": {
        "sender_genes": ["TGFB1", "IL10"],
        "receiver_genes": ["IL2RA", "TNFRSF18"],
        "threshold": 0.4,
        "notes": "Tregs send TGFβ/IL-10 and receive IL-2 survival signals",
    },
    "pDC": {
        "sender_genes": ["IFNA1", "IFNA2", "IFNB1"],
        "receiver_genes": [],
        "threshold": 0.33,
        "notes": "pDCs are primary IFN-I senders upon viral/TLR stimulation",
    },
    "NK": {
        "sender_genes": ["GZMB", "PRF1", "IFNG"],
        "receiver_genes": ["IL15RA", "IL12RB1"],
        "threshold": 0.33,
        "notes": "NK cells send cytotoxic effectors; receive IL-15/IL-12 survival",
    },
    "Tissue_Resident_Mac": {
        "sender_genes": ["SPP1", "CXCL10", "CXCL9"],
        "receiver_genes": ["CSF1R", "IL34"],
        "threshold": 0.33,
        "notes": "Tissue-resident Mφ receive CSF1; send SPP1 (osteopontin)",
    },
    "Monocyte": {
        "sender_genes": ["CCL2", "CCL7", "S100A8"],
        "receiver_genes": [],
        "threshold": 0.33,
        "notes": "Monocytes send CCL2 to recruit myeloid cells",
    },
    "GCB": {
        "sender_genes": ["LTA", "TNFSF13B"],
        "receiver_genes": ["CXCR4", "CXCR5"],
        "threshold": 0.4,
        "notes": "GCBs receive CXCL13 from FDCs via CXCR5",
    },
    "ILC2": {
        "sender_genes": ["IL13", "IL5", "AREG"],
        "receiver_genes": ["IL25R", "IL33R"],
        "threshold": 0.33,
        "notes": "ILC2 send IL-13/IL-5; receive alarmins IL-25/IL-33 from epithelium",
    },
    "Endothelial": {
        "sender_genes": ["CXCL12", "CXCL13", "BMP4"],
        "receiver_genes": ["KDR", "FLT1", "NOTCH1"],
        "threshold": 0.33,
        "notes": "Endothelium sends CXCL12; receives VEGF via KDR",
    },
}

# COMMOT obsm key patterns to detect pre-computed communication scores
COMMOT_SENDER_PREFIX = "commot-sum-sender"
COMMOT_RECEIVER_PREFIX = "commot-sum-receiver"

# Z-score threshold for "active" ligand/receptor
GENE_ACTIVE_THRESH: float = 1.5


def _try_import_commot():
    try:
        import commot as ct
        return ct
    except ImportError:
        return None


class CommunicationStream:
    """
    Ligand-receptor signaling stream.

    Annotates clusters based on their communication role (sender vs receiver)
    for known L-R pathways. Uses pre-computed COMMOT scores when available;
    falls back to marker gene expression inference otherwise.

    Parameters
    ----------
    gene_activities_key : str
        Key in adata.obsm containing precomputed COMMOT sender scores.
        Auto-detected if not specified (searches for commot-sum-sender-* keys).
    tier2_labels : dict[str, str] | None
        Tier 2 label hints per cluster_id.
    z_active_thresh : float
        Z-score threshold for calling a ligand/receptor gene "active".
    """

    STREAM_NAME = "communication"

    def __init__(
        self,
        tier2_labels: dict[str, str] | None = None,
        z_active_thresh: float = GENE_ACTIVE_THRESH,
    ):
        self.tier2_labels = tier2_labels or {}
        self.z_active_thresh = z_active_thresh

    def run(
        self,
        adata: "anndata.AnnData",
        cluster_key: str = "leiden",
    ) -> dict[str, StreamResult]:
        """
        Run communication stream on all clusters.

        Parameters
        ----------
        adata : AnnData
            May contain adata.obsm['commot-sum-sender-*'] (pre-computed).
            Falls back to using adata.X expression for known ligand/receptor genes.
        cluster_key : str

        Returns
        -------
        dict[str, StreamResult]
        """
        if cluster_key not in adata.obs.columns:
            logger.warning("cluster_key %r not in adata.obs — communication stream skipped", cluster_key)
            return {}

        # Mode 1: COMMOT pre-computed scores
        commot_data = self._load_commot_scores(adata)
        if commot_data is not None:
            logger.debug("Communication stream: using pre-computed COMMOT scores")
            return self._score_from_commot(adata, cluster_key, commot_data)

        # Mode 2: expression-based L-R inference
        lr_expr = self._load_lr_expression(adata)
        if lr_expr is not None:
            logger.debug("Communication stream: using ligand/receptor gene expression")
            return self._score_from_expression(adata, cluster_key, lr_expr)

        logger.info(
            "Communication stream: no COMMOT data and no L-R genes in adata. "
            "Returning unavailable."
        )
        return self._unavailable_results(adata, cluster_key)

    # ── COMMOT mode ───────────────────────────────────────────────────────────

    def _load_commot_scores(
        self,
        adata: "anndata.AnnData",
    ) -> "dict[str, pd.DataFrame] | None":
        """
        Detect and load COMMOT sender/receiver score matrices from adata.obsm.

        Returns dict with "sender" and/or "receiver" DataFrames, or None.
        """
        try:
            import pandas as pd
        except ImportError:
            return None

        sender_keys = [k for k in adata.obsm.keys() if k.startswith(COMMOT_SENDER_PREFIX)]
        receiver_keys = [k for k in adata.obsm.keys() if k.startswith(COMMOT_RECEIVER_PREFIX)]

        if not sender_keys and not receiver_keys:
            return None

        result: dict[str, "pd.DataFrame"] = {}
        for key in sender_keys:
            raw = adata.obsm[key]
            df = self._obsm_to_df(raw, adata, key)
            if df is not None:
                result["sender"] = df
                break  # use first found

        for key in receiver_keys:
            raw = adata.obsm[key]
            df = self._obsm_to_df(raw, adata, key)
            if df is not None:
                result["receiver"] = df
                break

        return result if result else None

    def _obsm_to_df(self, raw, adata, key: str) -> "pd.DataFrame | None":
        """Convert an obsm entry to a labeled DataFrame."""
        try:
            import pandas as pd
        except ImportError:
            return None

        if isinstance(raw, "pd.DataFrame" if False else object) and hasattr(raw, "columns"):
            return raw

        arr = np.asarray(raw)
        if arr.ndim == 1:
            arr = arr.reshape(-1, 1)

        # Try to get pathway names from adata.uns
        uns_key = key.replace("commot-sum-sender-", "commot_sender_pathway_names_") \
                     .replace("commot-sum-receiver-", "commot_receiver_pathway_names_")
        col_names = adata.uns.get(uns_key)
        if col_names is None:
            col_names = [f"pathway_{i}" for i in range(arr.shape[1])]

        return pd.DataFrame(arr, index=adata.obs_names, columns=list(col_names)[:arr.shape[1]])

    def _score_from_commot(
        self,
        adata: "anndata.AnnData",
        cluster_key: str,
        commot_data: "dict[str, pd.DataFrame]",
    ) -> dict[str, StreamResult]:
        """Score clusters using pre-computed COMMOT pathway scores."""
        clusters = adata.obs[cluster_key]
        unique_clusters = clusters.unique()
        results: dict[str, StreamResult] = {}

        sender_df = commot_data.get("sender")
        receiver_df = commot_data.get("receiver")

        for cid in unique_clusters:
            mask = (clusters == cid).values
            evidence: dict = {}

            if sender_df is not None:
                cluster_sender = sender_df.loc[mask].mean(axis=0)
                evidence["sender_pathways"] = {
                    col: round(float(v), 3)
                    for col, v in cluster_sender.nlargest(5).items()
                }
            if receiver_df is not None:
                cluster_receiver = receiver_df.loc[mask].mean(axis=0)
                evidence["receiver_pathways"] = {
                    col: round(float(v), 3)
                    for col, v in cluster_receiver.nlargest(5).items()
                }

            label, confidence = self._match_commot_to_signature(evidence, str(cid))
            results[str(cid)] = StreamResult(
                cluster_id=str(cid),
                stream=self.STREAM_NAME,
                label=label,
                confidence=confidence,
                available=True,
                evidence=evidence,
            )
        return results

    def _match_commot_to_signature(
        self,
        evidence: dict,
        cluster_id: str,
    ) -> tuple[str, float]:
        """Match COMMOT pathway evidence to known cell type communication profiles."""
        tier2_hint = self.tier2_labels.get(cluster_id, "")
        top_sender = set(evidence.get("sender_pathways", {}).keys())
        top_receiver = set(evidence.get("receiver_pathways", {}).keys())

        best_label, best_score = "", 0.0
        candidates = [tier2_hint] + [c for c in LR_SIGNATURES if c != tier2_hint] if tier2_hint else list(LR_SIGNATURES.keys())

        for ct, sig in LR_SIGNATURES.items():
            sender_overlap = len(set(sig["sender_genes"]) & {s.upper() for s in top_sender})
            receiver_overlap = len(set(sig["receiver_genes"]) & {r.upper() for r in top_receiver})
            n_expected = len(sig["sender_genes"]) + len(sig["receiver_genes"])
            if n_expected == 0:
                continue
            score = (sender_overlap + receiver_overlap) / n_expected
            if score > best_score:
                best_score = score
                best_label = ct

        confidence = min(1.0, best_score * 2.0)  # scale: 0.5 match → confidence 1.0
        if confidence < 0.2:
            return "unknown", confidence
        return best_label, confidence

    # ── expression mode ───────────────────────────────────────────────────────

    def _load_lr_expression(
        self,
        adata: "anndata.AnnData",
    ) -> "pd.DataFrame | None":
        """
        Extract expression of known L-R genes from adata.

        Uses the full CellPhoneDB L-R index (2,300+ genes) when pre-downloaded;
        falls back to the hardcoded LR_SIGNATURES dict (~40 genes).
        Returns cells × genes DataFrame for genes present in adata.var_names.
        """
        try:
            import pandas as pd
        except ImportError:
            return None

        all_lr_genes: set[str] = set()
        if LR_DATABASE:
            # Full CellPhoneDB index (Patch B)
            all_lr_genes.update(LR_DATABASE.keys())
        # Always include hardcoded signature genes (needed for scoring logic)
        for sig in LR_SIGNATURES.values():
            all_lr_genes.update(sig["sender_genes"])
            all_lr_genes.update(sig["receiver_genes"])

        present_genes = [g for g in all_lr_genes if g in adata.var_names]
        if not present_genes:
            return None

        gene_idx = [adata.var_names.get_loc(g) for g in present_genes]
        X = adata.X
        if hasattr(X, "toarray"):
            X = X.toarray()
        else:
            X = np.asarray(X)

        return pd.DataFrame(X[:, gene_idx], index=adata.obs_names, columns=present_genes)

    def _score_from_expression(
        self,
        adata: "anndata.AnnData",
        cluster_key: str,
        lr_expr: "pd.DataFrame",
    ) -> dict[str, StreamResult]:
        """Score clusters using L-R gene expression z-scores."""
        try:
            import pandas as pd
        except ImportError:
            return self._unavailable_results(adata, cluster_key)

        # Z-score each gene
        mean = lr_expr.mean(axis=0)
        std = lr_expr.std(axis=0).replace(0, 1)
        lr_zscores = (lr_expr - mean) / std

        clusters = adata.obs[cluster_key]
        unique_clusters = clusters.unique()
        results: dict[str, StreamResult] = {}

        for cid in unique_clusters:
            mask = (clusters == cid).values
            cluster_mean = lr_zscores.loc[mask].mean(axis=0)

            label, confidence, evidence = self._match_lr_expression(cluster_mean, str(cid))
            results[str(cid)] = StreamResult(
                cluster_id=str(cid),
                stream=self.STREAM_NAME,
                label=label,
                confidence=confidence,
                available=True,
                evidence=evidence,
            )
        return results

    def _match_lr_expression(
        self,
        cluster_mean: "pd.Series",
        cluster_id: str,
    ) -> tuple[str, float, dict]:
        """Match L-R expression profile to known signatures."""
        tier2_hint = self.tier2_labels.get(cluster_id, "")
        candidates = [tier2_hint] + [c for c in LR_SIGNATURES if c != tier2_hint] if tier2_hint else list(LR_SIGNATURES.keys())

        best_label, best_score = "", -1.0
        best_evidence: dict = {}

        for ct in candidates:
            sig = LR_SIGNATURES.get(ct)
            if sig is None:
                continue
            score, evidence = self._score_lr_signature(cluster_mean, sig, ct)
            if score > best_score:
                best_score = score
                best_label = ct
                best_evidence = evidence

        confidence = min(1.0, max(0.0, best_score))
        if confidence < 0.25:
            return "unknown", confidence, best_evidence
        return best_label, confidence, best_evidence

    def _score_lr_signature(
        self,
        cluster_mean: "pd.Series",
        sig: dict,
        ct: str,
    ) -> tuple[float, dict]:
        """Score a cluster's L-R expression against one cell-type signature."""
        n_tested, n_active = 0, 0
        evidence: dict = {"active_senders": [], "active_receivers": [], "missing": []}

        for gene in sig["sender_genes"] + sig["receiver_genes"]:
            if gene not in cluster_mean.index:
                evidence["missing"].append(gene)
                continue
            n_tested += 1
            z = float(cluster_mean[gene])
            is_active = z >= self.z_active_thresh
            if is_active:
                n_active += 1
                if gene in sig["sender_genes"]:
                    evidence["active_senders"].append(gene)
                else:
                    evidence["active_receivers"].append(gene)

        if n_tested == 0:
            return 0.0, evidence
        threshold = sig.get("threshold", 0.4)
        fraction = n_active / n_tested
        # Score: 1.0 if at or above threshold; scaled below
        score = min(1.0, fraction / max(threshold, 0.1))
        return score, evidence

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
                evidence={"reason": "no_lr_data"},
            )
            for cid in adata.obs[cluster_key].unique()
        }
