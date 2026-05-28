"""
Chromatin accessibility stream for Tier 3 multi-stream convergence.

Uses ATAC-seq gene activity scores (from Signac via rpy2, or pre-computed) to
verify cell type identity from chromatin state. Key loci checked:
  - FOXP3 locus accessible → Treg identity (fixes Phase 1 SCENIC miss)
  - AIRE locus accessible → mature mTEC identity (fixes Phase 1 SCENIC miss)
  - PSMB11 locus → cortical TEC identity
  - TBX21, GATA3, RORC → effector T / ILC identity

Reference:
  Stuart et al. (2021) Nature Methods 18:1272-1279
  doi: 10.1038/s41592-021-01282-5

License: Apache-2.0
"""
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np

from axiom_sc.tier3.stream_result import StreamResult

if TYPE_CHECKING:
    import anndata
    import pandas as pd

logger = logging.getLogger(__name__)

# Z-score threshold for "accessible" locus
ACCESSIBILITY_ACTIVE_THRESH: float = 1.0   # looser than regulon (1.5) — ATAC signal noisier
ACCESSIBILITY_INACTIVE_THRESH: float = -0.5

# Locus-to-cell-type mapping.
# Each key is a cell type label; value is a list of (gene_locus, required_direction) tuples.
# "high" = locus must be accessible; "low" = locus must be inaccessible.
# All conditions in a list are ANDed: every condition must be satisfied for PROVEN.
LOCUS_SIGNATURES: dict[str, list[tuple[str, str]]] = {
    "Treg": [
        ("FOXP3", "high"),   # FOXP3 locus must be open (Phase 1 SCENIC miss fixed here)
        ("IL2RA", "high"),   # CD25 locus
    ],
    "Mature_mTEC": [
        ("AIRE", "high"),    # AIRE locus open → post-Aire mTEC (Phase 1 SCENIC miss)
        ("PSMB11", "low"),   # PSMB11 closes as mTECs mature past cortical stage
    ],
    "Cortical_TEC": [
        ("PSMB11", "high"),  # Cortical TEC marker
    ],
    "Th1": [
        ("TBX21", "high"),
        ("GATA3", "low"),
        ("RORC", "low"),
    ],
    "Th2": [
        ("GATA3", "high"),
        ("TBX21", "low"),
        ("RORC", "low"),
    ],
    "Th17": [
        ("RORC", "high"),
        ("GATA3", "low"),
        ("TBX21", "low"),
    ],
    "pDC": [
        ("IRF7", "high"),
        ("SIGLEC1", "high"),
    ],
    "cDC1": [
        ("IRF8", "high"),
        ("BATF3", "high"),
    ],
    "ILC1": [
        ("TBX21", "high"),
        ("GATA3", "low"),
    ],
    "ILC2": [
        ("GATA3", "high"),
        ("RORC", "low"),
    ],
    "ILC3": [
        ("RORC", "high"),
        ("NCR2", "high"),
    ],
    "Naive_T": [
        ("TCF7", "high"),
        ("PDCD1", "low"),   # PD-1 locus closed in naïve T
    ],
    "Exhausted_T": [
        ("TOX", "high"),
        ("PDCD1", "high"),  # PD-1 locus open
        ("TCF7", "low"),
    ],
    "Progenitor_Exhausted_T": [
        ("TCF7", "high"),   # TCF7 locus still open
        ("TOX", "high"),    # TOX beginning to open
        ("PDCD1", "high"),
    ],
}

# Minimal set of loci that discriminate key ambiguous pairs
PRIORITY_LOCI: list[str] = [
    "FOXP3", "AIRE", "PSMB11", "TBX21", "GATA3", "RORC",
    "IRF7", "IRF8", "BATF3", "NCR2", "TOX", "TCF7", "PDCD1",
    "IL2RA", "SIGLEC1",
]


def _try_import_rpy2():
    """Returns (rpy2.robjects, True) or (None, False) if rpy2 not installed."""
    try:
        import rpy2.robjects as ro
        return ro, True
    except (ImportError, RuntimeError):
        return None, False


class ChromatinStream:
    """
    Chromatin accessibility stream using Signac gene activity scores.

    Two operational modes:
      1. Pre-computed: uses gene activity scores already stored in adata.obsm
         under key gene_activities_key (e.g., adata.obsm['gene_activities']).
      2. rpy2/Signac: calls R Signac::GeneActivity() via rpy2 (optional).
         Falls back to mode 1 if rpy2 or R is unavailable.

    When neither mode is available (no ATAC data at all), returns unavailable
    StreamResults for all clusters so convergence.py can skip this stream.

    Parameters
    ----------
    gene_activities_key : str
        Key in adata.obsm containing the gene activity score matrix.
        Columns should be gene symbols.
    use_rpy2 : bool
        If True, attempt Signac gene activity computation via rpy2.
        Ignored if rpy2 / R / Signac are not installed.
    tier2_labels : dict[str, str] | None
        Tier 2 verdicts per cluster. Used to prioritise loci to check.
    """

    STREAM_NAME = "chromatin"

    def __init__(
        self,
        gene_activities_key: str = "gene_activities",
        use_rpy2: bool = True,
        tier2_labels: dict[str, str] | None = None,
    ):
        self.gene_activities_key = gene_activities_key
        self.use_rpy2 = use_rpy2
        self.tier2_labels = tier2_labels or {}

    def run(
        self,
        adata: "anndata.AnnData",
        cluster_key: str = "leiden",
    ) -> dict[str, StreamResult]:
        """
        Run chromatin accessibility stream on all clusters.

        Parameters
        ----------
        adata : AnnData
            Should contain ATAC gene activity scores. May contain:
              - adata.obsm[gene_activities_key]: pre-computed gene activity matrix
              - adata.layers['atac_counts']: raw peak matrix (for Signac compute)
        cluster_key : str
            Column in adata.obs with cluster labels.

        Returns
        -------
        dict[str, StreamResult]
            Per-cluster StreamResult. available=False when no ATAC data found.
        """
        if cluster_key not in adata.obs.columns:
            logger.warning("cluster_key %r not in adata.obs — chromatin stream skipped", cluster_key)
            return self._unavailable_results(adata, cluster_key)

        gene_activities = self._load_gene_activities(adata)
        if gene_activities is None:
            logger.info("No gene activity scores available; chromatin stream unavailable.")
            return self._unavailable_results(adata, cluster_key)

        return self._score_clusters(adata, cluster_key, gene_activities)

    # ── gene activity loading ─────────────────────────────────────────────────

    def _load_gene_activities(
        self,
        adata: "anndata.AnnData",
    ) -> "pd.DataFrame | None":
        """
        Load or compute gene activity scores as a cells × genes DataFrame.

        Priority:
          1. adata.obsm[gene_activities_key] if present (pre-computed)
          2. Signac::GeneActivity() via rpy2 if requested
          3. None (stream unavailable)
        """
        # 1. Pre-computed scores
        if self.gene_activities_key in adata.obsm:
            ga = adata.obsm[self.gene_activities_key]
            return self._to_dataframe(ga, adata)

        # 2. Signac via rpy2
        if self.use_rpy2:
            ga = self._run_signac_gene_activity(adata)
            if ga is not None:
                return ga

        return None

    def _to_dataframe(
        self,
        gene_activities,
        adata: "anndata.AnnData",
    ) -> "pd.DataFrame | None":
        """Convert gene activity matrix to cells × genes DataFrame."""
        try:
            import pandas as pd
        except ImportError:
            logger.warning("pandas not installed; chromatin stream unavailable.")
            return None

        if hasattr(gene_activities, "toarray"):
            # scipy sparse
            arr = gene_activities.toarray()
        else:
            arr = np.asarray(gene_activities)

        # Try to get gene names from adata.uns or infer from matrix shape
        gene_names = adata.uns.get("gene_activity_genes")
        if gene_names is None and arr.shape[1] == len(adata.var_names):
            gene_names = list(adata.var_names)
        if gene_names is None:
            # Can't map genes — stream unavailable
            logger.warning(
                "gene_activities matrix has no associated gene names. "
                "Set adata.uns['gene_activity_genes'] or ensure columns match adata.var_names."
            )
            return None

        return pd.DataFrame(arr, index=adata.obs_names, columns=gene_names)

    def _run_signac_gene_activity(
        self,
        adata: "anndata.AnnData",
    ) -> "pd.DataFrame | None":
        """
        Compute gene activity matrix via R Signac::GeneActivity().
        Requires rpy2, R ≥ 4.0, Signac ≥ 1.9, BSgenome.Hsapiens.UCSC.hg38.

        Falls back gracefully when any of these are missing.
        """
        ro, has_rpy2 = _try_import_rpy2()
        if not has_rpy2:
            logger.debug("rpy2 not available; skipping Signac gene activity computation.")
            return None

        if "atac_counts" not in adata.layers and "peak_counts" not in adata.layers:
            logger.debug("No ATAC peak count layer found; skipping Signac computation.")
            return None

        try:
            import pandas as pd
            return self._signac_gene_activity_rpy2(adata, ro, pd)
        except Exception as exc:
            logger.warning("Signac gene activity computation failed: %s", exc)
            return None

    def _signac_gene_activity_rpy2(
        self,
        adata: "anndata.AnnData",
        ro,
        pd,
    ) -> "pd.DataFrame | None":
        """Internal: call Signac::GeneActivity via rpy2 and return as DataFrame."""
        import tempfile, os

        peak_layer = "atac_counts" if "atac_counts" in adata.layers else "peak_counts"
        peak_matrix = adata.layers[peak_layer]
        if hasattr(peak_matrix, "toarray"):
            peak_matrix = peak_matrix.toarray()

        # Write peak matrix and peak names to tempfiles for R
        with tempfile.TemporaryDirectory() as tmpdir:
            peaks_rds = os.path.join(tmpdir, "peaks.rds")
            barcodes_file = os.path.join(tmpdir, "barcodes.txt")
            genes_out = os.path.join(tmpdir, "gene_activity.csv")

            # Write barcodes
            np.savetxt(barcodes_file, adata.obs_names.to_numpy(), fmt="%s")

            # Load Signac and run GeneActivity via R
            ro.r(f"""
            library(Signac)
            library(Seurat)
            library(BSgenome.Hsapiens.UCSC.hg38)

            peak_matrix <- readRDS("{peaks_rds}")
            chrom_assay <- CreateChromatinAssay(
              counts = peak_matrix,
              sep = c(":", "-"),
              genome = BSgenome.Hsapiens.UCSC.hg38
            )
            obj <- CreateSeuratObject(counts = chrom_assay, assay = "peaks")
            gene_activity <- GeneActivity(obj)
            write.csv(as.matrix(gene_activity), file = "{genes_out}")
            """)

            if not os.path.exists(genes_out):
                return None

            ga_df = pd.read_csv(genes_out, index_col=0).T  # cells × genes
            return ga_df

    # ── scoring ───────────────────────────────────────────────────────────────

    def _score_clusters(
        self,
        adata: "anndata.AnnData",
        cluster_key: str,
        gene_activities: "pd.DataFrame",
    ) -> dict[str, StreamResult]:
        """Score each cluster against locus signatures."""
        try:
            import pandas as pd
        except ImportError:
            return self._unavailable_results(adata, cluster_key)

        clusters = adata.obs[cluster_key]
        unique_clusters = clusters.unique()

        # Z-score gene activity per locus across all cells
        gene_zscores = self._zscore_activities(gene_activities)

        results: dict[str, StreamResult] = {}
        for cid in unique_clusters:
            mask = (clusters == cid).values
            cluster_mean = gene_zscores.loc[mask].mean(axis=0)
            label, confidence, evidence = self._match_locus_signature(
                cluster_mean, str(cid)
            )
            results[str(cid)] = StreamResult(
                cluster_id=str(cid),
                stream=self.STREAM_NAME,
                label=label,
                confidence=confidence,
                available=True,
                evidence=evidence,
            )
        return results

    def _zscore_activities(self, gene_activities: "pd.DataFrame") -> "pd.DataFrame":
        """Z-score each gene's activity across all cells."""
        import pandas as pd

        mean = gene_activities.mean(axis=0)
        std = gene_activities.std(axis=0).replace(0, 1)  # avoid div/0 for constant genes
        return (gene_activities - mean) / std

    def _match_locus_signature(
        self,
        cluster_mean: "pd.Series",
        cluster_id: str,
    ) -> tuple[str, float, dict]:
        """
        Score cluster against all locus signatures and return best match.

        Returns (label, confidence, evidence_dict).
        """
        tier2_hint = self.tier2_labels.get(cluster_id, "")

        # Check tier2 hint first if available
        candidates = list(LOCUS_SIGNATURES.keys())
        if tier2_hint and tier2_hint in LOCUS_SIGNATURES:
            candidates = [tier2_hint] + [c for c in candidates if c != tier2_hint]

        best_label = ""
        best_score = -1.0
        best_evidence: dict = {}

        for cell_type, conditions in LOCUS_SIGNATURES.items():
            score, evidence = self._score_signature(cluster_mean, conditions, cell_type)
            if score > best_score:
                best_score = score
                best_label = cell_type
                best_evidence = evidence

        # Convert raw score to confidence [0, 1]
        # Perfect score = 1.0 (all conditions met), 0 = no conditions met
        confidence = min(1.0, max(0.0, best_score))

        # If no signature reached threshold, return unknown
        if confidence < 0.3:
            return "unknown", confidence, best_evidence

        return best_label, confidence, best_evidence

    def _score_signature(
        self,
        cluster_mean: "pd.Series",
        conditions: list[tuple[str, str]],
        cell_type: str,
    ) -> tuple[float, dict]:
        """
        Score how well a cluster matches a locus signature.

        Returns (score 0-1, evidence_dict).
        Score = fraction of testable conditions satisfied.
        """
        n_testable = 0
        n_satisfied = 0
        evidence: dict = {"loci_checked": {}, "loci_missing": []}

        for gene, direction in conditions:
            if gene not in cluster_mean.index:
                evidence["loci_missing"].append(gene)
                continue
            z = float(cluster_mean[gene])
            n_testable += 1
            if direction == "high":
                satisfied = z >= ACCESSIBILITY_ACTIVE_THRESH
            else:  # "low"
                satisfied = z <= ACCESSIBILITY_INACTIVE_THRESH
            n_satisfied += int(satisfied)
            evidence["loci_checked"][gene] = {
                "z": round(z, 2),
                "direction": direction,
                "satisfied": satisfied,
            }

        if n_testable == 0:
            return 0.0, evidence
        return n_satisfied / n_testable, evidence

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
                evidence={"reason": "no_atac_data"},
            )
            for cid in clusters
        }
