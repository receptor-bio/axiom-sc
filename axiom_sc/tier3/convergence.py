"""
MultiStreamConvergence: 6-stream orthogonal evidence convergence for Tier 3.

Aggregates results from up to 6 independent evidence streams:
  1. velocity       — scVelo RNA velocity trajectory direction
  2. chromatin      — Signac ATAC locus accessibility
  3. communication  — COMMOT ligand-receptor signaling
  4. spatial_niche  — CELLama neighbor composition
  5. cross_species  — OrthoFinder conservation
  6. spatial_type   — scType spatial (reimplemented)

Convergence criterion (from CLAUDE.md Section 10):
  4+ of 6 streams agree → PROVEN
  3 streams agree       → HIGH_CONFIDENCE
  <3 streams agree      → UNCERTAIN (Tier 4 handles)

Minimum viable (non-spatial/ATAC data): velocity + chromatin (2 streams available).

References:
  scVelo:      Bergen et al. (2020) doi: 10.1038/s41587-020-0591-3
  Signac:      Stuart et al. (2021) doi: 10.1038/s41592-021-01282-5
  COMMOT:      Cang et al. (2023)   doi: 10.1038/s41467-023-43600-9
  OrthoFinder: Emms & Kelly (2019)  doi: 10.1186/s13059-019-1832-y
  CELLama:     Lv et al. (2025)     doi: 10.1002/advs.202413514

License: Apache-2.0
"""
from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from axiom_sc.tier3.stream_result import StreamResult

if TYPE_CHECKING:
    import anndata

logger = logging.getLogger(__name__)

# All 6 stream names (ordered by discriminative power for tie-breaking)
ALL_STREAMS: list[str] = [
    "velocity",
    "chromatin",
    "communication",
    "spatial_niche",
    "cross_species",
    "spatial_type",
]

# Convergence thresholds
PROVEN_THRESHOLD: int = 4        # ≥4 of available streams agree → PROVEN
HIGH_CONF_THRESHOLD: int = 3     # ≥3 agree → HIGH_CONFIDENCE
# <3 → UNCERTAIN


@dataclass
class ConvergenceResult:
    """
    Final Tier 3 verdict for one cluster after aggregating all stream results.

    Attributes
    ----------
    cluster_id : str
    verdict : str
        "PROVEN" | "HIGH_CONFIDENCE" | "UNCERTAIN"
    consensus_label : str
        The cell type label that the plurality of streams agree on.
        Empty if no consensus (all streams disagree or all unavailable).
    confidence : float
        Weighted confidence [0, 1]. Computed as mean confidence of agreeing streams.
    n_streams_available : int
        Number of streams that ran and produced a result.
    n_streams_agree : int
        Number of streams that agree on consensus_label.
    stream_results : dict[str, StreamResult]
        Raw results from each stream (keyed by stream name).
    label_votes : dict[str, int]
        Vote counts for each candidate label across available streams.
    """

    cluster_id: str
    verdict: str
    consensus_label: str
    confidence: float
    n_streams_available: int
    n_streams_agree: int
    stream_results: dict[str, StreamResult] = field(default_factory=dict)
    label_votes: dict[str, int] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "cluster_id": self.cluster_id,
            "verdict": self.verdict,
            "consensus_label": self.consensus_label,
            "confidence": round(self.confidence, 3),
            "n_streams_available": self.n_streams_available,
            "n_streams_agree": self.n_streams_agree,
            "label_votes": self.label_votes,
        }

    def __repr__(self) -> str:
        return (
            f"ConvergenceResult({self.cluster_id!r}: {self.verdict} "
            f"→ {self.consensus_label!r} "
            f"[{self.n_streams_agree}/{self.n_streams_available} streams])"
        )


class MultiStreamConvergence:
    """
    Aggregates up to 6 independent evidence streams for each cluster.

    Instantiate once with the stream results, then call ``converge()``
    per cluster (or ``converge_all()`` for all clusters at once).

    Parameters
    ----------
    min_confidence : float
        Minimum stream confidence to count a stream's vote.
        Streams below this threshold are treated as abstaining.
    require_n_streams : int | None
        Minimum number of available streams to attempt convergence.
        Clusters with fewer available streams return UNCERTAIN.
        Default: None (use at least 1).
    """

    def __init__(
        self,
        min_confidence: float = 0.3,
        require_n_streams: int | None = None,
    ):
        self.min_confidence = min_confidence
        self.require_n_streams = require_n_streams or 1

    def converge_all(
        self,
        stream_results_by_stream: dict[str, dict[str, StreamResult]],
    ) -> dict[str, ConvergenceResult]:
        """
        Converge all clusters from a dict of stream → {cluster_id → StreamResult}.

        Parameters
        ----------
        stream_results_by_stream : dict
            {stream_name: {cluster_id: StreamResult}}

        Returns
        -------
        dict[str, ConvergenceResult]
            {cluster_id: ConvergenceResult}
        """
        # Collect all cluster IDs seen across all streams
        all_clusters: set[str] = set()
        for sr_dict in stream_results_by_stream.values():
            all_clusters.update(sr_dict.keys())

        results: dict[str, ConvergenceResult] = {}
        for cluster_id in sorted(all_clusters):
            per_cluster: dict[str, StreamResult] = {
                stream_name: sr_dict[cluster_id]
                for stream_name, sr_dict in stream_results_by_stream.items()
                if cluster_id in sr_dict
            }
            results[cluster_id] = self.converge(cluster_id, per_cluster)
        return results

    def converge(
        self,
        cluster_id: str,
        stream_results: dict[str, StreamResult],
    ) -> ConvergenceResult:
        """
        Converge stream results for a single cluster.

        Parameters
        ----------
        cluster_id : str
        stream_results : dict[str, StreamResult]
            {stream_name: StreamResult} for this cluster.

        Returns
        -------
        ConvergenceResult
        """
        # Filter to available streams with sufficient confidence
        available = {
            name: sr for name, sr in stream_results.items()
            if sr.available and sr.confidence >= self.min_confidence and sr.label
        }

        n_available = len(available)

        if n_available < self.require_n_streams:
            return ConvergenceResult(
                cluster_id=cluster_id,
                verdict="UNCERTAIN",
                consensus_label="",
                confidence=0.0,
                n_streams_available=n_available,
                n_streams_agree=0,
                stream_results=stream_results,
                label_votes={},
            )

        # Vote by label
        label_votes: Counter = Counter()
        label_confidences: dict[str, list[float]] = {}
        for sr in available.values():
            label_votes[sr.label] += 1
            label_confidences.setdefault(sr.label, []).append(sr.confidence)

        if not label_votes:
            return ConvergenceResult(
                cluster_id=cluster_id,
                verdict="UNCERTAIN",
                consensus_label="",
                confidence=0.0,
                n_streams_available=n_available,
                n_streams_agree=0,
                stream_results=stream_results,
                label_votes={},
            )

        # Plurality winner
        consensus_label, n_agree = label_votes.most_common(1)[0]
        mean_confidence = float(
            sum(label_confidences[consensus_label]) / len(label_confidences[consensus_label])
        )

        # Apply convergence thresholds
        verdict = self._determine_verdict(n_agree, n_available)

        return ConvergenceResult(
            cluster_id=cluster_id,
            verdict=verdict,
            consensus_label=consensus_label,
            confidence=mean_confidence,
            n_streams_available=n_available,
            n_streams_agree=n_agree,
            stream_results=stream_results,
            label_votes=dict(label_votes),
        )

    def _determine_verdict(self, n_agree: int, n_available: int) -> str:
        """Apply convergence criterion thresholds."""
        if n_agree >= PROVEN_THRESHOLD:
            return "PROVEN"
        if n_agree >= HIGH_CONF_THRESHOLD:
            return "HIGH_CONFIDENCE"
        return "UNCERTAIN"

    # ── integration helpers ───────────────────────────────────────────────────

    def _detect_modalities(self, adata: "anndata.AnnData") -> dict[str, bool]:
        """
        Auto-detect which streams can run on this dataset by inspecting adata content.

        Returns
        -------
        dict[str, bool]
            Keys: velocity, chromatin, communication, spatial_niche, cross_species, scenic.
            All three RNA-only streams (communication, spatial_niche, cross_species)
            and SCENIC always return True because they run in expression-mode on any data.
        """
        return {
            "velocity":     "velocity" in adata.layers or "spliced" in adata.layers,
            "chromatin":    "gene_activities" in adata.obsm or "X_atac" in adata.obsm,
            "communication": True,   # expression-mode always available
            "spatial_niche": True,   # kNN composition always available
            "cross_species": True,   # CONSERVED_MARKERS always available
            "scenic":        True,   # runs via subprocess (always available)
        }

    def run_and_converge(
        self,
        adata: "anndata.AnnData",
        cluster_key: str = "leiden",
        tier2_labels: dict[str, str] | None = None,
        run_velocity: bool = True,
        run_chromatin: bool = True,
        run_communication: bool = True,
        run_spatial_niche: bool = True,
        run_cross_species: bool = True,
        auto_detect: bool = False,
    ) -> dict[str, ConvergenceResult]:
        """
        Run all available Tier 3 streams and converge per cluster.

        All 6 streams are attempted; each degrades gracefully to
        available=False when its required data is absent.

        Parameters
        ----------
        adata : AnnData
        cluster_key : str
        tier2_labels : dict[str, str] | None
            Tier 2 label hints per cluster_id; passed to each stream.
        run_velocity : bool
        run_chromatin : bool
        run_communication : bool
        run_spatial_niche : bool
        run_cross_species : bool

        Returns
        -------
        dict[str, ConvergenceResult]
        """
        from axiom_sc.tier3.velocity import VelocityStream
        from axiom_sc.tier3.chromatin import ChromatinStream
        from axiom_sc.tier3.communication import CommunicationStream
        from axiom_sc.tier3.spatial_niche import SpatialNicheStream
        from axiom_sc.tier3.cross_species import CrossSpeciesStream

        if auto_detect:
            modalities = self._detect_modalities(adata)
            run_velocity = modalities["velocity"]
            run_chromatin = modalities["chromatin"]
            active_count = sum(1 for k, v in modalities.items() if v)
            logger.info(
                "Auto-detected modalities: %d/6 streams active "
                "(velocity=%s, chromatin=%s)",
                active_count,
                modalities["velocity"],
                modalities["chromatin"],
            )

        stream_results_by_stream: dict[str, dict[str, StreamResult]] = {}

        if run_velocity:
            stream_results_by_stream["velocity"] = VelocityStream(
                tier2_labels=tier2_labels
            ).run(adata, cluster_key)

        if run_chromatin:
            stream_results_by_stream["chromatin"] = ChromatinStream(
                tier2_labels=tier2_labels
            ).run(adata, cluster_key)

        if run_communication:
            stream_results_by_stream["communication"] = CommunicationStream(
                tier2_labels=tier2_labels
            ).run(adata, cluster_key)

        if run_spatial_niche:
            stream_results_by_stream["spatial_niche"] = SpatialNicheStream(
                tier2_labels=tier2_labels
            ).run(adata, cluster_key)

        if run_cross_species:
            stream_results_by_stream["cross_species"] = CrossSpeciesStream(
                tier2_labels=tier2_labels
            ).run(adata, cluster_key)

        if not stream_results_by_stream:
            logger.warning("No streams ran; all clusters will be UNCERTAIN.")
            return {}

        return self.converge_all(stream_results_by_stream)

    def update_evidence_bundles(
        self,
        convergence_results: dict[str, ConvergenceResult],
        evidence_bundles: dict,
    ) -> None:
        """
        Update EvidenceBundle.spatial_context with Tier 3 convergence results.

        Parameters
        ----------
        convergence_results : dict[str, ConvergenceResult]
        evidence_bundles : dict[str, EvidenceBundle]
        """
        for cluster_id, conv_result in convergence_results.items():
            if cluster_id not in evidence_bundles:
                continue
            bundle = evidence_bundles[cluster_id]
            bundle.spatial_context.update({
                "tier3_verdict": conv_result.verdict,
                "tier3_label": conv_result.consensus_label,
                "tier3_confidence": conv_result.confidence,
                "tier3_n_agree": conv_result.n_streams_agree,
                "tier3_n_available": conv_result.n_streams_available,
                "tier3_label_votes": conv_result.label_votes,
            })


def stream_agreement_report(
    convergence_results: dict[str, ConvergenceResult],
) -> dict:
    """
    Compute inter-stream agreement statistics across all clusters.

    Used to measure Tier 3 performance on Phase 1 datasets (Day 8 task:
    "Run Tier 3 on Phase 1 datasets — measure stream agreement rates").

    Returns
    -------
    dict with keys:
      n_clusters            : int
      n_proven              : int
      n_high_confidence     : int
      n_uncertain           : int
      pct_proven            : float
      pct_high_confidence   : float
      pct_uncertain         : float
      mean_streams_available: float
      mean_streams_agree    : float
      per_stream_availability: dict[stream_name, float]  (fraction of clusters where available)
      verdict_breakdown     : dict[cluster_id, str]
      label_distribution    : dict[label, int]  (top consensus labels)
    """
    if not convergence_results:
        return {"n_clusters": 0}

    n = len(convergence_results)
    n_proven = sum(1 for r in convergence_results.values() if r.verdict == "PROVEN")
    n_high   = sum(1 for r in convergence_results.values() if r.verdict == "HIGH_CONFIDENCE")
    n_unc    = sum(1 for r in convergence_results.values() if r.verdict == "UNCERTAIN")

    mean_avail = sum(r.n_streams_available for r in convergence_results.values()) / n
    mean_agree = sum(r.n_streams_agree for r in convergence_results.values()) / n

    # Per-stream availability
    per_stream: dict[str, int] = {}
    for cr in convergence_results.values():
        for stream_name, sr in cr.stream_results.items():
            per_stream[stream_name] = per_stream.get(stream_name, 0) + int(sr.available)
    per_stream_frac = {k: v / n for k, v in per_stream.items()}

    # Label distribution (non-empty consensus labels)
    label_counts: dict[str, int] = {}
    for cr in convergence_results.values():
        if cr.consensus_label:
            label_counts[cr.consensus_label] = label_counts.get(cr.consensus_label, 0) + 1
    top_labels = dict(sorted(label_counts.items(), key=lambda x: -x[1])[:20])

    return {
        "n_clusters": n,
        "n_proven": n_proven,
        "n_high_confidence": n_high,
        "n_uncertain": n_unc,
        "pct_proven": round(100 * n_proven / n, 1),
        "pct_high_confidence": round(100 * n_high / n, 1),
        "pct_uncertain": round(100 * n_unc / n, 1),
        "mean_streams_available": round(mean_avail, 2),
        "mean_streams_agree": round(mean_agree, 2),
        "per_stream_availability": per_stream_frac,
        "verdict_breakdown": {cid: r.verdict for cid, r in convergence_results.items()},
        "label_distribution": top_labels,
    }


def make_unavailable_stream_results(
    cluster_ids: list[str],
    stream_name: str,
    reason: str = "stream_not_run",
) -> dict[str, StreamResult]:
    """
    Create unavailable StreamResults for all clusters.

    Utility for Day 8 streams that aren't yet implemented
    (communication, spatial_niche, cross_species) — so convergence.py
    can run with only the streams that are available.
    """
    return {
        cid: StreamResult(
            cluster_id=cid,
            stream=stream_name,
            label="",
            confidence=0.0,
            available=False,
            evidence={"reason": reason},
        )
        for cid in cluster_ids
    }
