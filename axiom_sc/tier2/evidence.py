"""
EvidenceBundle: structured container for per-cluster evidence passed to Tier 2.

Holds marker gene z-scores and regulon AUC z-scores (post-SCENIC normalisation).
All z-scores are dataset-level z-scores (Phase 1 validated: always z-score per
regulon per dataset before any threshold comparison).

License: Apache 2.0
"""
from __future__ import annotations

from dataclasses import dataclass, field

# Phase 1 validated thresholds (hardcoded, do not change without re-validation)
Z_ACTIVE_THRESH: float = 1.5    # z ≥ 1.5 = gene/regulon active / high / present
Z_INACTIVE_THRESH: float = -0.5  # z ≤ -0.5 = gene/regulon inactive / low / absent


@dataclass
class EvidenceBundle:
    """
    All evidence available for one cluster at Tier 2 decision time.

    marker_genes : dict[gene_symbol -> z_score]
        Cluster-level marker gene z-scores (e.g. from scanpy rank_genes_groups).
        Positive z = enriched in cluster vs rest.

    regulons : dict[regulon_name -> z_score]
        Per-regulon AUC values, z-scored across all cells in the dataset
        (Phase 1 lesson: raw AUC values not comparable; must z-score per run).
        Regulon name = TF symbol (e.g. "FOXP3", "IRF7").

    cluster_id : str
        Identifier for this cluster (used in logging and output).

    tissue : str
        Tissue of origin — used to filter tissue-specific KG rules.

    spatial_context : dict
        Optional spatial metadata from Tier 3 (added by convergence.py).
    """

    cluster_id: str
    marker_genes: dict[str, float] = field(default_factory=dict)
    regulons: dict[str, float] = field(default_factory=dict)
    tissue: str = ""
    spatial_context: dict = field(default_factory=dict)

    # ── lookup helpers ────────────────────────────────────────────────────────

    def get_marker(self, gene: str) -> float | None:
        """Returns z-score from marker_genes, or None if not in evidence."""
        return self.marker_genes.get(gene)

    def get_regulon(self, regulon: str) -> float | None:
        """Returns z-score from regulons, or None if not in evidence."""
        return self.regulons.get(regulon)

    def lookup(self, gene: str, evidence_source: str) -> float | None:
        """
        Look up a gene value from the appropriate evidence source.
        For 'regulon': searches regulons dict.
        For all other sources (marker_genes, chromatin, etc.): searches marker_genes.
        """
        if evidence_source == "regulon":
            return self.regulons.get(gene)
        return self.marker_genes.get(gene)

    def lookup_any(self, gene: str) -> float | None:
        """
        Check regulons first, then marker_genes.
        Used by circuit rules where components may span both evidence types.
        Returns None only if gene is in neither dict.
        """
        val = self.regulons.get(gene)
        if val is not None:
            return val
        return self.marker_genes.get(gene)

    def satisfies_direction(self, value: float, direction: str) -> bool:
        """
        Test whether a numeric z-score satisfies the rule direction threshold.

        Directions:
          high / present / active  → z ≥ Z_ACTIVE_THRESH  (1.5)
          low / absent / inactive  → z ≤ Z_INACTIVE_THRESH (-0.5)
        """
        if direction in ("high", "present", "active"):
            return value >= Z_ACTIVE_THRESH
        if direction in ("low", "absent", "inactive"):
            return value <= Z_INACTIVE_THRESH
        raise ValueError(f"Unknown direction: {direction!r}")

    def __repr__(self) -> str:
        return (
            f"EvidenceBundle(cluster={self.cluster_id!r}, "
            f"markers={len(self.marker_genes)}, regulons={len(self.regulons)})"
        )
