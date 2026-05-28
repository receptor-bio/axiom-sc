"""
Tier 4 RAG (Retrieval-Augmented Generation) over the AXIOM KG.

Retrieves the most relevant KG cell types given a cluster's marker gene profile.
Uses marker-gene overlap scoring (no vector DB required — KG is the knowledge base).

The RAG context is injected into the LLM prompt so the model reasons over the
structured biological knowledge already in the KG, not just its training priors.

License: Apache 2.0
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class RAGResult:
    """Top-K cell types retrieved by marker overlap for one cluster."""

    cluster_id: str
    top_types: list[dict] = field(default_factory=list)  # [{cell_type, score, key_markers}]

    def as_prompt_text(self) -> str:
        """Format for injection into LLM prompt."""
        if not self.top_types:
            return "  (no KG entries matched cluster markers)"
        lines = []
        for entry in self.top_types:
            ct = entry["cell_type"]
            score = entry["overlap_score"]
            markers = ", ".join(entry["key_markers"][:6])
            lines.append(f"  {ct} (overlap={score:.2f}): key markers — {markers}")
        return "\n".join(lines)


class KGRetriever:
    """
    Retrieves relevant cell types from the KG given a cluster's marker profile.

    Scoring: Jaccard overlap between query markers and KG rule genes,
    weighted by rule confidence (high=1.0, medium=0.7, low=0.4).
    Returns top-K entries sorted by weighted overlap score.
    """

    CONFIDENCE_WEIGHT = {"high": 1.0, "medium": 0.7, "low": 0.4}

    def __init__(self, kg_path: str | Path):
        self._index: dict[str, dict] = {}  # cell_type → {genes, weighted_genes}
        self._load(Path(kg_path))

    def _load(self, kg_path: Path) -> None:
        with open(kg_path) as f:
            rules = json.load(f)
        for rule in rules:
            if rule.get("status") != "ACTIVE":
                continue
            if rule.get("rule_type") == "negative":
                continue  # negatives exclude types — don't boost retrieval
            ct = rule["cell_type"]
            if ct not in self._index:
                self._index[ct] = {"genes": set(), "weighted": {}}
            weight = self.CONFIDENCE_WEIGHT.get(rule.get("confidence", "low"), 0.4)
            for gene in rule.get("gene_or_regulon", []):
                self._index[ct]["genes"].add(gene)
                self._index[ct]["weighted"][gene] = max(
                    self._index[ct]["weighted"].get(gene, 0.0), weight
                )

    def retrieve(
        self,
        marker_genes: dict[str, float],
        top_k: int = 5,
        min_overlap: float = 0.0,
    ) -> list[dict]:
        """
        Return top-K cell types sorted by weighted marker overlap.

        Parameters
        ----------
        marker_genes : dict[gene -> z_score]
        top_k : int
        min_overlap : float  minimum score to include in results
        """
        query_genes = set(marker_genes.keys())
        scores = []
        for ct, idx in self._index.items():
            kg_genes = idx["genes"]
            matched = query_genes & kg_genes
            if not matched:
                continue
            weighted_score = sum(idx["weighted"].get(g, 0.0) for g in matched)
            normalized = weighted_score / max(len(kg_genes), 1)
            if normalized < min_overlap:
                continue
            # Key markers = intersection sorted by query z-score descending
            key_markers = sorted(matched, key=lambda g: -marker_genes.get(g, 0.0))
            scores.append({
                "cell_type": ct,
                "overlap_score": round(normalized, 4),
                "key_markers": key_markers,
                "n_kg_genes": len(kg_genes),
                "n_matched": len(matched),
            })
        scores.sort(key=lambda x: -x["overlap_score"])
        return scores[:top_k]

    def retrieve_for_cluster(
        self,
        cluster_id: str,
        marker_genes: dict[str, float],
        top_k: int = 5,
    ) -> RAGResult:
        """Convenience wrapper returning a RAGResult."""
        top = self.retrieve(marker_genes, top_k=top_k)
        return RAGResult(cluster_id=cluster_id, top_types=top)

    @property
    def n_cell_types(self) -> int:
        return len(self._index)
