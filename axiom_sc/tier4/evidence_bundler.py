"""
Tier 4 EvidenceBundler: formats multi-tier evidence into LLM-ready text.

CRITICAL DIFFERENCE FROM CASSIA:
  CASSIA sends marker genes to the LLM (single evidence stream).
  AXIOM-SC Tier 4 sends the FULL evidence bundle (Tier 1-3 results +
  velocity direction + chromatin state + spatial context + AXIOM rules fired).
  LLMs reason over all evidence simultaneously, not just marker genes.

License: Apache 2.0
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from axiom_sc.tier2.evidence import EvidenceBundle
from axiom_sc.tier4.prompts import build_user_prompt, SYSTEM_PROMPT


@dataclass
class BundledEvidence:
    """
    LLM-ready evidence package for one cluster.

    Produced by EvidenceBundler.bundle() — passed to LLMEnsemble.annotate().
    """

    cluster_id: str
    system_prompt: str
    user_prompt: str
    raw_fields: dict = field(default_factory=dict)  # structured fields for logging


class EvidenceBundler:
    """
    Assembles structured evidence from Tiers 1-3 into a prompt-ready text block.

    Parameters
    ----------
    kg_retriever : optional KGRetriever
        If provided, injects RAG context from the KG.
    """

    def __init__(self, kg_retriever=None):
        self._rag = kg_retriever

    def bundle(
        self,
        evidence: EvidenceBundle,
        tier1_predictions: dict[str, float] | None = None,
        tier2_rules_fired: list[dict] | None = None,
        tier2_summary: str = "",
        tier3_stream_results: list | None = None,
    ) -> BundledEvidence:
        """
        Build the LLM evidence package for one cluster.

        Parameters
        ----------
        evidence : EvidenceBundle
            Tier 2 evidence (marker genes + regulons + spatial context).
        tier1_predictions : dict[label -> confidence], optional
        tier2_rules_fired : list of rule firing dicts, optional
        tier2_summary : str
            Human-readable Tier 2 outcome (e.g. "1 PROVEN, 3 UNCERTAIN, 2 CONTRADICTED").
        tier3_stream_results : list[StreamResult], optional
        """
        t1_preds = tier1_predictions or {}
        t2_rules = tier2_rules_fired or []
        t3_results = tier3_stream_results or []

        # Tier 3 parsing
        t3_streams_available = [r.stream for r in t3_results if r.available]
        t3_agreement = self._format_t3_agreement(t3_results)
        t3_details = self._format_t3_details(t3_results)

        # Spatial context from EvidenceBundle or Tier 3 spatial niche
        spatial_ctx = self._extract_spatial_context(evidence, t3_results)

        # RAG context
        rag_ctx = ""
        if self._rag is not None:
            rag_result = self._rag.retrieve_for_cluster(
                cluster_id=evidence.cluster_id,
                marker_genes=evidence.marker_genes,
                top_k=5,
            )
            rag_ctx = rag_result.as_prompt_text()

        user_prompt = build_user_prompt(
            cluster_id=evidence.cluster_id,
            tier1_predictions=t1_preds,
            tier2_rules_fired=t2_rules,
            tier2_summary=tier2_summary,
            tier3_streams_available=t3_streams_available,
            tier3_agreement=t3_agreement,
            tier3_details=t3_details,
            marker_genes=evidence.marker_genes,
            regulons=evidence.regulons,
            spatial_context=spatial_ctx,
            rag_context=rag_ctx,
        )

        raw_fields = {
            "cluster_id": evidence.cluster_id,
            "n_markers": len(evidence.marker_genes),
            "n_regulons": len(evidence.regulons),
            "n_tier2_rules": len(t2_rules),
            "n_tier3_streams": len(t3_streams_available),
            "has_rag": self._rag is not None,
        }

        return BundledEvidence(
            cluster_id=evidence.cluster_id,
            system_prompt=SYSTEM_PROMPT,
            user_prompt=user_prompt,
            raw_fields=raw_fields,
        )

    # ── private helpers ───────────────────────────────────────────────────────

    def _format_t3_agreement(self, t3_results: list) -> str:
        available = [r for r in t3_results if r.available and r.label]
        if not available:
            return "no streams available"
        labels = [r.label for r in available]
        majority = max(set(labels), key=labels.count)
        n_agree = labels.count(majority)
        return f"{n_agree}/{len(available)} streams → {majority!r}"

    def _format_t3_details(self, t3_results: list) -> str:
        lines = []
        for r in t3_results:
            status = f"conf={r.confidence:.2f}" if r.available else "unavailable"
            label = r.label or "(no label)"
            lines.append(f"  {r.stream}: {label} [{status}]")
        return "\n".join(lines) if lines else ""

    def _extract_spatial_context(self, evidence: EvidenceBundle, t3_results: list) -> str:
        if evidence.spatial_context:
            ctx = evidence.spatial_context
            if isinstance(ctx, dict):
                parts = []
                if "niche_label" in ctx:
                    parts.append(f"niche={ctx['niche_label']}")
                if "tissue" in ctx:
                    parts.append(f"tissue={ctx['tissue']}")
                return "  " + ", ".join(parts) if parts else "  (spatial dict — no label)"
            return f"  {ctx}"
        for r in t3_results:
            if r.stream == "spatial_niche" and r.available and r.label:
                return f"  spatial_niche stream: {r.label} (conf={r.confidence:.2f})"
        return ""
