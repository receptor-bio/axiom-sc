"""
Tier 4 LLM prompt templates for cell type annotation.

All prompts are versioned and deterministic: same evidence bundle → same prompt text.
Prompt versioning enables reproducibility audits — different PROMPT_VERSION strings
indicate the LLM saw different instructions, allowing comparison of annotation runs.

License: Apache 2.0
"""
from __future__ import annotations

PROMPT_VERSION = "4.0"

SYSTEM_PROMPT = """\
You are an expert computational biologist annotating clusters from a single-cell
RNA-seq experiment. You receive multi-omic evidence assembled by the AXIOM-SC pipeline
from upstream tiers. Your task: determine the most likely cell type, or declare a
potential novel cell state if no known type fits the evidence.

Key rules:
1. A cell type claim requires positive evidence (marker expression or regulon activity).
2. If Tier 2 KG has CONTRADICTED a cell type, do NOT propose it under any circumstances.
3. AXIOM-SC uses proof-by-contradiction: one hard mechanistic violation disqualifies a type.
4. If you are uncertain between two cell types, say so explicitly with mechanistic reasoning.
5. Respond ONLY in the required JSON format — no prose outside the JSON block.
"""

USER_TEMPLATE = """\
## Cluster: {cluster_id}

### Tier 1 Evidence (AXIOMTier1 MLP Ensemble)
Top predictions (label: confidence):
{tier1_predictions}

### Tier 2 Evidence (AXIOM Knowledge Graph — proof-by-contradiction)
Rules fired:
{tier2_rules_fired}
Tier 2 summary: {tier2_summary}

### Tier 3 Evidence (Multi-stream Orthogonal Convergence)
Streams run: {tier3_streams_available}
Stream agreement: {tier3_agreement}
{tier3_details}

### Top Marker Genes (cluster vs rest, sorted by z-score descending)
{marker_genes}

### Active Regulons (SCENIC AUC z-score ≥ 1.5)
{active_regulons}

### Spatial / Niche Context
{spatial_context}

### RAG — Closest KG cell types by marker overlap
{rag_context}

Based on ALL the above evidence, identify the cell type.
Return ONLY a JSON object with these exact keys:
{{
  "cell_type": "<best cell type label or 'NOVEL_CANDIDATE'>",
  "confidence": <float 0.0-1.0>,
  "reasoning": "<mechanistic justification referencing specific evidence above>",
  "alternatives": ["<second choice>", "<third choice>"],
  "novel_candidate": <true|false>
}}
"""

NOVEL_FOLLOWUP_TEMPLATE = """\
The cluster {cluster_id} was declared a NOVEL_CANDIDATE by the LLM ensemble.
Based on the active circuits and closest known attractors, generate a mechanistic
characterization:

Active TF circuits: {active_circuits}
Closest known velocity attractor: {velocity_sink}
Spatial niche description: {spatial_context}
Cross-species conservation (top hits): {cross_species}

Return ONLY a JSON object:
{{
  "nearest_ontology_term": "<closest Cell Ontology term>",
  "proposed_name": "<descriptive working name for this candidate state>",
  "perturbation_predictions": [
    {{"gene": "<TF/gene>", "perturbation": "KO|OE", "predicted_effect": "<description>"}}
  ],
  "confidence_in_novelty": <float 0.0-1.0>
}}
"""


def format_tier1_predictions(predictions: dict[str, float]) -> str:
    """Format Tier 1 top predictions as indented bullet list."""
    if not predictions:
        return "  (no Tier 1 predictions available)"
    items = sorted(predictions.items(), key=lambda x: -x[1])[:5]
    return "\n".join(f"  {label}: {conf:.3f}" for label, conf in items)


def format_tier2_rules(rules_fired: list[dict]) -> str:
    """Format Tier 2 rule firings as a compact table."""
    if not rules_fired:
        return "  (Tier 2 not run or no rules fired)"
    lines = []
    for rf in rules_fired:
        verdict = rf.get("verdict", "?")
        icon = {"PASS": "✓", "FAIL": "✗", "NOT_TESTABLE": "—"}.get(verdict, "?")
        rule_id = rf.get("rule_id", "?")
        cell_type = rf.get("cell_type", "?")
        basis = rf.get("mechanistic_basis", "")[:80]
        lines.append(f"  {icon} {cell_type} / {rule_id}: {basis}")
    return "\n".join(lines)


def format_marker_genes(marker_genes: dict[str, float], top_n: int = 20) -> str:
    """Format top N marker genes by z-score."""
    if not marker_genes:
        return "  (no marker gene data)"
    top = sorted(marker_genes.items(), key=lambda x: -x[1])[:top_n]
    return "  " + ", ".join(f"{g}({z:.1f})" for g, z in top)


def format_active_regulons(regulons: dict[str, float], thresh: float = 1.5) -> str:
    """Format regulons above the active threshold."""
    active = {k: v for k, v in regulons.items() if v >= thresh}
    if not active:
        return "  (no active regulons above z=1.5)"
    top = sorted(active.items(), key=lambda x: -x[1])[:10]
    return "  " + ", ".join(f"{r}(z={z:.1f})" for r, z in top)


def build_user_prompt(
    cluster_id: str,
    tier1_predictions: dict[str, float],
    tier2_rules_fired: list[dict],
    tier2_summary: str,
    tier3_streams_available: list[str],
    tier3_agreement: str,
    tier3_details: str,
    marker_genes: dict[str, float],
    regulons: dict[str, float],
    spatial_context: str,
    rag_context: str,
) -> str:
    """Assemble the full user prompt from structured evidence fields."""
    return USER_TEMPLATE.format(
        cluster_id=cluster_id,
        tier1_predictions=format_tier1_predictions(tier1_predictions),
        tier2_rules_fired=format_tier2_rules(tier2_rules_fired),
        tier2_summary=tier2_summary or "(none)",
        tier3_streams_available=", ".join(tier3_streams_available) or "none",
        tier3_agreement=tier3_agreement or "n/a",
        tier3_details=tier3_details or "  (no Tier 3 details)",
        marker_genes=format_marker_genes(marker_genes),
        active_regulons=format_active_regulons(regulons),
        spatial_context=spatial_context or "  (no spatial context)",
        rag_context=rag_context or "  (no RAG results)",
    )
