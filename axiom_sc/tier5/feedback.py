"""
Tier 5 → Tier 2 KG feedback loop.

Validated novel attractors auto-generate KG rule candidates for expert review.
These candidates enter the review queue (axiom_sc/kg/review_cli.py) and — after
researcher approval — are merged to the KG with a publication PMID.

This feedback loop is what distinguishes AXIOM-SC from static annotation tools:
discoveries made in Tier 5 strengthen Tier 2's mechanistic knowledge base,
turning novel cell states into first-class annotatable cell types.

License: Apache 2.0
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from axiom_sc.tier5.candidate_state import CandidateState


@dataclass
class KGCandidate:
    """
    A KG rule candidate generated from a validated Tier 5 discovery.
    Schema matches oracle_kg_v0.2.0.json (PENDING_REVIEW status).
    """

    cell_type: str
    rule_id: str
    rule_type: str      # "positive" | "circuit"
    evidence_source: str
    gene_or_regulon: list[str]
    direction: str
    paired_with: list[str]
    incompatible_with: list[str]
    mechanistic_basis: str
    pmid: str           # "NEEDS_REVIEW" until publication
    confidence: str     # "low" until validated
    tissue_context: list[str]
    source_db: str
    status: str         # always "PENDING_REVIEW"
    added_in_version: str
    candidate_id: str   # links back to the CandidateState


class DiscoveryFeedback:
    """
    Converts validated CandidateState entries into KG rule candidates.

    Usage
    -----
    feedback = DiscoveryFeedback()
    candidates = feedback.generate_candidates(validated_state)
    feedback.write_candidates(candidates, "kg_data/candidates/tier5_batch.json")
    """

    VERSION = "0.2.0"

    def generate_candidates(
        self,
        state: CandidateState,
        proposed_cell_type: str | None = None,
    ) -> list[KGCandidate]:
        """
        Generate KG rule candidates from a validated CandidateState.

        Only VALIDATED or PUBLISHED states should be passed here.
        CANDIDATE states raise a ValueError.

        Parameters
        ----------
        state : CandidateState
        proposed_cell_type : str, optional
            Cell type label to use in rules. Falls back to state.proposed_name,
            then to f"Novel_{state.candidate_id[:8]}".
        """
        if state.status == "CANDIDATE":
            raise ValueError(
                f"CandidateState {state.candidate_id} is still CANDIDATE. "
                "Validate the discovery before generating KG candidates."
            )

        cell_type = (
            proposed_cell_type
            or state.proposed_name
            or f"Novel_{state.candidate_id[:8]}"
        )

        candidates: list[KGCandidate] = []

        # 1. Positive rules — one per active TF circuit
        for i, tf in enumerate(state.active_circuits[:6], start=1):
            candidates.append(KGCandidate(
                cell_type=cell_type,
                rule_id=f"{cell_type.upper()[:8]}_POS_{i:03d}",
                rule_type="positive",
                evidence_source="regulon",
                gene_or_regulon=[tf],
                direction="active",
                paired_with=[],
                incompatible_with=[],
                mechanistic_basis=(
                    f"{tf} regulon was found active (SCENIC AUC z ≥ 1.5) in "
                    f"validated novel attractor {state.candidate_id}. "
                    f"Mechanistic role requires literature confirmation."
                ),
                pmid="NEEDS_REVIEW",
                confidence="low",
                tissue_context=[state.spatial_context] if state.spatial_context else [],
                source_db="AXIOM_Tier5",
                status="PENDING_REVIEW",
                added_in_version=self.VERSION,
                candidate_id=state.candidate_id,
            ))

        # 2. Circuit rule — top 2 TF circuits combined (if ≥ 2 active)
        if len(state.active_circuits) >= 2:
            top2 = state.active_circuits[:2]
            label_part = "_".join(tf[:4].upper() for tf in top2)
            candidates.append(KGCandidate(
                cell_type=cell_type,
                rule_id=f"{cell_type.upper()[:8]}_CIRCUIT_001",
                rule_type="circuit",
                evidence_source="regulon",
                gene_or_regulon=top2,
                direction="active",
                paired_with=[],
                incompatible_with=[],
                mechanistic_basis=(
                    f"Co-activation of {' and '.join(top2)} regulons defines the "
                    f"transcriptional identity of novel attractor {state.candidate_id}. "
                    f"Both were active in SCENIC at z ≥ 1.5. Requires experimental validation."
                ),
                pmid="NEEDS_REVIEW",
                confidence="low",
                tissue_context=[state.spatial_context] if state.spatial_context else [],
                source_db="AXIOM_Tier5",
                status="PENDING_REVIEW",
                added_in_version=self.VERSION,
                candidate_id=state.candidate_id,
            ))

        return candidates

    def write_candidates(
        self,
        candidates: list[KGCandidate],
        output_path: str | Path,
        append: bool = True,
    ) -> Path:
        """
        Write KG candidates to a JSON file for the review queue.

        Parameters
        ----------
        candidates : list[KGCandidate]
        output_path : str | Path
        append : bool
            If True and file exists, merge with existing candidates.
        """
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)

        existing: list[dict] = []
        if append and path.exists():
            with open(path) as f:
                data = json.load(f)
                existing = data.get("candidates", data) if isinstance(data, dict) else data

        new_dicts = [_candidate_to_dict(c) for c in candidates]
        # Deduplicate by rule_id
        existing_ids = {r["rule_id"] for r in existing}
        merged = existing + [d for d in new_dicts if d["rule_id"] not in existing_ids]

        payload = {
            "metadata": {
                "source": "AXIOM_Tier5",
                "generated_at": datetime.now(tz=timezone.utc).isoformat(),
                "n_candidates": len(merged),
            },
            "candidates": merged,
        }
        with open(path, "w") as f:
            json.dump(payload, f, indent=2)
        return path


def _candidate_to_dict(c: KGCandidate) -> dict:
    return {
        "cell_type": c.cell_type,
        "rule_id": c.rule_id,
        "rule_type": c.rule_type,
        "evidence_source": c.evidence_source,
        "gene_or_regulon": c.gene_or_regulon,
        "direction": c.direction,
        "paired_with": c.paired_with,
        "incompatible_with": c.incompatible_with,
        "mechanistic_basis": c.mechanistic_basis,
        "pmid": c.pmid,
        "confidence": c.confidence,
        "tissue_context": c.tissue_context,
        "source_db": c.source_db,
        "status": c.status,
        "added_in_version": c.added_in_version,
        "candidate_id": c.candidate_id,
    }
